#!/usr/bin/env python3
"""Matched long-read versus short-read HML-2 comparison.

The historical implementation reread the large short-read table several times,
treated only the literal value ``Intact`` as translatable, and counted some
structural summaries at the raw-row level.  This version scans each catalog
once, collapses copies and haplotypes to one person-by-locus cell, retains an
explicit unknown state, and bootstraps people rather than individual rows.

The short-read input is filtered during streaming to people and loci present in
the long-read catalog, keeping memory use small even for the ~1 GB source TSV.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


ORFS = ("gag", "pro", "pol", "env")
PRODUCTS = ("gag", "gag_pro", "gag_pro_pol", "env")
PRODUCT_LABELS = {
    "gag": "Gag",
    "gag_pro": "Gag-Pro",
    "gag_pro_pol": "Gag-Pro-Pol",
    "env": "Env-region ORF",
}
FEATURES = ("provirus", "multi_copy", "solo_ltr")
FEATURE_LABELS = {
    "provirus": "Any provirus",
    "multi_copy": "Multi-copy state",
    "solo_ltr": "Solo-LTR state",
}
INTACT = {"intact", "intact_fs_end", "frameshift_at_end"}
BROKEN = {"nonsense", "frameshift", "insertion", "deletion", "start_lost"}
UNKNOWN = {"", ".", "na", "n/a", "undetermined", "protein_missing"}
ALL_BITS = (1 << len(PRODUCTS)) - 1
PERSON_RE = re.compile(r"^(?:HG|NA)\d+$")


@dataclass
class Cell:
    """Compact aggregation of every row for one person and one locus."""

    hap_mask: int = 0
    recognized_structure: bool = False
    unknown_structure: bool = False
    has_provirus: bool = False
    has_multi_copy: bool = False
    has_solo_ltr: bool = False
    positive_products: int = 0
    negative_products: int = 0
    unknown_products: int = 0

    def product_state(self, index: int, require_diploid_negative: bool = True) -> int | None:
        bit = 1 << index
        if self.positive_products & bit:
            return 1
        if self.unknown_products & bit:
            return None
        if not (self.negative_products & bit):
            return None
        if require_diploid_negative and (self.hap_mask & 0b11) != 0b11:
            return None
        return 0

    def feature_state(self, feature: str, require_diploid_negative: bool = True) -> int | None:
        value = {
            "provirus": self.has_provirus,
            "multi_copy": self.has_multi_copy,
            "solo_ltr": self.has_solo_ltr,
        }[feature]
        if value:
            return 1
        if self.unknown_structure or not self.recognized_structure:
            return None
        if require_diploid_negative and (self.hap_mask & 0b11) != 0b11:
            return None
        return 0


def normalize_locus(value: str) -> str:
    locus = value.strip()
    locus = re.sub(r"HML-2_7p22\.1[ab]$", "HML-2_7p22.1", locus)
    return locus


def is_sex_chromosome_locus(locus: str) -> bool:
    short = locus.removeprefix("HML-2_")
    return short.startswith("X") or short.startswith("Y")


def haplotype_bit(value: str) -> int:
    text = value.strip().lower()
    if text in {"h1", "hap1", "mat", "maternal", "1"}:
        return 0b01
    if text in {"h2", "hap2", "pat", "paternal", "2"}:
        return 0b10
    return 0b100


def status_state(value: str | None) -> int | None:
    text = (value or "").strip().lower()
    if text in INTACT:
        return 1
    if text in BROKEN:
        return 0
    if text in UNKNOWN:
        return None
    # Fail closed for an unrecognized vocabulary item.
    return None


def product_states(row: dict[str, str]) -> list[int | None]:
    states = {orf: status_state(row.get(orf)) for orf in ORFS}

    def combine(required: Sequence[str]) -> int | None:
        values = [states[name] for name in required]
        if any(value == 0 for value in values):
            return 0
        if all(value == 1 for value in values):
            return 1
        return None

    return [
        combine(("gag",)),
        combine(("gag", "pro")),
        combine(("gag", "pro", "pol")),
        combine(("env",)),
    ]


def structure_class(value: str | None) -> str:
    text = (value or "").strip().lower()
    if "provirus_from_multi" in text or "multi" in text:
        return "multi_copy"
    if "provirus" in text:
        return "provirus"
    if "solo-ltr" in text or "solo_ltr" in text:
        return "solo_ltr"
    if text == "absent" or "explicit_absent" in text:
        return "absent"
    if "fragment" in text:
        return "fragment"
    # Insertion_Absent is a generated placeholder in the long-read catalog and
    # is not authorized as an observed empty allele.
    return "unknown"


def add_row(cell: Cell, row: dict[str, str]) -> None:
    cell.hap_mask |= haplotype_bit(row.get("Haplotype", ""))
    structure = structure_class(row.get("Structure"))
    if structure == "unknown":
        cell.unknown_structure = True
        cell.unknown_products |= ALL_BITS
        return

    cell.recognized_structure = True
    if structure in {"provirus", "multi_copy"}:
        cell.has_provirus = True
        cell.has_multi_copy |= structure == "multi_copy"
        for index, state in enumerate(product_states(row)):
            bit = 1 << index
            if state == 1:
                cell.positive_products |= bit
            elif state == 0:
                cell.negative_products |= bit
            else:
                cell.unknown_products |= bit
    else:
        cell.has_solo_ltr |= structure == "solo_ltr"
        cell.negative_products |= ALL_BITS


def iter_rows(path: Path) -> Iterable[dict[str, str]]:
    with (gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"Input has no header: {path}")
        required = {"Locus", "ID", "Haplotype", "Structure", *ORFS}
        missing = required - set(reader.fieldnames)
        if missing:
            raise RuntimeError(f"{path} lacks columns: {sorted(missing)}")
        yield from reader


def load_catalog(
    path: Path,
    allowed_people: set[str] | None = None,
    allowed_loci: set[str] | None = None,
) -> tuple[dict[tuple[str, str], Cell], set[str], set[str], int]:
    cells: dict[tuple[str, str], Cell] = {}
    people: set[str] = set()
    loci: set[str] = set()
    retained_rows = 0
    for row in iter_rows(path):
        if "analysis_include" in row and row["analysis_include"] != "1":
            continue
        person = (row.get("ID") or "").strip()
        if not PERSON_RE.fullmatch(person):
            continue
        locus = normalize_locus(row.get("Locus") or "")
        if not locus:
            continue
        if allowed_people is not None and person not in allowed_people:
            continue
        if allowed_loci is not None and locus not in allowed_loci:
            continue
        people.add(person)
        loci.add(locus)
        cell = cells.setdefault((person, locus), Cell())
        add_row(cell, row)
        retained_rows += 1
    return cells, people, loci, retained_rows


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def two_way_bootstrap_ratio(
    numerators: np.ndarray,
    denominators: np.ndarray,
    boot: int,
    seed: int,
) -> tuple[float, float, float, str]:
    """Ratio with a two-way person/locus bootstrap.

    At an observed 0% or 100% boundary, resampling cannot create unseen events;
    a Wilson interval is therefore reported as a finite descriptive bound.
    """
    denominator = float(denominators.sum())
    estimate = float(numerators.sum() / denominator) if denominator else math.nan
    if not denominator:
        return estimate, math.nan, math.nan, "not_estimable"
    successes = int(numerators.sum())
    total = int(denominator)
    if successes in {0, total}:
        low, high = wilson(successes, total)
        return estimate, low, high, "Wilson_boundary"
    if boot <= 0:
        return estimate, math.nan, math.nan, "none"
    rng = np.random.default_rng(seed)
    values = np.empty(boot, dtype=float)
    n_people, n_loci = numerators.shape
    for index in range(boot):
        sampled_people = rng.integers(0, n_people, size=n_people)
        sampled_loci = rng.integers(0, n_loci, size=n_loci)
        sampled_num = numerators[np.ix_(sampled_people, sampled_loci)]
        sampled_den = denominators[np.ix_(sampled_people, sampled_loci)]
        den = sampled_den.sum()
        values[index] = sampled_num.sum() / den if den else math.nan
    low, high = np.nanpercentile(values, (2.5, 97.5))
    return estimate, float(low), float(high), "person_locus_bootstrap"


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def write_tsv(path: Path, rows: Sequence[dict], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def analyze(
    lr_cells: dict[tuple[str, str], Cell],
    sr_cells: dict[tuple[str, str], Cell],
    people: Sequence[str],
    loci: Sequence[str],
    boot: int,
    seed: int,
) -> tuple[list[dict], list[dict], list[dict]]:
    person_index = {person: index for index, person in enumerate(people)}
    locus_index = {locus: index for index, locus in enumerate(loci)}
    metrics: list[dict] = []
    per_locus: list[dict] = []
    unknowns: list[dict] = []

    def append_metric(
        domain: str,
        feature: str,
        label: str,
        metric: str,
        numerator: np.ndarray,
        denominator: np.ndarray,
        metric_seed: int,
    ) -> None:
        estimate, low, high, method = two_way_bootstrap_ratio(
            numerator, denominator, boot, metric_seed
        )
        metrics.append({
            "domain": domain,
            "feature": feature,
            "label": label,
            "metric": metric,
            "estimate": estimate,
            "ci_low": low,
            "ci_high": high,
            "numerator": int(numerator.sum()),
            "denominator": int(denominator.sum()),
            "ci_method": method,
            "bootstrap_unit": "person_and_locus" if "bootstrap" in method else "fixed_panel",
        })

    for product_index, product in enumerate(PRODUCTS):
        shape = (len(people), len(loci))
        recovered = np.zeros(shape, dtype=np.int8)
        sr_callable = np.zeros(shape, dtype=np.int8)
        lr_positive = np.zeros(shape, dtype=np.int8)
        discordant = np.zeros(shape, dtype=np.int8)
        lr_negative_sr_callable = np.zeros(shape, dtype=np.int8)
        locus_counts = {locus: [0, 0, 0] for locus in loci}  # recovered, LR+, SR unknown
        for person in people:
            pidx = person_index[person]
            for locus in loci:
                lidx = locus_index[locus]
                require_diploid = not is_sex_chromosome_locus(locus)
                lr = lr_cells.get((person, locus))
                sr = sr_cells.get((person, locus))
                lr_state = lr.product_state(product_index, require_diploid) if lr else None
                sr_state = sr.product_state(product_index, require_diploid) if sr else None
                if lr_state == 1:
                    lr_positive[pidx, lidx] = 1
                    locus_counts[locus][1] += 1
                    if sr_state is None:
                        locus_counts[locus][2] += 1
                    else:
                        sr_callable[pidx, lidx] = 1
                        if sr_state == 1:
                            recovered[pidx, lidx] = 1
                            locus_counts[locus][0] += 1
                elif lr_state == 0 and sr_state is not None:
                    lr_negative_sr_callable[pidx, lidx] = 1
                    if sr_state == 1:
                        discordant[pidx, lidx] = 1
        append_metric(
            "coding_product", product, PRODUCT_LABELS[product],
            "short_read_callability_among_long_read_positive",
            sr_callable, lr_positive, seed + product_index,
        )
        append_metric(
            "coding_product", product, PRODUCT_LABELS[product],
            "conditional_positive_concordance",
            recovered, sr_callable, seed + 20 + product_index,
        )
        append_metric(
            "coding_product", product, PRODUCT_LABELS[product],
            "overall_workflow_recovery",
            recovered, lr_positive, seed + 40 + product_index,
        )
        append_metric(
            "coding_product", product, PRODUCT_LABELS[product],
            "long_read_negative_short_read_positive_discordance",
            discordant, lr_negative_sr_callable, seed + 60 + product_index,
        )
        for locus, (recovered, lr_positive, sr_unknown) in locus_counts.items():
            callable_total = lr_positive - sr_unknown
            low_locus, high_locus = wilson(recovered, callable_total)
            per_locus.append({
                "feature": product,
                "locus": locus,
                "lr_positive_people": lr_positive,
                "sr_callable_lr_positive_people": callable_total,
                "sr_recovered_people": recovered,
                "sr_unknown_people": sr_unknown,
                "conditional_positive_concordance": recovered / callable_total if callable_total else math.nan,
                "overall_workflow_recovery": recovered / lr_positive if lr_positive else math.nan,
                "ci_low": low_locus,
                "ci_high": high_locus,
                "missed_callable_people": callable_total - recovered,
            })

    for feature_index, feature in enumerate(FEATURES):
        shape = (len(people), len(loci))
        recovered = np.zeros(shape, dtype=np.int8)
        sr_callable = np.zeros(shape, dtype=np.int8)
        lr_positive = np.zeros(shape, dtype=np.int8)
        discordant = np.zeros(shape, dtype=np.int8)
        lr_negative_sr_callable = np.zeros(shape, dtype=np.int8)
        for person in people:
            pidx = person_index[person]
            for locus in loci:
                lidx = locus_index[locus]
                require_diploid = not is_sex_chromosome_locus(locus)
                lr = lr_cells.get((person, locus))
                sr = sr_cells.get((person, locus))
                lr_state = lr.feature_state(feature, require_diploid) if lr else None
                sr_state = sr.feature_state(feature, require_diploid) if sr else None
                if lr_state == 1:
                    lr_positive[pidx, lidx] = 1
                    if sr_state is not None:
                        sr_callable[pidx, lidx] = 1
                        if sr_state == 1:
                            recovered[pidx, lidx] = 1
                elif lr_state == 0 and sr_state is not None:
                    lr_negative_sr_callable[pidx, lidx] = 1
                    if sr_state == 1:
                        discordant[pidx, lidx] = 1
        append_metric(
            "carrier_structural_state", feature, FEATURE_LABELS[feature],
            "short_read_callability_among_long_read_positive",
            sr_callable, lr_positive, seed + 100 + feature_index,
        )
        append_metric(
            "carrier_structural_state", feature, FEATURE_LABELS[feature],
            "conditional_positive_concordance",
            recovered, sr_callable, seed + 120 + feature_index,
        )
        append_metric(
            "carrier_structural_state", feature, FEATURE_LABELS[feature],
            "overall_workflow_recovery",
            recovered, lr_positive, seed + 140 + feature_index,
        )
        append_metric(
            "carrier_structural_state", feature, FEATURE_LABELS[feature],
            "long_read_negative_short_read_positive_discordance",
            discordant, lr_negative_sr_callable, seed + 160 + feature_index,
        )

    for person in people:
        for locus in loci:
            lr = lr_cells.get((person, locus))
            sr = sr_cells.get((person, locus))
            if lr is None or sr is None:
                unknowns.append({
                    "person": person,
                    "locus": locus,
                    "long_read_cell": "present" if lr else "missing",
                    "short_read_cell": "present" if sr else "missing",
                })
    return metrics, per_locus, unknowns


def make_figure(metrics: Sequence[dict], per_locus: Sequence[dict], output_base: Path) -> None:
    """Optional aggregate plot; the manuscript reports these counts in Table 2."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coding = [row for row in metrics if row["domain"] == "coding_product"
              and row["metric"] == "overall_workflow_recovery"]
    structural = [row for row in metrics if row["domain"] == "carrier_structural_state"
                  and row["metric"] == "overall_workflow_recovery"]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.0))
    for ax, rows, title, letter in zip(
        axes, (coding, structural), ("Coding sequences", "Structural states"), ("A", "B")
    ):
        x = np.arange(len(rows))
        y = np.array([float(row["estimate"]) for row in rows]) * 100
        ax.vlines(x, 0, y, color="#333333", linewidth=1.5)
        ax.scatter(x, y, s=40, color="#333333", zorder=3)
        ax.set_xticks(x, [row["label"] for row in rows], fontsize=8)
        ax.set_ylim(0, 118)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_ylabel("Recovered by short reads (%)")
        ax.set_title(title, loc="left", fontsize=11)
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.16, 1.07, letter, transform=ax.transAxes, fontsize=14, fontweight="bold")
        for i, row in enumerate(rows):
            ax.text(i, y[i] + 3, f"{y[i]:.1f}%\n{int(row['numerator']):,}/{int(row['denominator']):,}",
                    ha="center", va="bottom", fontsize=8)
    fig.subplots_adjust(top=0.87, left=0.08, right=0.99, bottom=0.16, wspace=0.32)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(output_base.with_suffix(f".{extension}"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--long-read", type=Path)
    parser.add_argument("--short-read", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--skip-hashes", action="store_true")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--plot-only", action="store_true",
                        help="Rebuild the figure from existing compact result TSVs")
    args = parser.parse_args()

    if args.plot_only:
        def read_result(path: Path) -> list[dict]:
            with (gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")) as handle:
                return list(csv.DictReader(handle, delimiter="\t"))

        metrics = read_result(args.output_dir / "results" / "summary_metrics.tsv")
        per_locus = read_result(args.output_dir / "results" / "per_locus_coding_recovery.tsv")
        for row in metrics:
            for key in ("estimate", "ci_low", "ci_high"):
                row[key] = float(row[key])
            for key in ("numerator", "denominator"):
                row[key] = int(row[key])
        for row in per_locus:
            for key in (
                "conditional_positive_concordance", "overall_workflow_recovery",
                "ci_low", "ci_high",
            ):
                row[key] = float(row[key])
            for key in (
                "lr_positive_people", "sr_callable_lr_positive_people",
                "sr_recovered_people", "sr_unknown_people", "missed_callable_people",
            ):
                row[key] = int(row[key])
        make_figure(
            metrics, per_locus,
            args.output_dir / "figures" / "hml2_longread_shortread_validated_comparison",
        )
        print(f"Rebuilt figure from compact results in {args.output_dir}")
        return

    if args.long_read is None or args.short_read is None:
        parser.error("--long-read and --short-read are required unless --plot-only is used")

    print("Loading long-read catalog once...")
    lr_cells, lr_people, lr_loci, lr_rows = load_catalog(args.long_read)
    print(f"  {len(lr_people)} people; {len(lr_loci)} loci; {lr_rows:,} retained rows")
    print("Streaming short-read catalog once, retaining only long-read candidates...")
    sr_cells, sr_people, sr_loci, sr_rows = load_catalog(
        args.short_read, allowed_people=lr_people, allowed_loci=lr_loci
    )
    people = sorted(lr_people & sr_people)
    loci = sorted(lr_loci & sr_loci)
    if not people or not loci:
        raise RuntimeError("No shared people or loci between inputs")
    long_read_only_loci = sorted(lr_loci - sr_loci)
    print(f"Matched analysis: {len(people)} people x {len(loci)} loci; {sr_rows:,} retained short-read rows")
    print(f"Long-read loci absent from the legacy short-read universe: {len(long_read_only_loci)}")

    metrics, per_locus, unknowns = analyze(
        lr_cells, sr_cells, people, loci, args.bootstrap, args.seed
    )
    results = args.output_dir / "results"
    derived = args.output_dir / "derived"
    figures = args.output_dir / "figures"
    metric_columns = (
        "domain", "feature", "label", "metric", "estimate", "ci_low", "ci_high",
        "numerator", "denominator", "ci_method", "bootstrap_unit",
    )
    locus_columns = (
        "feature", "locus", "lr_positive_people", "sr_callable_lr_positive_people",
        "sr_recovered_people", "sr_unknown_people", "conditional_positive_concordance",
        "overall_workflow_recovery", "ci_low", "ci_high", "missed_callable_people",
    )
    write_tsv(results / "summary_metrics.tsv", metrics, metric_columns)
    write_tsv(results / "per_locus_coding_recovery.tsv", per_locus, locus_columns)
    write_tsv(
        derived / "missing_matched_cells.tsv", unknowns,
        ("person", "locus", "long_read_cell", "short_read_cell"),
    )
    write_tsv(
        derived / "locus_universe.tsv",
        [
            {
                "locus": locus,
                "long_read_catalog": "yes",
                "legacy_short_read_catalog": "yes" if locus in sr_loci else "no",
                "analysis_role": "matched" if locus in sr_loci else "long_read_only_not_represented",
            }
            for locus in sorted(lr_loci)
        ],
        ("locus", "long_read_catalog", "legacy_short_read_catalog", "analysis_role"),
    )
    provenance = []
    for role, path, retained in (
        ("long_read_catalog", args.long_read, lr_rows),
        ("short_read_catalog", args.short_read, sr_rows),
    ):
        provenance.append({
            "role": role,
            "path_basename": path.name,
            "bytes": path.stat().st_size,
            "sha256": "SKIPPED" if args.skip_hashes else sha256_file(path),
            "retained_rows": retained,
            "shared_people": len(people),
            "shared_loci": len(loci),
        })
    write_tsv(
        derived / "input_provenance.tsv", provenance,
        ("role", "path_basename", "bytes", "sha256", "retained_rows", "shared_people", "shared_loci"),
    )
    if not args.skip_plots:
        make_figure(metrics, per_locus, figures / "hml2_longread_shortread_validated_comparison")
    print(f"Wrote validated comparison to {args.output_dir}")


if __name__ == "__main__":
    main()
