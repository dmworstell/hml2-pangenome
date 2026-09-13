#!/usr/bin/env python3
"""Validate refreshed TSD calls and rebuild Supplementary Figure S11."""

from __future__ import annotations

import csv
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from figure_style import BLUE, GRAY, INK, ORANGE, SKY, apply_style, finish_axis, panel_title, save_figure


ROOT = Path(__file__).resolve().parents[1]
NEW_TABLE = ROOT / "inputs/orf_analysis/tsd_refresh_20260803/combined_hml2_orf_analysis.tsv"
PRIOR_TABLE = ROOT / "inputs/orf_analysis/current_v3r1/combined_hml2_orf_analysis.CNV_WEIGHTED.v3r1.tsv"
ANNOTATION_TABLE = (
    ROOT
    / "results/biological_orf_annotation_20260802/"
    "combined_hml2_orf_analysis.CNV_WEIGHTED.BIOLOGICALLY_ANNOTATED.v3.tsv"
)
DIRECT_EVIDENCE = (
    ROOT
    / "inputs/orf_analysis/tsd_refresh_20260803/review/"
    "candidate_tsd_direct_evidence_join.tsv"
)
SAMPLE_METADATA = ROOT / "inputs/sampling_frame/human_sample_frame.tsv"
OUT = ROOT / "manuscript/tsd_analysis"
FIGURE_OUT = ROOT / "manuscript/figures/tsd_refresh"

UNAVAILABLE = {"", "NA", "N/A", "NONE", "UNKNOWN", "UNOBSERVED", "-", "ERROR", "MAP_ERROR"}
PRIMARY_3Q_PAIR = ("GAGGT", "GAGGT")
ALTERNATE_3Q_PAIR = ("GAGGT", "GAGAT")


def normalized(value: str | None) -> str:
    clean = str(value or "").strip().upper()
    return "NONE" if clean in UNAVAILABLE else clean


def display_locus(value: str) -> str:
    if value == "HML-2_acro_type1":
        return "Acrocentric Type I"
    if value == "HML-2_acro_type2":
        return "Acrocentric Type II"
    return value.removeprefix("HML-2_").removesuffix("_hg38").removesuffix("_new")


def admissible_pair(left: str, right: str) -> bool:
    if left == "NONE" or right == "NONE":
        return False
    if len(left) != len(right) or len(left) not in {4, 5, 6}:
        return False
    if not re.fullmatch(r"[ACGT]+", left + right):
        return False
    return sum(a != b for a, b in zip(left, right)) <= 2


def read_annotation() -> dict[tuple[str, str], dict[str, str]]:
    annotations: dict[tuple[str, str], dict[str, str]] = {}
    with ANNOTATION_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["Locus"], row["ID_Full"])
            if key in annotations:
                raise ValueError(f"duplicate biological annotation key: {key}")
            annotations[key] = {
                "analysis_include": row["analysis_include"],
                "analysis_exclusion_reason": row["analysis_exclusion_reason"],
                "physical_locus_assignment": row["physical_locus_assignment"],
            }
    return annotations


def read_prior_pairs() -> dict[tuple[str, str], tuple[str, str]]:
    pairs: dict[tuple[str, str], tuple[str, str]] = {}
    with PRIOR_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["Locus"], row["ID_Full"])
            if key in pairs:
                raise ValueError(f"duplicate prior TSD key: {key}")
            pairs[key] = (normalized(row["5'_TSD"]), normalized(row["3'_TSD"]))
    return pairs


def read_superpopulations() -> dict[str, str]:
    result: dict[str, str] = {}
    with SAMPLE_METADATA.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            superpopulation = row["superpopulation"].strip()
            if superpopulation in {"AFR", "AMR", "EAS", "EUR", "SAS"}:
                result[row["sample_id"]] = superpopulation
    return result


def read_direct_evidence() -> tuple[dict[tuple[str, str], dict[str, str]], Counter]:
    rows: dict[tuple[str, str], dict[str, str]] = {}
    counts: Counter = Counter()
    with DIRECT_EVIDENCE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["candidate_locus"], row["ID_Full"])
            if key in rows:
                raise ValueError(f"duplicate direct-evidence key: {key}")
            rows[key] = row
            supported = row["review_flag"] == "NONE"
            counts["directly_supported" if supported else "flagged"] += 1
            if row["changed_or_new"] == "YES":
                counts["changed_directly_supported" if supported else "changed_flagged"] += 1
    return rows, counts


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    annotations = read_annotation()
    prior_pairs = read_prior_pairs()
    superpopulations = read_superpopulations()
    direct_evidence, direct_counts = read_direct_evidence()

    counts: Counter = Counter()
    loci: set[str] = set()
    seen: set[tuple[str, str]] = set()
    per_locus: dict[str, Counter] = defaultdict(Counter)
    per_superpopulation: dict[str, Counter] = defaultdict(Counter)
    directly_supported_mismatch: dict[str, Counter] = defaultdict(Counter)
    accepted_lengths: Counter = Counter()
    changed_rows: list[dict[str, object]] = []
    accepted_rows: list[dict[str, str]] = []
    direct_spot_rows: list[dict[str, object]] = []
    unmatched_annotation_keys: list[tuple[str, str]] = []

    with NEW_TABLE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Locus", "ID_Full", "ID", "Structure", "5'_TSD", "3'_TSD",
            "present_day_tsd_mismatch_count", "tsd_observable_state",
            "tsd_call_status", "tsd_pair_admission", "TSD_Orientation",
        }
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"refreshed TSD table lacks columns: {missing}")

        for row in reader:
            key = (row["Locus"], row["ID_Full"])
            if key in seen:
                raise ValueError(f"duplicate refreshed TSD key: {key}")
            seen.add(key)
            counts["rows"] += 1
            loci.add(row["Locus"])
            if row["Locus"] == "HML-2_8q24.3b":
                raise ValueError("retired 8q24.3b alias remains in refreshed TSD table")

            annotation = annotations.get(key)
            if annotation is None:
                unmatched_annotation_keys.append(key)
                include = True
            else:
                include = annotation["analysis_include"].strip().lower() in {
                    "1", "true", "yes"
                }
                if not include:
                    counts["excluded:" + annotation["analysis_exclusion_reason"]] += 1
            if not include:
                continue
            counts["analysis_rows"] += 1

            left = normalized(row["5'_TSD"])
            right = normalized(row["3'_TSD"])
            status = row["tsd_call_status"]
            admission = row["tsd_pair_admission"]
            counts["status:" + status] += 1
            counts["admission:" + admission] += 1
            per_locus[row["Locus"]]["analysis_rows"] += 1
            per_locus[row["Locus"]]["status:" + status] += 1

            prior_left, prior_right = prior_pairs.get(key, ("<ABSENT>", "<ABSENT>"))
            prior_admissible = admissible_pair(prior_left, prior_right)
            if prior_admissible:
                counts["prior_admissible_pairs"] += 1
            if (left, right) != (prior_left, prior_right):
                counts["changed_pair"] += 1
            if admission != "PAIRED_ACCEPTED":
                if prior_admissible:
                    counts["prior_admissible_not_retained"] += 1
                continue

            if left == "NONE" or right == "NONE":
                raise ValueError(f"accepted pair is unavailable: {key}")
            if len(left) != len(right) or len(left) not in {4, 5, 6}:
                raise ValueError(f"accepted pair has invalid lengths: {key} {left}/{right}")
            if not re.fullmatch(r"[ACGT]+", left + right):
                raise ValueError(f"accepted pair contains non-ACGT characters: {key}")
            mismatch = sum(a != b for a, b in zip(left, right))
            if mismatch > 2:
                raise ValueError(f"accepted pair exceeds two mismatches: {key}")
            if str(mismatch) != row["present_day_tsd_mismatch_count"].strip():
                raise ValueError(f"recorded mismatch count disagrees: {key}")

            accepted_lengths[len(left)] += 1
            counts["accepted_pairs"] += 1
            if prior_admissible:
                counts["accepted_pair_also_admissible_before_refresh"] += 1
            else:
                counts["newly_admissible_pair_after_refresh"] += 1
            counts[f"accepted_mismatch:{mismatch}"] += 1
            per_locus[row["Locus"]]["accepted_pairs"] += 1
            per_locus[row["Locus"]][f"mismatch:{mismatch}"] += 1
            if mismatch:
                per_locus[row["Locus"]]["different_pairs"] += 1

            direct = direct_evidence.get(key)
            if direct is not None and direct["review_flag"] == "NONE":
                observed = (
                    normalized(direct["evidence_5prime_tsd"]),
                    normalized(direct["evidence_3prime_tsd"]),
                )
                if observed != (left, right):
                    raise ValueError(f"direct evidence disagrees with accepted pair: {key}")
                counts["accepted_directly_supported"] += 1
                directly_supported_mismatch[row["Locus"]]["accepted_pairs"] += 1
                directly_supported_mismatch[row["Locus"]][f"mismatch:{mismatch}"] += 1
                if mismatch:
                    directly_supported_mismatch[row["Locus"]]["different_pairs"] += 1
                direct_spot_rows.append({
                    "locus": row["Locus"],
                    "ID_Full": row["ID_Full"],
                    "structure": row["Structure"],
                    "prior_5prime_tsd": prior_left,
                    "prior_3prime_tsd": prior_right,
                    "refreshed_5prime_tsd": left,
                    "refreshed_3prime_tsd": right,
                    "mismatches": mismatch,
                    "changed_or_new": direct["changed_or_new"],
                    "evidence_match": direct["evidence_match"],
                    "evidence_5prime_tsd": normalized(direct["evidence_5prime_tsd"]),
                    "evidence_3prime_tsd": normalized(direct["evidence_3prime_tsd"]),
                    "left_flank_exact": direct["left_flank_exact"],
                    "right_flank_exact": direct["right_flank_exact"],
                    "boundary_owner": direct["boundary_owner"],
                })

            if prior_left == "NONE" or prior_right == "NONE" or key not in prior_pairs:
                counts["newly_recovered_accepted_pair"] += 1

            if (
                row["Locus"] == "HML-2_3q12.3"
                and direct is not None
                and direct["review_flag"] == "NONE"
            ):
                sample = row["ID"]
                superpopulation = superpopulations.get(sample)
                if superpopulation:
                    pair_class = (
                        "GAGGT/GAGGT"
                        if (left, right) == PRIMARY_3Q_PAIR
                        else "GAGGT/GAGAT"
                        if (left, right) == ALTERNATE_3Q_PAIR
                        else "other"
                    )
                    per_superpopulation[superpopulation][pair_class] += 1

            accepted_rows.append(row)
            if (left, right) != (prior_left, prior_right):
                changed_rows.append({
                    "locus": row["Locus"],
                    "ID_Full": row["ID_Full"],
                    "structure": row["Structure"],
                    "prior_5prime_tsd": prior_left,
                    "prior_3prime_tsd": prior_right,
                    "refreshed_5prime_tsd": left,
                    "refreshed_3prime_tsd": right,
                    "mismatches": mismatch,
                    "direct_evidence": "YES" if direct and direct["review_flag"] == "NONE" else "NO",
                })

    if counts["rows"] != 63750 or len(loci) != 103:
        raise ValueError(f"unexpected refreshed domain: rows={counts['rows']} loci={len(loci)}")
    if "HML-2_8q24.3c" not in loci:
        raise ValueError("8q24.3c is absent from refreshed TSD table")

    locus_rows: list[dict[str, object]] = []
    for locus, locus_counts in sorted(per_locus.items()):
        paired = locus_counts["accepted_pairs"]
        different = locus_counts["different_pairs"]
        supported = directly_supported_mismatch[locus]["accepted_pairs"]
        supported_different = directly_supported_mismatch[locus]["different_pairs"]
        locus_rows.append({
            "locus": locus,
            "display_locus": display_locus(locus),
            "analysis_rows": locus_counts["analysis_rows"],
            "paired_accepted": paired,
            "paired_with_different_copies": different,
            "fraction_different": different / paired if paired else "",
            "directly_supported_pairs": supported,
            "directly_supported_different_pairs": supported_different,
            "directly_supported_fraction_different": (
                supported_different / supported if supported else ""
            ),
        })

    write_tsv(
        OUT / "per_locus_tsd_summary.tsv",
        locus_rows,
        [
            "locus", "display_locus", "analysis_rows", "paired_accepted",
            "paired_with_different_copies", "fraction_different",
            "directly_supported_pairs", "directly_supported_different_pairs",
            "directly_supported_fraction_different",
        ],
    )

    superpopulation_rows: list[dict[str, object]] = []
    for superpopulation in ("AFR", "AMR", "EAS", "EUR", "SAS"):
        row_counts = per_superpopulation[superpopulation]
        total = sum(row_counts.values())
        for pair_class in ("GAGGT/GAGGT", "GAGGT/GAGAT", "other"):
            count = row_counts[pair_class]
            superpopulation_rows.append({
                "superpopulation": superpopulation,
                "pair_class": pair_class,
                "count": count,
                "total": total,
                "fraction": count / total if total else "",
            })
    write_tsv(
        OUT / "threeq12_3_tsd_pairs_by_superpopulation.tsv",
        superpopulation_rows,
        ["superpopulation", "pair_class", "count", "total", "fraction"],
    )

    # Deterministic spot checks cover unchanged and newly recovered calls, equal and
    # mismatched pairs, multiple structures, and the manually anchored GRCh38 call.
    grouped_spots: list[dict[str, object]] = []
    used: set[tuple[str, str]] = set()
    for target_locus in (
        "HML-2_3q21.2", "HML-2_4q32.3", "HML-2_5p12",
        "HML-2_3q12.3", "HML-2_8q24.3c", "HML-2_7p22.1",
    ):
        for row in direct_spot_rows:
            key = (str(row["locus"]), str(row["ID_Full"]))
            if row["locus"] == target_locus and key not in used:
                grouped_spots.append(row)
                used.add(key)
                break
    for mismatch in (0, 1, 2):
        for row in direct_spot_rows:
            key = (str(row["locus"]), str(row["ID_Full"]))
            if row["mismatches"] == mismatch and key not in used:
                grouped_spots.append(row)
                used.add(key)
                break
    for structure in ("Provirus", "Solo-LTR", "Tandem"):
        for row in direct_spot_rows:
            key = (str(row["locus"]), str(row["ID_Full"]))
            if structure.lower() in str(row["structure"]).lower() and key not in used:
                grouped_spots.append(row)
                used.add(key)
                break
    grouped_spots = grouped_spots[:10]
    write_tsv(
        OUT / "spot_checked_recovered_tsd_calls.tsv",
        grouped_spots,
        [
            "locus", "ID_Full", "structure", "prior_5prime_tsd", "prior_3prime_tsd",
            "refreshed_5prime_tsd", "refreshed_3prime_tsd", "mismatches",
            "changed_or_new", "evidence_match", "evidence_5prime_tsd",
            "evidence_3prime_tsd", "left_flank_exact", "right_flank_exact",
            "boundary_owner",
        ],
    )

    plot_rows = [
        row
        for row in locus_rows
        if int(row["directly_supported_pairs"]) >= 20
        and int(row["directly_supported_different_pairs"]) > 0
    ]
    plot_rows.sort(
        key=lambda row: (
            float(row["directly_supported_fraction_different"]),
            str(row["display_locus"]),
        )
    )

    apply_style()
    fig_height = max(3.6, 0.32 * len(plot_rows) + 1.2)
    fig, ax = plt.subplots(figsize=(7.1, fig_height))
    labels = [str(row["display_locus"]) for row in plot_rows]
    fractions = [float(row["directly_supported_fraction_different"]) for row in plot_rows]
    y = range(len(plot_rows))
    ax.barh(y, fractions, color=ORANGE, height=0.66)
    ax.set_yticks(list(y), labels)
    ax.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.set_xlabel("directly supported TSD pairs with different 5′ and 3′ copies")
    panel_title(ax, "A", "Paired TSD differences by locus")
    for index, row in enumerate(plot_rows):
        paired = int(row["directly_supported_pairs"])
        different = int(row["directly_supported_different_pairs"])
        ax.text(
            float(row["directly_supported_fraction_different"])
            + max(fractions, default=0.01) * 0.018,
            index,
            f"{different}/{paired}",
            va="center",
            fontsize=7.4,
            color=INK,
        )
    ax.set_xlim(0, min(1.05, max(0.05, max(fractions, default=0.01) * 1.24)))
    finish_axis(ax, grid="x")
    fig.tight_layout()
    save_figure(
        fig,
        FIGURE_OUT / "Figure_S11A_tsd_differences.png",
        FIGURE_OUT / "Figure_S11A_tsd_differences.pdf",
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.1, 4.2))
    order = ("AFR", "AMR", "EAS", "EUR", "SAS")
    categories = (
        ("GAGGT/GAGGT", GRAY),
        ("GAGGT/GAGAT", ORANGE),
        ("other", SKY),
    )
    bottoms = [0.0] * len(order)
    for category, color in categories:
        values = []
        for superpopulation in order:
            row_counts = per_superpopulation[superpopulation]
            total = sum(row_counts.values())
            values.append(row_counts[category] / total if total else 0.0)
        ax.bar(order, values, bottom=bottoms, color=color, width=0.72, label=category)
        bottoms = [left + right for left, right in zip(bottoms, values)]
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.set_ylabel("fraction of directly supported 3q12.3 TSD pairs")
    ax.set_xlabel("superpopulation")
    panel_title(ax, "B", "3q12.3 TSD pairs by superpopulation")
    ax.legend(title="5′ / 3′ TSD", bbox_to_anchor=(1.02, 1), loc="upper left")
    finish_axis(ax, grid="y")
    fig.tight_layout()
    save_figure(
        fig,
        FIGURE_OUT / "Figure_S11B_3q12_3_superpopulation.png",
        FIGURE_OUT / "Figure_S11B_3q12_3_superpopulation.pdf",
    )
    plt.close(fig)

    summary = {
        "schema": "hml2.tsd-analysis.v1",
        "source_rows": counts["rows"],
        "source_locus_labels": len(loci),
        "analysis_rows": counts["analysis_rows"],
        "accepted_pairs": counts["accepted_pairs"],
        "prior_admissible_pairs_in_current_analysis_rows": counts["prior_admissible_pairs"],
        "accepted_pairs_also_admissible_before_refresh": counts[
            "accepted_pair_also_admissible_before_refresh"
        ],
        "newly_admissible_pairs_after_refresh": counts["newly_admissible_pair_after_refresh"],
        "prior_admissible_pairs_not_retained": counts["prior_admissible_not_retained"],
        "accepted_pairs_by_length": dict(sorted(accepted_lengths.items())),
        "accepted_pairs_by_mismatch_count": {
            str(value): counts[f"accepted_mismatch:{value}"] for value in (0, 1, 2)
        },
        "newly_recovered_accepted_pairs": counts["newly_recovered_accepted_pair"],
        "accepted_pairs_with_direct_junction_support": counts["accepted_directly_supported"],
        "direct_evidence_review": dict(sorted(direct_counts.items())),
        "changed_pair_rows": counts["changed_pair"],
        "annotation_keys_not_matched": len(unmatched_annotation_keys),
        "analysis_exclusions": {
            key.removeprefix("excluded:"): value
            for key, value in sorted(counts.items())
            if key.startswith("excluded:")
        },
        "threeq12_3_superpopulation_totals": {
            superpopulation: sum(per_superpopulation[superpopulation].values())
            for superpopulation in ("AFR", "AMR", "EAS", "EUR", "SAS")
        },
        "plotted_loci": len(plot_rows),
        "spot_checks": len(grouped_spots),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "tsd_analysis_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    for source_name, manuscript_name in (
        ("per_locus_tsd_summary.tsv", "Figure_S11A_source_data.tsv"),
        ("threeq12_3_tsd_pairs_by_superpopulation.tsv", "Figure_S11B_source_data.tsv"),
        ("spot_checked_recovered_tsd_calls.tsv", "Figure_S11_spot_checks.tsv"),
        ("tsd_analysis_summary.json", "Figure_S11_analysis_summary.json"),
    ):
        shutil.copy2(OUT / source_name, ROOT / "manuscript" / manuscript_name)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
