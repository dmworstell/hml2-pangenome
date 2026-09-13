#!/usr/bin/env python3
"""Read-only 7p22.1 coordinate, copy-number, and proxy audit.

Writes only beside this script.  The production phase-1 pipeline is not used or
modified.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
CATALOG = ROOT / "HML2_ProjectResources/data/catalog/combined_hml2_orf_analysis.tsv"
TRUTH = ROOT / "manuscript_figures/python/analysis/catalog_functional_screen/results/sample_contrast_truth.tsv"
PREDICTOR = ROOT / "manuscript_figures/python/analysis/sevenp22_multicopy_tag_predictor"
VCF = ROOT / "HML2_ProjectResources/data/onekg/windows/HML-2_7p22.1.snps.vcf.gz"
LOCUS = "HML-2_7p22.1"
PUBLIC = re.compile(r"(?:HG|NA)\d+")
PART = re.compile(r"_MULTI_part(\d+)$")

# Catalog split coordinates and the authoritative complete-array panel interval.
OLD_PARTS = {
    "part1": (4_590_936, 4_600_393),
    "part2": (4_582_432, 4_591_890),
}
ARRAY = (4_606_788, 4_624_749)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    old_union = (min(x[0] for x in OLD_PARTS.values()), max(x[1] for x in OLD_PARTS.values()))
    shift = ARRAY[0] - old_union[0]
    assert ARRAY[1] - old_union[1] == shift == 24_356
    assert old_union[1] - old_union[0] + 1 == ARRAY[1] - ARRAY[0] + 1 == 17_962
    assert OLD_PARTS["part2"][1] - OLD_PARTS["part1"][0] + 1 == 955

    coordinate_rows: list[dict[str, object]] = []
    for label, (start, end) in sorted(OLD_PARTS.items()):
        coordinate_rows.append({
            "record": label,
            "catalog_interval_hg38": f"chr7:{start}-{end}(-)",
            "resolved_interval_hg38": f"chr7:{start + shift}-{end + shift}(-)",
            "constant_offset_bp": shift,
            "role": "tiled array unit; not a competing locus",
        })
    coordinate_rows.append({
        "record": "array_union",
        "catalog_interval_hg38": f"chr7:{old_union[0]}-{old_union[1]}(-)",
        "resolved_interval_hg38": f"chr7:{ARRAY[0]}-{ARRAY[1]}(-)",
        "constant_offset_bp": shift,
        "role": "authoritative complete-array anchor",
    })
    write_tsv(OUT / "coordinate_resolution.tsv", coordinate_rows)

    cells: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    all_haps: set[tuple[str, str]] = set()
    multi_source_rows = 0
    with CATALOG.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample, hap = row["ID"].strip(), row["Haplotype"].strip()
            if not PUBLIC.fullmatch(sample):
                continue
            all_haps.add((sample, hap))
            if row["Locus"] == LOCUS:
                cells[(sample, hap)].append(row)
                if row["Structure"] == "Provirus_from_Multi":
                    multi_source_rows += 1

    hap_rows: list[dict[str, object]] = []
    for sample, hap in sorted(all_haps):
        rows = cells.get((sample, hap), [])
        parts = []
        for row in rows:
            match = PART.search(row["ID_Full"])
            if match and row["Structure"] == "Provirus_from_Multi":
                parts.append(int(match.group(1)))
        if parts:
            copy_number, state = max(parts), "tandem_array"
        elif any(row["Structure"] == "Provirus" and "_asmdup" not in row["ID_Full"] for row in rows):
            copy_number, state = 1, "single_provirus"
        elif any(row["Structure"] == "Solo-LTR" for row in rows):
            copy_number, state = 0, "solo_ltr"
        else:
            copy_number, state = 0, ("absent" if not rows else "other_retained")
        hap_rows.append({
            "sample": sample,
            "haplotype": hap,
            "array_state": state,
            "array_copy_number": copy_number,
            "multi_vs_single_binary": int(copy_number >= 2),
            "high_copy_ge3_binary": int(copy_number >= 3),
        })
    write_tsv(OUT / "haplotype_copy_number_truth.tsv", hap_rows)

    hap_distribution = Counter(int(row["array_copy_number"]) for row in hap_rows)
    person_copy = defaultdict(int)
    person_multi_haps = defaultdict(int)
    person_high = defaultdict(int)
    for row in hap_rows:
        sample = str(row["sample"])
        person_copy[sample] += int(row["array_copy_number"])
        person_multi_haps[sample] += int(row["multi_vs_single_binary"])
        person_high[sample] = max(person_high[sample], int(row["high_copy_ge3_binary"]))

    distribution_rows: list[dict[str, object]] = []
    for value, n in sorted(hap_distribution.items()):
        distribution_rows.append({"target": "haplotype_array_copy_number", "value": value, "n": n})
    for value, n in sorted(Counter(person_copy.values()).items()):
        distribution_rows.append({"target": "person_total_array_units", "value": value, "n": n})
    for value, n in sorted(Counter(person_multi_haps.values()).items()):
        distribution_rows.append({"target": "person_multi_copy_haplotype_dosage", "value": value, "n": n})
    for value, n in sorted(Counter(person_high.values()).items()):
        distribution_rows.append({"target": "person_has_haplotype_copy_number_ge3", "value": value, "n": n})
    write_tsv(OUT / "copy_number_target_distributions.tsv", distribution_rows)

    # Freeze the already leakage-controlled validation metrics into one compact audit.
    proxy_rows: list[dict[str, object]] = []
    for filename, analysis, target, metric in (
        ("carrier_metrics.tsv", "nested_family_oof", "any_multi_copy_haplotype", "balanced_accuracy"),
        ("dosage_metrics.tsv", "nested_family_oof", "multi_copy_haplotype_dosage_0_1_2", "dosage_balanced_accuracy"),
        ("five_marker_carrier_metrics.tsv", "nested_family_oof_five_marker", "any_multi_copy_haplotype", "balanced_accuracy"),
        ("five_marker_dosage_metrics.tsv", "nested_family_oof_five_marker", "multi_copy_haplotype_dosage_0_1_2", "dosage_balanced_accuracy"),
    ):
        rows = load_rows(PREDICTOR / "results" / filename)
        row = next(r for r in rows if r["analysis"] == analysis and r["subset"] == "all")
        proxy_rows.append({
            "model": "single_flanking_snp" if "five_marker" not in analysis else "five_flanking_snp_panel",
            "validation": analysis,
            "target": target,
            "n": row["n_total"],
            "metric": metric,
            "estimate": row[metric],
            "ci_low": row.get(f"{metric}_ci_low", ""),
            "ci_high": row.get(f"{metric}_ci_high", ""),
            "screen_gate": "fail",
        })
    write_tsv(OUT / "proxy_validation_resolution.tsv", proxy_rows)

    summary = {
        "status": "coordinate blocker resolved; proxy performance blocker remains",
        "coordinate_resolution": {
            "authoritative_complete_array_hg38": f"chr7:{ARRAY[0]}-{ARRAY[1]}",
            "catalog_split_union_hg38": f"chr7:{old_union[0]}-{old_union[1]}",
            "constant_shift_bp": shift,
            "array_length_bp": ARRAY[1] - ARRAY[0] + 1,
            "shared_ltr_overlap_bp": 955,
            "interpretation": "the two catalog coordinates are tiled units of one reverse-strand tandem array, not competing loci",
        },
        "truth": {
            "people": len(person_copy),
            "haplotypes": len(hap_rows),
            "multi_copy_source_rows": multi_source_rows,
            "multi_copy_haplotypes": sum(int(row["multi_vs_single_binary"]) for row in hap_rows),
            "multi_copy_carrier_people": sum(value > 0 for value in person_multi_haps.values()),
            "max_haplotype_copy_number": max(hap_distribution),
            "max_person_total_array_units": max(person_copy.values()),
            "high_copy_ge3_people": sum(person_high.values()),
        },
        "decision": {
            "single_and_five_marker_proxy": "failed ancestry-portable validation; do not use for primary public phenotype screen",
            "direct_long_read_contrast": "usable now as multi-versus-single and as ordinal copy number 0-6",
            "copy_number_policy": "retain exact haplotype copy number 0,1,2,3,4,6; do not collapse to 0/1/2 multi-copy-haplotype dosage except for the separate public-panel proxy target",
        },
        "provenance": {
            "catalog_sha256": sha256(CATALOG),
            "vcf_sha256": sha256(VCF),
            "predictor_dosage_metrics_sha256": sha256(PREDICTOR / "results/dosage_metrics.tsv"),
            "predictor_carrier_metrics_sha256": sha256(PREDICTOR / "results/carrier_metrics.tsv"),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
