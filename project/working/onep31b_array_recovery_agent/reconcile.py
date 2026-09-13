#!/usr/bin/env python3
"""Read-only reconciliation of HML-2_1p31.1b tandem arrays."""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "HML2_project_data"
OUT = ROOT / "project/working/onep31b_array_recovery_agent"
CATALOG = ROOT / "HML2_ProjectResources/data/catalog/combined_hml2_orf_analysis.tsv"
STATES = ROOT / "manuscript_figures/python/analysis/catalog_functional_screen/results/locus_haplotype_states.tsv"
TRUTH = ROOT / "manuscript_figures/python/analysis/catalog_functional_screen/results/sample_contrast_truth.tsv"
LOCUS = "HML-2_1p31.1b"
PART = re.compile(r"_MULTI_part(\d+)$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_fasta(stem: str) -> Path:
    candidates = [
        DATA / "processed_loci" / LOCUS / f"{stem}.fa",
        DATA / "hgsvc3-2024-02-23-mc-chm13" / LOCUS / f"{stem}.fa",
    ]
    matches = [path for path in candidates if path.exists()]
    if len(matches) != 1:
        raise RuntimeError(f"expected one FASTA for {stem}, found {matches}")
    return matches[0]


def fasta_record(path: Path) -> tuple[str, int]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        header = handle.readline().rstrip("\n")
        length = sum(len(line.strip()) for line in handle)
    return header, length


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    population = {}
    with TRUTH.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            population.setdefault(row["sample"], (row["superpopulation"], row["population"]))
    rows = []
    with CATALOG.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            match = PART.search(row["ID_Full"])
            if row["Locus"] == LOCUS and match:
                row["copy_index"] = int(match.group(1))
                rows.append(row)

    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault((row["ID"], row["Haplotype"]), []).append(row)

    source_by_cell = {}
    for key, cell_rows in grouped.items():
        cell_rows.sort(key=lambda row: int(row["copy_index"]))
        expected = list(range(1, len(cell_rows) + 1))
        observed = [int(row["copy_index"]) for row in cell_rows]
        if observed != expected:
            raise RuntimeError(f"non-contiguous copy indices for {key}: {observed}")
        stem = PART.sub("", cell_rows[0]["ID_Full"])
        fasta = find_fasta(stem)
        header, length = fasta_record(fasta)
        if "_MULTI" not in header:
            raise RuntimeError(f"array FASTA lacks _MULTI tag: {fasta}")
        source_by_cell[key] = (fasta, header, length)

    detail_fields = [
        "locus", "sample", "haplotype", "superpopulation", "population", "copy_number", "copy_index",
        "id_full", "source_identifier", "structure",
        "gag", "missense_gag", "pro", "missense_pro", "pol", "missense_pol",
        "env", "missense_env", "np9", "missense_np9", "rec",
        "fasta_path", "fasta_header", "fasta_sequence_length",
    ]
    with (OUT / "copy_specific_orf_reconciliation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=detail_fields, delimiter="\t")
        writer.writeheader()
        for row in sorted(rows, key=lambda x: (x["ID"], x["Haplotype"], int(x["copy_index"]))):
            key = (row["ID"], row["Haplotype"])
            fasta, header, length = source_by_cell[key]
            writer.writerow({
                "locus": LOCUS,
                "sample": row["ID"],
                "haplotype": row["Haplotype"],
                "superpopulation": population.get(row["ID"], ("UNK", "UNK"))[0],
                "population": population.get(row["ID"], ("UNK", "UNK"))[1],
                "copy_number": len(grouped[key]),
                "copy_index": row["copy_index"],
                "id_full": row["ID_Full"],
                "source_identifier": row["Source_Identifier"],
                "structure": row["Structure"],
                "gag": row["gag"], "missense_gag": row["missense_gag"],
                "pro": row["pro"], "missense_pro": row["missense_pro"],
                "pol": row["pol"], "missense_pol": row["missense_pol"],
                "env": row["env"], "missense_env": row["missense_env"],
                "np9": row["np9"], "missense_np9": row["missense_np9"],
                "rec": row["rec"],
                "fasta_path": str(fasta), "fasta_header": header,
                "fasta_sequence_length": length,
            })

    cell_fields = [
        "locus", "sample", "haplotype", "superpopulation", "population", "copy_number", "copy_indices",
        "array_unit_rows", "all_units_structure", "compatible_env_units",
        "compatible_np9_units", "broken_gag_units", "broken_pro_units",
        "broken_pol_units", "fasta_path", "fasta_sequence_length",
    ]
    with (OUT / "haplotype_array_reconciliation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cell_fields, delimiter="\t")
        writer.writeheader()
        for key, cell_rows in sorted(grouped.items()):
            fasta, _header, length = source_by_cell[key]
            writer.writerow({
                "locus": LOCUS, "sample": key[0], "haplotype": key[1],
                "superpopulation": population.get(key[0], ("UNK", "UNK"))[0],
                "population": population.get(key[0], ("UNK", "UNK"))[1],
                "copy_number": len(cell_rows),
                "copy_indices": ";".join(str(row["copy_index"]) for row in cell_rows),
                "array_unit_rows": len(cell_rows),
                "all_units_structure": ";".join(sorted({row["Structure"] for row in cell_rows})),
                "compatible_env_units": sum(row["env"].lower() in {"intact", "intact_fs_end"} for row in cell_rows),
                "compatible_np9_units": sum(row["np9"].lower() in {"intact", "intact_fs_end"} for row in cell_rows),
                "broken_gag_units": sum(row["gag"].lower() not in {"intact", "intact_fs_end"} for row in cell_rows),
                "broken_pro_units": sum(row["pro"].lower() not in {"intact", "intact_fs_end"} for row in cell_rows),
                "broken_pol_units": sum(row["pol"].lower() not in {"intact", "intact_fs_end"} for row in cell_rows),
                "fasta_path": str(fasta), "fasta_sequence_length": length,
            })

    state_row = None
    with STATES.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["locus"] == LOCUS:
                state_row = row
                break
    if state_row is None:
        raise RuntimeError("downstream locus state row missing")

    copy_counts = Counter(len(cell_rows) for cell_rows in grouped.values())
    source_workbook = ROOT / "project/manuscript/source_snapshot/HML2_Supplementary_Data_S1.xlsx"
    summary = {
        "catalog_path": str(CATALOG),
        "catalog_sha256": sha256(CATALOG),
        "source_workbook": str(source_workbook),
        "source_workbook_sha256": sha256(source_workbook),
        "source_workbook_sheet": "ORF_analysis_full",
        "n_array_haplotypes": len(grouped),
        "n_array_unit_rows": len(rows),
        "copy_number_counts": dict(sorted(copy_counts.items())),
        "max_haplotype_copy_number": max(copy_counts),
        "all_array_units_structure": sorted({row["Structure"] for row in rows}),
        "all_array_units_env_compatible": all(row["env"].lower() in {"intact", "intact_fs_end"} for row in rows),
        "all_array_units_np9_compatible": all(row["np9"].lower() in {"intact", "intact_fs_end"} for row in rows),
        "all_array_units_gag_broken": all(row["gag"].lower() not in {"intact", "intact_fs_end"} for row in rows),
        "all_array_units_pro_broken": all(row["pro"].lower() not in {"intact", "intact_fs_end"} for row in rows),
        "all_array_units_pol_broken": all(row["pol"].lower() not in {"intact", "intact_fs_end"} for row in rows),
        "downstream_state_counts": {key: state_row[key] for key in ("Multi-copy", "Provirus", "Solo-LTR", "Fragment", "Absent")},
        "corrected_state_counts": {"Multi-copy": 9, "Provirus": 0, "Solo-LTR": 469, "Fragment": 106, "Absent": 0},
    }
    import json
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
