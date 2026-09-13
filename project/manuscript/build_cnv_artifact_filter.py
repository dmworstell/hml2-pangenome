#!/usr/bin/env python3
"""Materialize the manuscript CNV artifact exclusions and filtered ORF input."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SOURCE_ORF = (
    PROJECT
    / "inputs/orf_analysis/current_v3r1/"
    "combined_hml2_orf_analysis.CNV_WEIGHTED.v3r1.tsv"
)
WEIGHTS = (
    PROJECT
    / "working/cnv_depth_weight_finalization_v1/results_44/"
    "final_assayed_candidate_weights.v1.tsv"
)
RECEIPT = (
    PROJECT
    / "working/cnv_recovery_status/final_v8_results/"
    "LANE_G_CNV_V8_CONSUMER_RELEASE_RECEIPT.v1.json"
)
SUPPLEMENT = (
    PROJECT / "manuscript/supplement/Table_S4_CNV_assembly_artifact_qc.tsv"
)
EXCLUDED_ROWS = (
    PROJECT / "manuscript/supplement/Table_S_CNV_excluded_ORF_records.tsv"
)
FILTERED_ORF = (
    PROJECT
    / "inputs/orf_analysis/current_v3r1/"
    "combined_hml2_orf_analysis.CNV_WEIGHTED.v3r1.ARTIFACT_FILTERED.tsv.gz"
)
MANIFEST = FILTERED_ORF.with_suffix(".manifest.json")
EXPECTED_SOURCE_SHA256 = (
    "1a3f0462b541fa14c5e1f882137d323b3e5418daa88813dc9bcbc9205d43259d"
)
ASM_DUP_ID = "HG00658_pat_hprc_r2_v1.0.1_HML-2_7p22.1_asmdup"
ASM_DUP_AUTHORITY = (
    PROJECT
    / "working/sevenp22_proxy_resolution_agent/haplotype_copy_number_truth.tsv"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if sha256(SOURCE_ORF) != EXPECTED_SOURCE_SHA256:
        raise ValueError("terminal v3r1 ORF input changed")

    receipt = json.loads(RECEIPT.read_text())
    id_by_candidate = receipt["authority_exact_controller_replay_map"]

    with WEIGHTS.open(newline="") as handle:
        weights = list(csv.DictReader(handle, delimiter="\t"))
    if len(weights) != 44:
        raise ValueError(f"expected 44 terminal CNV authority rows, found {len(weights)}")

    supplement_fields = [
        "candidate_key",
        "ID_Full",
        "sample",
        "locus",
        "region",
        "authority_class",
        "measurement_source",
        "raw_body_median",
        "raw_local_flank_median",
        "sample_global_diploid_target_flank_baseline",
        "assembly_artifact_probability",
        "later_duplication_probability",
        "authenticated_segdup_probability",
        "dominant_state",
        "source_hard_call",
        "manuscript_analysis_decision",
        "decision_scope",
    ]
    supplement_rows = []
    excluded_ids: set[str] = set()
    for row in weights:
        candidate = row["candidate_key"]
        id_full = id_by_candidate[candidate]
        artifact = row["dominant_state"] == "assembly_artifact"
        if artifact:
            excluded_ids.add(id_full)
        supplement_rows.append(
            {
                "candidate_key": candidate,
                "ID_Full": id_full,
                "sample": row["sample"],
                "locus": row["locus"],
                "region": row["region"],
                "authority_class": row["authority_class"],
                "measurement_source": row["measurement_source"],
                "raw_body_median": row["raw_body_median"],
                "raw_local_flank_median": row["raw_local_flank_median"],
                "sample_global_diploid_target_flank_baseline": row[
                    "sample_global_diploid_target_flank_baseline"
                ],
                "assembly_artifact_probability": row["assembly_artifact"],
                "later_duplication_probability": row["later_duplication"],
                "authenticated_segdup_probability": row["authenticated_segdup"],
                "dominant_state": row["dominant_state"],
                "source_hard_call": row["hard_call"],
                "manuscript_analysis_decision": (
                    "EXCLUDE_EXACT_ASSEMBLY_RECORD"
                    if artifact
                    else "RETAIN_EXACT_ASSEMBLY_RECORD"
                ),
                "decision_scope": (
                    "Exact ID_Full only; never a locus-wide exclusion"
                ),
            }
        )

    # The terminal candidate-tip model does not include the previously
    # adjudicated HG00658 paternal 7p22.1 assembly-gap duplicate.  Exact
    # coordinate truth sets that haplotype to CN0, so this one hard artifact is
    # added to the same manuscript exclusion ledger rather than silently
    # persisting through the newer table.
    excluded_ids.add(ASM_DUP_ID)
    supplement_rows.append(
        {
            "candidate_key": "HG00658|HML-2_7p22.1|asmdup",
            "ID_Full": ASM_DUP_ID,
            "sample": "HG00658",
            "locus": "HML-2_7p22.1",
            "region": "asmdup",
            "authority_class": "EXACT_COORDINATE_ARRAY_TRUTH",
            "measurement_source": str(ASM_DUP_AUTHORITY.relative_to(PROJECT.parent)),
            "raw_body_median": "",
            "raw_local_flank_median": "",
            "sample_global_diploid_target_flank_baseline": "",
            "assembly_artifact_probability": "1",
            "later_duplication_probability": "0",
            "authenticated_segdup_probability": "0",
            "dominant_state": "assembly_artifact",
            "source_hard_call": "True",
            "manuscript_analysis_decision": "EXCLUDE_EXACT_ASSEMBLY_RECORD",
            "decision_scope": "Exact ID_Full only; never a locus-wide exclusion",
        }
    )

    SUPPLEMENT.parent.mkdir(parents=True, exist_ok=True)
    with SUPPLEMENT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=supplement_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(supplement_rows)

    excluded_fields = [
        "ID_Full",
        "Locus",
        "ID",
        "Haplotype",
        "Structure",
        "observation_state",
        "observation_state_reason",
        "observation_source_kind",
        "copy_observation_class",
        "dominant_copy_state",
        "assembly_artifact_probability",
    ]
    excluded_probability = {
        row["ID_Full"]: row["assembly_artifact_probability"]
        for row in supplement_rows
        if row["manuscript_analysis_decision"] == "EXCLUDE_EXACT_ASSEMBLY_RECORD"
    }
    found_ids: set[str] = set()
    excluded_rows: list[dict[str, str]] = []
    source_rows = 0
    retained_rows = 0

    FILTERED_ORF.parent.mkdir(parents=True, exist_ok=True)
    with SOURCE_ORF.open(newline="") as source, FILTERED_ORF.open("wb") as raw:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None or "ID_Full" not in reader.fieldnames:
            raise ValueError("ORF table lacks ID_Full")
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(
                    text,
                    fieldnames=reader.fieldnames,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in reader:
                    source_rows += 1
                    id_full = row["ID_Full"]
                    if id_full in excluded_ids:
                        found_ids.add(id_full)
                        excluded_rows.append(
                            {
                                **{
                                    field: row.get(field, "")
                                    for field in excluded_fields[:-2]
                                },
                                "dominant_copy_state": "assembly_artifact",
                                "assembly_artifact_probability": excluded_probability[
                                    id_full
                                ],
                            }
                        )
                        continue
                    writer.writerow(row)
                    retained_rows += 1
                text.flush()

    missing_ids = sorted(excluded_ids - found_ids)
    if missing_ids:
        raise ValueError(f"artifact IDs absent from v3r1 ORF table: {missing_ids}")

    with EXCLUDED_ROWS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=excluded_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(excluded_rows)

    payload = {
        "schema": "hml2.manuscript-cnv-artifact-filter.v1",
        "policy": (
            "Exclude each terminal CNV candidate whose dominant_state is "
            "assembly_artifact plus the coordinate-adjudicated HG00658 7p22.1 "
            "assembly-gap duplicate. Every exclusion applies only to the exact "
            "ID_Full; valid copies at the same locus and sample remain."
        ),
        "source_rows": source_rows,
        "retained_rows": retained_rows,
        "excluded_orf_rows": len(excluded_rows),
        "excluded_unique_ID_Full": len(excluded_ids),
        "source_hard_call_status": (
            "34 candidate-tip exclusions use final materialized-assay weights; "
            "1 exact-coordinate 7p22.1 exclusion is a hard structural call"
        ),
        "source_weight_status": "MIXED_WEIGHTED_AND_COORDINATE_ADJUDICATED",
        "source_calibration_status": (
            "34 numerical artifact probabilities remain provisional; the "
            "HG00658 7p22.1 duplicate is coordinate-adjudicated"
        ),
        "inputs": {
            str(SOURCE_ORF.relative_to(PROJECT.parent)): sha256(SOURCE_ORF),
            str(WEIGHTS.relative_to(PROJECT.parent)): sha256(WEIGHTS),
            str(RECEIPT.relative_to(PROJECT.parent)): sha256(RECEIPT),
            str(ASM_DUP_AUTHORITY.relative_to(PROJECT.parent)): sha256(
                ASM_DUP_AUTHORITY
            ),
        },
        "outputs": {
            str(FILTERED_ORF.relative_to(PROJECT.parent)): sha256(FILTERED_ORF),
            str(SUPPLEMENT.relative_to(PROJECT.parent)): sha256(SUPPLEMENT),
            str(EXCLUDED_ROWS.relative_to(PROJECT.parent)): sha256(EXCLUDED_ROWS),
        },
    }
    MANIFEST.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(FILTERED_ORF)
    print(SUPPLEMENT)
    print(EXCLUDED_ROWS)


if __name__ == "__main__":
    main()
