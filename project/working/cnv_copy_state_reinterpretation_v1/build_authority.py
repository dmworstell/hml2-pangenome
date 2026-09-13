#!/usr/bin/env python3
"""Build the deterministic retained-CNV copy-state reinterpretation authority."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from copy_state_authority import (
    STATE_ARTIFACT,
    STATE_TRUE,
    calibrate_positive_controls,
    score_candidate,
)


LANE = Path(__file__).resolve().parent
ROOT = LANE.parents[1]
SNAPSHOT = LANE / "source_snapshot"
RESULTS = LANE / "results"

DEPTH = SNAPSHOT / "cnv_depth_summary.tsv"
TARGETS = SNAPSHOT / "cnv_targets.tsv"
SEG_EVIDENCE = ROOT / "handoff/HML2_PARALOG_SYNTENY_PHYSICAL_COPY_EVIDENCE_20260718.tsv"
ARRAY_TRUTH = ROOT / "working/array_structural_reaudit_v1/haplotype_copy_truth.tsv"
ARRAY_ALL = ROOT / "working/catalog_cnv_array_inventory_agent/assembly_array_haplotypes.tsv"

EXPECTED_SHA256 = {
    DEPTH: "566cec9fc714c5b3f37f6cd892d800c87bcd46cd9b5b66fbe5a17497d0d12923",
    TARGETS: "bb2e05103d8303c6aea278cd8fb2052aaf2edbcc9ae38e0632426cd03867fca1",
    SEG_EVIDENCE: "3c1f26f19535b730e8d889b1f9917fb0d1f30611aa1639be79d8077170d387f3",
    ARRAY_TRUTH: "33af9333e053a200d70c85725b8e4ca5406ac11fa0ca01edf433dc1bee6299b9",
}

HG00423_KEY = "HG00423|HML-2_1q22|region2"
HG00423_ALT1 = "HG00423_pat_hprc_r2_v1.0.1_HML-2_1q22_alt1"
HG00423_ALT2 = "HG00423_pat_hprc_r2_v1.0.1_HML-2_1q22_alt2"
PRESERVED_WEIGHT_SHA256 = "597de702cbc995d39dd21ff51c5e549d4b8f83326d02ce708afb6dfd11d6bcda"
PRESERVED_CALIBRATION_SHA256 = "8b56d58036a3728069971bb167e5fc474c15f2e13bf8f86a32e7df7d917c50d8"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def key(row: dict[str, str]) -> str:
    return f"{row['sample']}|{row['locus']}|region{row['region']}"


def pair(row: dict[str, str]) -> tuple[str, str]:
    return row["sample"], row["locus"]


def is_unmeasured(row: dict[str, str]) -> bool:
    return row["low_baseline"] == "True" and row["hapcopies_all"] == ""


def dump_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def molecular_state_law(
    *, candidate_tip: str, parent_tip_candidates: list[str]
) -> dict[str, Any]:
    """Return the required sequence coupling for the two copy states."""
    return {
        "schema": "hml2.copy-state-molecular-coupling.v1",
        "independent_root_model_authorized": False,
        "observed_candidate_bytes_are_biological_truth_without_error_model": False,
        "states": {
            STATE_ARTIFACT: {
                "candidate_tip": candidate_tip,
                "candidate_bytes_role": "measurement_or_assembly_error_evidence_only",
                "biological_lineage_created": False,
                "included_in_biological_sequence_likelihood": False,
                "second_insertion_root_created": False,
            },
            STATE_TRUE: {
                "candidate_tip": candidate_tip,
                "candidate_bytes_role": "error-prone_observation_of_biological_daughter_lineage",
                "biological_lineage_created": True,
                "second_insertion_root_created": False,
                "parent_tip_candidates": parent_tip_candidates,
                "parent_assignment_latent": len(parent_tip_candidates) != 1,
                "duplication_birth_base_constraint": (
                    "daughter_birth_sequence_equals_parent_sequence_at_duplication_time"
                ),
                "duplication_event_sequence_error_sensitivity": {
                    "enabled": False,
                    "support_authority": None,
                    "activation_rule": (
                        "May be enabled only as an explicit sensitivity after independent evidence for "
                        "a rare duplication-event copying error is supplied; it is not a free default rate."
                    ),
                },
                "required_marginalization": [
                    "shared_ancestral_bases_at_duplication_birth",
                    "duplication_time",
                    "parent_post_duplication_mutation_history",
                    "daughter_post_duplication_mutation_history",
                    "post_duplication_gene_conversion_history",
                    "parent_sequence_measurement_error",
                    "daughter_sequence_measurement_error",
                ],
                "present_day_sequence_differences_rule": (
                    "Differences are generated only after the identical birth state through post-duplication "
                    "mutation or gene conversion and observation error, marginalized jointly with ancestry "
                    "and duplication time."
                ),
            },
        },
    }


def main() -> None:
    for path, expected in EXPECTED_SHA256.items():
        observed = sha256(path)
        if observed != expected:
            raise SystemExit(f"source digest changed: {path}: {observed}")

    target_pairs: set[tuple[str, str]] = set()
    for line in TARGETS.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 2:
            raise SystemExit("cnv_targets.tsv must remain a no-header two-column table")
        target_pairs.add((fields[0], fields[1]))
    if len(target_pairs) != 131:
        raise SystemExit("current target universe changed")

    depth_rows = read_tsv(DEPTH)
    historical_unmeasured = [row for row in depth_rows if is_unmeasured(row)]
    current_rows = [row for row in depth_rows if pair(row) in target_pairs]
    current_pairs = {pair(row) for row in current_rows}
    current_unmeasured = [row for row in current_rows if is_unmeasured(row)]
    if (len(depth_rows), len(historical_unmeasured), len(current_rows),
            len(current_pairs), len(current_unmeasured)) != (416, 68, 326, 120, 44):
        raise SystemExit("retained historical/current row counts changed")

    segdup_loci = {
        "HML-2_" + row["catalog_locus_id"]
        for row in read_tsv(SEG_EVIDENCE)
        if row["evidence_scope"] == "family_relationship"
        and "SEGDUP" in row["observed_structure"]
    }
    if len(segdup_loci) != 8:
        raise SystemExit("authenticated segmental-duplication locus set changed")

    # The all-locus catalog table supplies exact sample/locus array membership;
    # the independent truth table remains separately digest-bound above.  Only
    # rows with explicit accepted array membership are fixed as arrays.
    array_pairs = {
        (row["sample"], row["locus"])
        for row in read_tsv(ARRAY_ALL)
        if int(row["copy_number"]) >= 2 and row["part_number_complete"] == "1"
    }

    # Real, external controls only: authenticated segmental-duplication copies
    # whose retained local-flank measurement is present and historically
    # measurable.  Candidate rows do not enter calibration.  Body is treated
    # as cross-map ambiguous for both these duplicated controls and candidates.
    controls = [
        (key(row), row) for row in current_rows
        if row["locus"] in segdup_loci
        and not is_unmeasured(row)
        and row["flank_median_all"] != ""
    ]
    calibration = calibrate_positive_controls(
        controls, body_crossmap_ambiguous=True
    )
    if calibration.n_positive_controls != 61:
        raise SystemExit("real-control calibration row count changed")

    group_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in current_rows:
        group_rows[pair(row)].append(row)

    records: list[dict[str, Any]] = []
    for row in current_unmeasured:
        candidate_key = key(row)
        sample_locus = pair(row)
        base = {
            "candidate_key": candidate_key,
            "sample": row["sample"],
            "locus": row["locus"],
            "region": int(row["region"]),
            "historical_low_baseline": True,
            "retained_body_median_all": None if row["body_median_all"] == "" else float(row["body_median_all"]),
            "retained_local_flank_median_all": None if row["flank_median_all"] == "" else float(row["flank_median_all"]),
            "retained_sample_ref_cov": float(row["sample_ref_cov"]),
            "called_cn_is_catalog_row_count_not_depth_estimate": int(row["called_cn"]),
        }
        if row["locus"] in segdup_loci:
            base.update({
                "authority_class": "AUTHENTICATED_SEGMENTAL_DUPLICATION_FIXED",
                "assay_use_count_in_copy_state_weight": 0,
                "state_likelihoods": {"AUTHENTICATED_SEGMENTAL_DUPLICATION": 1.0},
                "copy_state_rule": "Known structural class is not unset by a low retained CNV baseline.",
            })
            records.append(base)
            continue
        if sample_locus in array_pairs:
            base.update({
                "authority_class": "AUTHENTICATED_ARRAY_FIXED",
                "assay_use_count_in_copy_state_weight": 0,
                "state_likelihoods": {"AUTHENTICATED_ARRAY": 1.0},
                "copy_state_rule": "Explicit accepted array annotation is not unset by a low retained CNV baseline.",
            })
            records.append(base)
            continue

        measured_tips = [
            key(other) for other in group_rows[sample_locus]
            if not is_unmeasured(other)
        ]
        all_candidate_tips = [
            key(other) for other in group_rows[sample_locus]
            if is_unmeasured(other)
            and other["locus"] not in segdup_loci
            and pair(other) not in array_pairs
        ]
        if candidate_key == HG00423_KEY:
            measured_tips = [HG00423_ALT1]
            candidate_tip = HG00423_ALT2
        else:
            candidate_tip = candidate_key

        sequence_law = molecular_state_law(
            candidate_tip=candidate_tip,
            parent_tip_candidates=measured_tips,
        )

        if row["flank_median_all"] == "":
            base.update({
                "authority_class": "UNINFORMATIVE_NO_LOCAL_FLANK_MEASUREMENT",
                "assay_use_count_in_copy_state_weight": 1,
                "body_crossmap_ambiguous": True,
                "state_likelihoods": None,
                "copy_state_rule": (
                    "A cross-map-ambiguous body without a retained local-flank measurement "
                    "cannot produce a normalized copy-state likelihood without inventing evidence or a prior."
                ),
                "state_contract": {
                    STATE_ARTIFACT: {
                        "candidate_tip_biological": False,
                        "molecular_coupling": sequence_law["states"][STATE_ARTIFACT],
                    },
                    STATE_TRUE: {
                        "candidate_tip_biological": True,
                        "molecular_coupling": sequence_law["states"][STATE_TRUE],
                    },
                },
                "molecular_coupling": sequence_law,
                "max_supported_candidate_tips": all_candidate_tips,
            })
            records.append(base)
            continue

        scored = score_candidate(
            candidate_key=candidate_key,
            row=row,
            calibration=calibration,
            always_biological_tips=measured_tips,
            candidate_tip=candidate_tip,
            terminal_host=row["sample"],
            body_crossmap_ambiguous=True,
        )
        sensitivity_favors = {
            name: max(weights, key=weights.get)
            for name, weights in scored["state_likelihoods"].items()
        }
        for state in (STATE_ARTIFACT, STATE_TRUE):
            scored["state_contract"][state]["molecular_coupling"] = (
                sequence_law["states"][state]
            )
        scored["molecular_coupling"] = sequence_law
        base.update(scored)
        base.update({
            "authority_class": "NON_SEGDUP_DUPLICATION_LATENT_CNV_WEIGHTED",
            "assay_use_count_in_copy_state_weight": 1,
            "sensitivity_favors": sensitivity_favors,
            "sensitivity_direction_stable": len(set(sensitivity_favors.values())) == 1,
            "max_supported_candidate_tips": all_candidate_tips,
        })
        records.append(base)

    class_counts = Counter(record["authority_class"] for record in records)
    if class_counts != Counter({
        "NON_SEGDUP_DUPLICATION_LATENT_CNV_WEIGHTED": 26,
        "UNINFORMATIVE_NO_LOCAL_FLANK_MEASUREMENT": 9,
        "AUTHENTICATED_SEGMENTAL_DUPLICATION_FIXED": 9,
    }):
        raise SystemExit(f"current classification counts changed: {class_counts}")

    hg = next(record for record in records if record["candidate_key"] == HG00423_KEY)
    if hg["metrics"]["region_support_ratio"] != 1 / 41:
        raise SystemExit("HG00423 retained unique-flank support changed")

    calibration_json = {
        "schema": "hml2.retained-cnv-real-control-calibration.v1",
        "calibration": calibration.as_dict(),
        "control_authority": {
            "class": "authenticated_segmental_duplication_physical_copies",
            "evidence_path": str(SEG_EVIDENCE.relative_to(ROOT)),
            "evidence_sha256": EXPECTED_SHA256[SEG_EVIDENCE],
            "candidate_rows_excluded": True,
            "body_crossmap_ambiguous": True,
            "calibration_feature": "local_flank_ratio_to_sample_diploid_target_flank_baseline",
        },
        "likelihood_model": {
            "state_centers": {STATE_ARTIFACT: 0.0, STATE_TRUE: 1.0},
            "kernel": "Gaussian radial likelihood with shared empirical scale",
            "normalized_weights": "likelihood divided by sum of the two state likelihoods",
            "prior_applied": False,
            "sensitivity_scales": [
                "narrow_q75_lower_residual",
                "central_rmse_lower_residual",
                "wide_q90_lower_residual",
            ],
        },
    }
    authority_json = {
        "schema": "hml2.retained-cnv-copy-state-authority.v1",
        "source_universes": {
            "historical_summary_rows": 416,
            "historical_unmeasured_rows": 68,
            "current_target_pairs": 131,
            "current_target_pairs_with_summary": 120,
            "current_summary_rows": 326,
            "current_unmeasured_rows": 44,
        },
        "classification_counts": dict(sorted(class_counts.items())),
        "cnv_owner": "this_authority_exactly_once_for_non_segdup_candidate_state_likelihood",
        "downstream_rule": (
            "Consumers marginalize state-conditional biology with central_state_likelihoods and must not "
            "reuse retained CNV values as a prior, likelihood, filter, or downstream weight."
        ),
        "molecular_coupling_rule": {
            "artifact_state": (
                "Candidate sequence bytes are measurement or assembly-error evidence only and create no "
                "second biological lineage."
            ),
            "true_duplication_state": (
                "Parent and daughter share identical base content at duplication birth; present differences "
                "require joint marginalization over ancestral bases, duplication time, post-duplication "
                "mutation/gene conversion, and sequence measurement error."
            ),
            "rare_duplication_event_error_default_enabled": False,
            "independent_root_model_authorized": False,
        },
        "records": records,
    }

    RESULTS.mkdir(exist_ok=True)
    dump_json(RESULTS / "real_control_calibration.v1.json", calibration_json)
    dump_json(RESULTS / "current_unmeasured_copy_state_authority.v1.json", authority_json)
    dump_json(RESULTS / "hg00423_1q22_copy_state_authority.v1.json", hg)

    weight_fields = [
        "candidate_key", "sample", "locus", "region", "authority_class",
        "body_median_all", "local_flank_median_all", "sample_ref_cov",
        "region_support_ratio", "artifact_weight_central", "true_duplication_weight_central",
        "central_likelihood_favors", "sensitivity_direction_stable",
    ]
    with (RESULTS / "current_unmeasured_copy_state_weights.v1.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=weight_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            metrics = record.get("metrics", {})
            weights = record.get("central_state_likelihoods", {})
            writer.writerow({
                "candidate_key": record["candidate_key"],
                "sample": record["sample"],
                "locus": record["locus"],
                "region": record["region"],
                "authority_class": record["authority_class"],
                "body_median_all": record["retained_body_median_all"],
                "local_flank_median_all": record["retained_local_flank_median_all"],
                "sample_ref_cov": record["retained_sample_ref_cov"],
                "region_support_ratio": metrics.get("region_support_ratio", ""),
                "artifact_weight_central": weights.get(STATE_ARTIFACT, ""),
                "true_duplication_weight_central": weights.get(STATE_TRUE, ""),
                "central_likelihood_favors": record.get("central_likelihood_favors", ""),
                "sensitivity_direction_stable": record.get("sensitivity_direction_stable", ""),
            })

    if sha256(RESULTS / "current_unmeasured_copy_state_weights.v1.tsv") != PRESERVED_WEIGHT_SHA256:
        raise SystemExit("molecular amendment changed normalized CNV weights")
    if sha256(RESULTS / "real_control_calibration.v1.json") != PRESERVED_CALIBRATION_SHA256:
        raise SystemExit("molecular amendment changed CNV calibration")

    amendment = {
        "schema": "hml2.copy-state-molecular-coupling-amendment.v1",
        "applies_to": {
            "current_authority_path": "results/current_unmeasured_copy_state_authority.v1.json",
            "current_authority_sha256": sha256(
                RESULTS / "current_unmeasured_copy_state_authority.v1.json"
            ),
            "hg00423_projection_path": "results/hg00423_1q22_copy_state_authority.v1.json",
            "hg00423_projection_sha256": sha256(
                RESULTS / "hg00423_1q22_copy_state_authority.v1.json"
            ),
        },
        "preserved_numeric_authorities": {
            "normalized_cnv_weights_sha256": PRESERVED_WEIGHT_SHA256,
            "real_control_calibration_sha256": PRESERVED_CALIBRATION_SHA256,
            "cnv_weights_changed": False,
        },
        "state_law": hg["molecular_coupling"],
        "consumer_requirements": {
            "artifact_candidate_enters_biological_tree": False,
            "true_duplication_tips_are_independent_roots": False,
            "true_duplication_birth_sequences_are_identical": True,
            "present_sequence_differences_are_observed_truth_without_error_model": False,
            "marginalize_molecular_coupling_with_copy_state": True,
        },
        "exact4_outputs_modified": False,
        "broad_cnv_rerun_or_submission_performed": False,
    }
    dump_json(RESULTS / "MOLECULAR_COUPLING_AMENDMENT.v1.json", amendment)

    current_keys = {key(row) for row in current_unmeasured}
    with (RESULTS / "historical_unmeasured_inventory.v1.tsv").open("w", newline="") as handle:
        fields = ["candidate_key", "sample", "locus", "region", "current_target", "current_authority_class"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        by_key = {record["candidate_key"]: record for record in records}
        for row in historical_unmeasured:
            row_key = key(row)
            writer.writerow({
                "candidate_key": row_key,
                "sample": row["sample"],
                "locus": row["locus"],
                "region": row["region"],
                "current_target": str(row_key in current_keys).upper(),
                "current_authority_class": by_key.get(row_key, {}).get(
                    "authority_class", "STALE_NOT_IN_CURRENT_131_TARGET_PAIR_UNIVERSE"
                ),
            })

    result_files = sorted(path for path in RESULTS.iterdir() if path.name != "RECEIPT.v1.json")
    receipt = {
        "schema": "hml2.retained-cnv-copy-state-reinterpretation-receipt.v1",
        "source_sha256": {
            str(path.relative_to(ROOT)): expected for path, expected in EXPECTED_SHA256.items()
        },
        "additional_array_membership_sha256": sha256(ARRAY_ALL),
        "result_sha256": {path.name: sha256(path) for path in result_files},
        "counts": {
            "historical_unmeasured": 68,
            "current_unmeasured": 44,
            "authenticated_fixed": 9,
            "latent_normalized": 26,
            "uninformative_missing_local_flank": 9,
            "real_calibration_controls": 61,
        },
        "rerun_or_submission_performed": False,
    }
    dump_json(RESULTS / "RECEIPT.v1.json", receipt)


if __name__ == "__main__":
    main()
