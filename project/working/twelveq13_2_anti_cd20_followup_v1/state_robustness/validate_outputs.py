#!/usr/bin/env python3
"""Independent contracts for the 12q13.2 anti-CD20 robustness package."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RESULTS = HERE / "results"
ANTI = ROOT / "project/working/direct_anti_cd20_viability_screen_v1/results"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    expected_files = {
        "accepted_orf_state_catalog.tsv",
        "adjusted_model_status.tsv",
        "adjusted_state_models.tsv",
        "all_individual_rows.tsv",
        "dose_rank_trend_tests.tsv",
        "individual_leave_one_out.tsv",
        "individual_leave_one_out_summary.tsv",
        "population_leave_one_out.tsv",
        "population_leave_one_out_summary.tsv",
        "raw_state_distributions.tsv",
        "source_manifest.tsv",
        "structural_omnibus_tests.tsv",
        "summary.json",
        "unadjusted_pairwise_tests.tsv",
    }
    observed_files = {path.name for path in RESULTS.iterdir() if path.is_file()}
    assert expected_files.issubset(observed_files), expected_files - observed_files

    people = pd.read_csv(RESULTS / "all_individual_rows.tsv", sep="\t")
    assert len(people) == 49 and people["sample"].nunique() == 49
    assert people["structural_class"].value_counts().to_dict() == {
        "internal_element_present": 33,
        "retained_noninternal_reference": 11,
        "biological_locus_absent": 5,
    }
    assert people["dose"].value_counts().sort_index().to_dict() == {0: 16, 1: 26, 2: 7}
    assert people["pro_orf_state"].value_counts().to_dict() == {
        "positive_carrier": 25,
        "third_state": 16,
        "known_negative_only": 8,
    }
    assert people["env_orf_state"].value_counts().to_dict() == {
        "positive_carrier": 32,
        "retained_noninternal_element": 11,
        "locus_absent": 5,
        "known_negative_only": 1,
    }

    internal = people["structural_class"] == "internal_element_present"
    retained = people["structural_class"] == "retained_noninternal_reference"
    absent = people["structural_class"] == "biological_locus_absent"
    assert (people.loc[internal, ["retained", "internal", "provirus"]] == 1).all().all()
    assert (people.loc[internal, "dose"] >= 1).all()
    assert (people.loc[retained, "retained"] == 1).all()
    assert (people.loc[retained, "internal"] == 0).all()
    assert people.loc[retained, "provirus"].isna().all()
    assert (people.loc[retained, "dose"] == 0).all()
    assert (people.loc[absent, ["retained", "internal", "provirus", "dose"]] == 0).all().all()
    assert people.loc[~internal, ["gag_orf_exposure", "pro_orf_exposure", "env_orf_exposure"]].isna().all().all()
    assert people.loc[internal, "gag_orf_exposure"].eq(1).all()
    assert people.loc[internal, "ordered_route_exposure"].eq(0).all()

    outcome_counts = {
        "Obin_media": 45,
        "Obin_serum": 49,
        "Ofat_media": 48,
        "Ofat_serum": 49,
        "Ritux_serum": 49,
    }
    assert {column: int(people[column].notna().sum()) for column in outcome_counts} == outcome_counts

    adjusted = pd.read_csv(RESULTS / "adjusted_state_models.tsv", sep="\t")
    status = pd.read_csv(RESULTS / "adjusted_model_status.tsv", sep="\t")
    assert len(adjusted) == 40 and len(status) == 40 and status["status"].eq("fitted").all()
    assert set(adjusted["outcome_id"]) == set(outcome_counts)
    assert adjusted.groupby("outcome_id")["model_id"].nunique().eq(8).all()
    for column in [
        "hc3_p",
        "superpopulation_fe_hc3_p",
        "permutation_p",
        "hc3_bh_all_adjusted_state_models",
        "permutation_bh_all_adjusted_state_models",
    ]:
        assert adjusted[column].between(0, 1).all(), column

    # Exact numerical recovery of the already-audited primary internal marker.
    focal = pd.read_csv(ANTI / "focal_results.tsv", sep="\t")
    focal = focal[focal["marker_id"] == "12q13.2::internal_presence::all_people"]
    recovered = adjusted[adjusted["model_id"] == "internal_vs_pooled_noninternal"].merge(
        focal[["outcome_id", "beta", "hc3_se", "hc3_p"]], on="outcome_id", suffixes=("_new", "_old")
    )
    assert len(recovered) == 5
    for column in ["beta", "hc3_se", "hc3_p"]:
        assert np.allclose(recovered[f"{column}_new"], recovered[f"{column}_old"], atol=1e-12, rtol=0)

    # Exact numerical recovery of both existing three-class contrasts.
    old_structural = pd.read_csv(ANTI / "twelveq13_posthoc_structural_contrasts.tsv", sep="\t")
    old_structural["model_id"] = old_structural["comparison"].map(
        {
            "internal_vs_retained_noninternal": "structural_internal_vs_retained_noninternal",
            "internal_vs_biological_absence": "structural_internal_vs_true_absence",
        }
    )
    recovered = adjusted.merge(
        old_structural[["outcome_id", "model_id", "beta", "hc3_se", "hc3_p"]],
        on=["outcome_id", "model_id"],
        suffixes=("_new", "_old"),
    )
    assert len(recovered) == 10
    for column in ["beta", "hc3_se", "hc3_p"]:
        assert np.allclose(recovered[f"{column}_new"], recovered[f"{column}_old"], atol=1e-12, rtol=0)

    pairwise = pd.read_csv(RESULTS / "unadjusted_pairwise_tests.tsv", sep="\t")
    assert len(pairwise) == 35
    assert pairwise["permutation_p"].between(0, 1).all()
    exact = pairwise[pairwise["permutation_method"] == "exact_all_label_allocations"]
    assert len(exact) >= 10 and exact["permutation_mc_se"].eq(0).all()

    rank = pd.read_csv(RESULTS / "dose_rank_trend_tests.tsv", sep="\t")
    omnibus = pd.read_csv(RESULTS / "structural_omnibus_tests.tsv", sep="\t")
    assert len(rank) == 5 and len(omnibus) == 5
    assert rank["permutation_p"].between(0, 1).all()
    assert omnibus["permutation_p"].between(0, 1).all()

    individual = pd.read_csv(RESULTS / "individual_leave_one_out.tsv", sep="\t")
    individual_summary = pd.read_csv(RESULTS / "individual_leave_one_out_summary.tsv", sep="\t")
    population = pd.read_csv(RESULTS / "population_leave_one_out.tsv", sep="\t")
    population_summary = pd.read_csv(RESULTS / "population_leave_one_out_summary.tsv", sep="\t")
    assert len(individual) == 1787 and len(individual_summary) == 40
    assert len(population) == 260 and len(population_summary) == 40
    assert individual["status"].eq("fitted").all()
    assert population["status"].eq("fitted").all()

    manifest = pd.read_csv(RESULTS / "source_manifest.tsv", sep="\t")
    assert len(manifest) == 7
    for row in manifest.itertuples(index=False):
        path = ROOT / row.path
        assert path.is_file() and sha256(path) == row.sha256

    print("PASS: 12q13.2 state semantics, 40 adjusted cells, exact recoveries, and robustness outputs validated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
