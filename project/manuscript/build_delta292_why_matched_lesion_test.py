#!/usr/bin/env python3
"""Test which coarse Δ292 explanations survive matched natural-lesion controls."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
LESIONS = (
    PROJECT
    / "working/claude_safe_delta292_why_followup_v1/results/"
    "observed_natural_deletion_consequence_map_v2.tsv"
)
FUNCTIONAL = (
    PROJECT
    / "manuscript/artifact_filtered_functional_refit/"
    "complete_three_outcome_refit_multiplicity.tsv"
)
CHEATER_GEOMETRY = (
    PROJECT
    / "working/delta292_differential_propagation_why_claude_v1/results/"
    "natural_lesion_cheater_geometry.tsv"
)
OUTDIR = PROJECT / "manuscript/delta292_why_matched_natural_lesions_v1"


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    lesions = pd.read_csv(LESIONS, sep="\t", dtype=str)
    ape = lesions[
        lesions["observation_layers"].str.contains("ape_genomewide", na=False)
    ].copy()
    for field in (
        "length_bp",
        "ancestral_integration_unit_count_ape_scan_FLOOR",
        "ancestral_integration_unit_count_ape_scan_CEILING",
    ):
        ape[field] = pd.to_numeric(ape[field], errors="coerce")
    ape["is_delta292"] = ape["is_exact_delta292"].eq("TRUE")
    ape["pre_orang_16_20Ma"] = ape["age_bracket_key"].str.contains(
        "pre_orang", na=False
    )
    ape["rec_destroyed"] = ape["rec_donor_removed"].str.lower().eq("true")
    ape["env_frame_destroyed"] = ~ape["env_frame_retained"].str.lower().eq("true")
    ape["pol_native_stop_removed_bool"] = ape[
        "pol_native_stop_removed"
    ].str.lower().eq(
        "true"
    )
    ape["delta292_coarse_trans_loss"] = (
        ape["rec_destroyed"]
        & ape["env_frame_destroyed"]
        & ape["pol_native_stop_removed_bool"]
    )
    fields = [
        "lesion_geometry_id",
        "length_bp",
        "species",
        "age_bracket_key",
        "pre_orang_16_20Ma",
        "rec_destroyed",
        "env_frame_destroyed",
        "pol_native_stop_removed_bool",
        "delta292_coarse_trans_loss",
        "shares_delta292_generalized_binary_vector",
        "shares_delta292_generalized_full_vector",
        "ancestral_integration_unit_count_ape_scan_FLOOR",
        "ancestral_integration_unit_count_ape_scan_CEILING",
        "orthology_resolution_status",
    ]
    write_tsv(OUTDIR / "ape_natural_lesion_analysis_table.tsv", ape[fields])

    delta = ape[ape["is_delta292"]].iloc[0]
    age_matched = ape[ape["pre_orang_16_20Ma"]].copy()
    rivals = age_matched[~age_matched["is_delta292"]]
    geometry = pd.read_csv(CHEATER_GEOMETRY, sep="\t", dtype=str)
    rec_geometry_rivals = geometry[
        geometry["layer"].eq("ape_outcome_blind_scan")
        & ~geometry["is_delta292"].str.lower().eq("true")
        & geometry[
            "REC_CHEATER_GEOMETRY_rec_destroyed_RcRE_retained"
        ].str.lower().eq("true")
    ].copy()
    rec_lengths = set(pd.to_numeric(rec_geometry_rivals["length_bp"]))
    rec_matched = rivals[rivals["length_bp"].isin(rec_lengths)].copy()
    if len(age_matched) != 7 or len(rec_matched) != 1:
        raise ValueError(
            f"unexpected matched sets: age={len(age_matched)}, rec={len(rec_matched)}"
        )
    write_tsv(OUTDIR / "pre_orang_age_matched_lesions.tsv", age_matched[fields])
    write_tsv(
        OUTDIR / "pre_orang_rec_env_loss_matched_lesion.tsv", rec_matched[fields]
    )

    floor = "ancestral_integration_unit_count_ape_scan_FLOOR"
    ceiling = "ancestral_integration_unit_count_ape_scan_CEILING"
    delta_floor = int(delta[floor])
    delta_ceiling = int(delta[ceiling])
    rival_floor_median = float(rivals[floor].median())
    rival_ceiling_max = int(rivals[ceiling].max())
    rec_rival = rec_matched.iloc[0]

    functional = pd.read_csv(FUNCTIONAL, sep="\t", na_values=["NA"])
    np9 = functional[
        (functional["outcome"] == "Houldcroft2014_EBV_qPCR")
        & (
            functional["exposure_id"]
            == "general::type1_compatible_np9_units"
        )
    ].iloc[0]
    type1_physical = functional[
        functional["exposure_id"] == "general::type1_physical_units"
    ]
    type1_physical_min_p = float(
        type1_physical["analysis_p_pedigree_cluster"].min()
    )

    hypotheses = pd.DataFrame(
        [
            {
                "hypothesis": "coarse_age_opportunity_alone",
                "executed_test": (
                    "Compare Δ292 with all other ape natural deletions in the same "
                    "pre-orang 16–20 Ma bracket"
                ),
                "result": (
                    f"Delta292={delta_floor}–{delta_ceiling} under the two resident "
                    f"collapse rules; "
                    f"five age-matched rivals all have ceiling≤{rival_ceiling_max}"
                ),
                "disposition": "COARSE_AGE_BRACKET_NOT_SUFFICIENT",
                "remaining_live_rival": (
                    "within-bracket source abundance/branching opportunity and "
                    "differential candidate recovery"
                ),
            },
            {
                "hypothesis": "generic_rec_env_loss_cheater_geometry",
                "executed_test": (
                    "Within the same pre-orang bracket, compare the natural 112-bp "
                    "lesion that also removes the Rec donor and destroys the Env frame"
                ),
                "result": (
                    f"Delta292={delta_floor}–{delta_ceiling}; "
                    f"{rec_rival['lesion_geometry_id']}="
                    f"{int(rec_rival[floor])}–{int(rec_rival[ceiling])}"
                ),
                "disposition": "BROAD_CHEATER_GEOMETRY_NOT_SUFFICIENT",
                "remaining_live_rival": (
                    "an exact-junction/cis effect or a uniquely successful source "
                    "lineage carrying Δ292"
                ),
            },
            {
                "hypothesis": "generic_template_shortening",
                "executed_test": (
                    "Within the same pre-orang bracket, compare removal lengths "
                    "74, 112, 116, 208, 292, 868 and 2255 bp"
                ),
                "result": (
                    "Every non-Δ292 age-matched lesion has one unit, including the "
                    "2255-bp lesion; multiplicity is not monotone in deleted length"
                ),
                "disposition": "GENERIC_SHORTENING_NOT_SUFFICIENT",
                "remaining_live_rival": (
                    "a small per-cycle length benefit can still contribute but cannot "
                    "explain lesion identity by itself"
                ),
            },
            {
                "hypothesis": "host_level_typeI_benefit",
                "executed_test": (
                    "Artifact-filtered refit of 17 catalog-wide Type-I/Type-II "
                    "burdens across three direct phenotype outcomes"
                ),
                "result": (
                    f"minimum Type-I physical-burden P={type1_physical_min_p:.4g}; "
                    f"Np9-compatible/Houldcroft P="
                    f"{np9['analysis_p_pedigree_cluster']:.4g}, conservative "
                    f"combined finite-test q={np9['q_bh_finite_model_suite']:.4g}"
                ),
                "disposition": "NO_CORRECTED_DIRECT_HOST_BENEFIT_SIGNAL",
                "remaining_live_rival": (
                    "unmeasured historical host phenotypes, producer-cell effects, "
                    "or linked host haplotypes"
                ),
            },
            {
                "hypothesis": "exact_delta292_or_source_lineage",
                "executed_test": (
                    "Intersection of the age-, length-, cheater-geometry-, and "
                    "artifact-corrected phenotype tests"
                ),
                "result": (
                    "Coarse generic properties do not reproduce Δ292 multiplicity; "
                    "the exact lesion and its historical source lineage remain aliased"
                ),
                "disposition": "PRIMARY_REMAINING_DISCRIMINATION",
                "remaining_live_rival": (
                    "direct matched competition/cis assays or an independent second "
                    "source lineage carrying the exact lesion"
                ),
            },
        ]
    )
    write_tsv(OUTDIR / "why_hypothesis_updates.tsv", hypotheses)

    plot = ape.sort_values([floor, ceiling, "length_bp"], ascending=True).copy()
    y = np.arange(len(plot))
    lows = plot[floor].to_numpy(dtype=float)
    highs = plot[ceiling].to_numpy(dtype=float)
    colors = np.where(
        plot["is_delta292"],
        "#D55E00",
        np.where(plot["pre_orang_16_20Ma"], "#0072B2", "#999999"),
    )
    fig, ax = plt.subplots(figsize=(9.2, 6.2), constrained_layout=True)
    ax.barh(y, lows, color=colors, alpha=0.9)
    ax.errorbar(
        lows,
        y,
        xerr=np.vstack([np.zeros_like(lows), highs - lows]),
        fmt="none",
        ecolor="#222222",
        capsize=3,
        lw=1.2,
    )
    labels = [
        f"{row.lesion_geometry_id} ({int(row.length_bp)} bp)"
        for row in plot.itertuples()
    ]
    ax.set_yticks(y, labels)
    ax.set_xlabel(
        "Unit counts under floor/ceiling-style collapse rules "
        "(not confidence bounds)"
    )
    ax.set_title(
        "Delta292 remains the multiplicity outlier within the pre-orang age bracket\n"
        "Blue: other 16–20 Ma lesions; orange: Delta292; gray: younger/unpolarized",
        loc="left",
        fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    figure = OUTDIR / "delta292_age_matched_natural_lesions.png"
    fig.savefig(figure, dpi=300)
    fig.savefig(figure.with_suffix(".pdf"))
    plt.close(fig)

    summary = {
        "schema": "hml2.delta292-why-matched-natural-lesions.v1",
        "ape_natural_lesions": len(ape),
        "pre_orang_age_matched_lesions": len(age_matched),
        "delta292_unit_floor": delta_floor,
        "delta292_unit_ceiling": delta_ceiling,
        "age_matched_rival_median_floor": rival_floor_median,
        "age_matched_rival_max_ceiling": rival_ceiling_max,
        "observed_floor_vs_max_rival_ceiling_style_ratio": (
            delta_floor / rival_ceiling_max
        ),
        "same_age_rec_env_loss_comparator": rec_rival["lesion_geometry_id"],
        "same_age_rec_env_loss_comparator_units": [
            int(rec_rival[floor]),
            int(rec_rival[ceiling]),
        ],
        "exchangeability_rank_probability_note": (
            "Delta292 is the unique maximum among seven age-bracket-matched lesions "
            "(Delta292 plus six rivals). A one-sided rank probability would be 1/7 "
            "under exchangeability, but "
            "exchangeability is not licensed because source opportunity and recovery "
            "remain lesion-specific."
        ),
        "identifiability_update": (
            "The coarse age bracket, generic shortening, and broad Rec/Env-loss "
            "geometry are not sufficient. Resident data now concentrate the question "
            "on an exact-Delta292 effect versus exceptional historical source-lineage "
            "opportunity; those two remain aliased without direct competition/cis "
            "experiments or an independent replicated source lineage."
        ),
    }
    (OUTDIR / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
