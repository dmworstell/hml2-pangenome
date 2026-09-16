#!/usr/bin/env python3
"""Refit catalog-wide functional burdens after exact assembly-artifact removal."""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(__file__).resolve().parents[1]
CATALOG = (
    PROJECT
    / "results/short_orf_rule_correction_20260915/"
    "combined_hml2_orf_analysis.RESOLVED.SHORT_ORF_CORRECTED.tsv"
)
OLD_MATRIX = (
    PROJECT
    / "working/direct_ebv_fitness_screen_v1/results/person_level_direct_matrix.tsv"
)
OLD_MANDAGE_MODELS = (
    PROJECT / "working/direct_ebv_fitness_screen_v1/results/model_results.tsv"
)
OLD_SECONDARY_MODELS = (
    PROJECT
    / "working/direct_ebv_fitness_screen_v1/results/"
    "secondary_outcome_model_results.tsv"
)
HOLD = (
    PROJECT
    / "working/data_search_direct_cellular_phenotypes_v1/derived/"
    "Houldcroft2014_EBV_qPCR_unique_direct292_join.tsv"
)
IM = (
    PROJECT
    / "working/data_search_direct_cellular_phenotypes_v1/derived/"
    "Im2012_LCL_intrinsic_growth_direct292_join.tsv"
)
EXCLUSIONS = PROJECT / "manuscript/supplement/Table_S_CNV_excluded_ORF_records.tsv"
SEVENP22 = (
    PROJECT / "working/sevenp22_proxy_resolution_agent/haplotype_copy_number_truth.tsv"
)
ONEQ22 = PROJECT / "manuscript/supplement/Table_S9_1q22_Gag_artifact_corrected_truth.tsv"
OUTDIR = PROJECT / "manuscript/artifact_filtered_functional_refit"
PUBLIC_ID = re.compile(r"^(?:HG|NA)\d+$")
COMPATIBLE = {"Intact", "Frameshift_at_end", "Intact_FS_End"}
PROVIRUS = {"Provirus", "Provirus_from_Multi"}

GENERAL_FEATURES = [
    "general::type1_physical_units",
    "general::type2_physical_units",
    "general::type1_loci_present",
    "general::type2_loci_present",
    "general::type1_fraction",
    "general::type1_type2_log_ratio",
    "general::type1_minus_type2_units",
    "general::type1_compatible_gag_units",
    "general::type1_compatible_pro_units",
    "general::type1_compatible_pol_units",
    "general::type1_compatible_env_units",
    "general::type1_compatible_np9_units",
    "general::type2_compatible_gag_units",
    "general::type2_compatible_pro_units",
    "general::type2_compatible_pol_units",
    "general::type2_compatible_env_units",
    "general::type2_compatible_rec_units",
]
LOCUS_FEATURES = [
    "HML-2_4p16.3a::structural::solo_ltr_vs_absent",
    "HML-2_15q25.2::structural::solo_ltr_vs_absent",
    "HML-2_1q22::orf::orf_gag",
    "HML-2_7p22.1::structural::multi_vs_single",
]


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def load_provirus() -> pd.DataFrame:
    rows = [
        row
        for row in load_catalog_rows()
        if row["observation_state"] == "PRESENT"
        and row["Structure"] in PROVIRUS
        and row["provirus_type"] in {"type1", "type2"}
    ]
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("no provirus rows in filtered catalog")
    return frame


def load_catalog_rows() -> list[dict[str, str]]:
    with CATALOG.open(newline="", encoding="utf-8") as handle:
        public_rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if PUBLIC_ID.fullmatch(row["ID"])
        ]
    excluded = {
        "non_HML2_HML11_sequence_identity": 1168,
        "alias_duplicate_of_8q24.3c": 584,
        "assembly_artifact_not_supported_by_CNV_depth": 35,
        "duplicate_catalog_label_for_same_assembled_interval": 78,
    }
    observed = {
        reason: sum(
            row["analysis_include"] == "0"
            and row["analysis_exclusion_reason"] == reason
            for row in public_rows
        )
        for reason in excluded
    }
    if observed != excluded:
        raise ValueError(f"unexpected biological exclusions: {observed}")
    rows = [row for row in public_rows if row["analysis_include"] == "1"]
    if len(rows) != 59_656:
        raise ValueError(f"unexpected analysis row count: {len(rows)}")
    return rows


def compatible(frame: pd.DataFrame, field: str) -> pd.Series:
    return frame[field].isin(COMPATIBLE)


def derive_exposures(prov: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    result = pd.DataFrame({"sample": samples}).set_index("sample")
    for provirus_type in ("type1", "type2"):
        sub = prov[prov["provirus_type"] == provirus_type].copy()
        units = sub.groupby("ID").size().reindex(samples, fill_value=0)
        loci = sub.groupby("ID")["Locus"].nunique().reindex(samples, fill_value=0)
        result[f"general::{provirus_type}_physical_units"] = units
        result[f"general::{provirus_type}_loci_present"] = loci
        for gene in ("gag", "pro", "env"):
            values = (
                sub.assign(_compatible=compatible(sub, gene))
                .groupby("ID")["_compatible"]
                .sum()
                .reindex(samples, fill_value=0)
            )
            result[f"general::{provirus_type}_compatible_{gene}_units"] = values
        pol_route = (
            compatible(sub, "gag")
            & compatible(sub, "pro")
            & compatible(sub, "pol")
        )
        result[f"general::{provirus_type}_compatible_pol_units"] = (
            sub.assign(_compatible=pol_route)
            .groupby("ID")["_compatible"]
            .sum()
            .reindex(samples, fill_value=0)
        )
        accessory = "np9" if provirus_type == "type1" else "rec"
        result[f"general::{provirus_type}_compatible_{accessory}_units"] = (
            sub.assign(_compatible=compatible(sub, accessory))
            .groupby("ID")["_compatible"]
            .sum()
            .reindex(samples, fill_value=0)
        )

    t1 = result["general::type1_physical_units"]
    t2 = result["general::type2_physical_units"]
    result["general::type1_fraction"] = t1 / np.maximum(t1 + t2, 1)
    result["general::type1_type2_log_ratio"] = np.log2((t1 + 0.5) / (t2 + 0.5))
    result["general::type1_minus_type2_units"] = t1 - t2
    return result.reset_index()[["sample", *GENERAL_FEATURES]]


def derive_locus_exposures(
    rows: list[dict[str, str]], samples: list[str]
) -> pd.DataFrame:
    by_cell: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    haps: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        by_cell[(row["Locus"], row["ID"], row["Haplotype"])].append(row)
        haps[row["ID"]].add(row["Haplotype"])

    result = pd.DataFrame({"sample": samples}).set_index("sample")
    for locus in ("HML-2_4p16.3a", "HML-2_15q25.2"):
        values = {}
        for sample in samples:
            states = []
            for hap in sorted(haps[sample]):
                cell = by_cell.get((locus, sample, hap), [])
                present = [
                    row for row in cell if row["observation_state"] == "PRESENT"
                ]
                if any(row["Structure"] == "Solo-LTR" for row in present):
                    states.append("focal")
                elif present:
                    states.append("third")
                elif any(
                    row["observation_state"] == "UNKNOWN_TECHNICAL" for row in cell
                ):
                    states.append("missing")
                else:
                    states.append("reference")
            values[sample] = (
                sum(state == "focal" for state in states)
                if all(state in {"focal", "reference"} for state in states)
                else np.nan
            )
        result[f"{locus}::structural::solo_ltr_vs_absent"] = pd.Series(values)

    oneq = pd.read_csv(ONEQ22, sep="\t")
    result["HML-2_1q22::orf::orf_gag"] = (
        oneq.set_index("sample")["gag_compatible_dosage"].reindex(samples)
    )
    seven = pd.read_csv(SEVENP22, sep="\t")
    if seven.duplicated(["sample", "haplotype"]).any():
        raise ValueError("duplicate 7p22.1 sample-haplotype calls")
    if not seven.groupby("sample").size().eq(2).all():
        raise ValueError("7p22.1 dosage requires two haplotype records per donor")
    seven["multi"] = (seven["array_copy_number"] >= 2).astype(float).where(
        seven["array_copy_number"].notna()
    )
    result["HML-2_7p22.1::structural::multi_vs_single"] = (
        seven.groupby("sample")["multi"].sum(min_count=2).reindex(samples)
    )
    return result.reset_index()[["sample", *LOCUS_FEATURES]]


def design_matrix(
    data: pd.DataFrame, exposure: str, adjustment: str | None, ancestry: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    fields = ["outcome", exposure, "sex", ancestry, "pedigree_component"]
    if adjustment:
        fields.append(adjustment)
    d = data[fields].dropna().copy()
    if len(d) < 20 or d[exposure].nunique() < 2:
        raise ValueError("not_estimable")
    columns = [
        pd.Series(1.0, index=d.index, name="intercept"),
        d[exposure].astype(float).rename("x"),
    ]
    if adjustment:
        columns.append(d[adjustment].astype(float).rename("adjustment"))
    for field in (ancestry, "sex"):
        dummy = pd.get_dummies(d[field].astype(str), prefix=field, drop_first=True)
        columns.extend(dummy[column].astype(float) for column in dummy.columns)
    design = pd.concat(columns, axis=1)
    x = design.to_numpy(dtype=float)
    y = d["outcome"].to_numpy(dtype=float)
    clusters = d["pedigree_component"].astype(str).to_numpy()
    return x, y, clusters, len(d)


def cluster_fit(
    data: pd.DataFrame, exposure: str, adjustment: str | None, ancestry: str
) -> dict[str, object]:
    try:
        x, y, clusters, n = design_matrix(data, exposure, adjustment, ancestry)
    except ValueError:
        return {"status": "not_estimable", "n": 0}
    rank = np.linalg.matrix_rank(x)
    if rank != x.shape[1]:
        return {"status": "singular", "n": n}
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    bread = np.linalg.inv(x.T @ x)
    unique_clusters = np.unique(clusters)
    meat = np.zeros((x.shape[1], x.shape[1]))
    for cluster in unique_clusters:
        idx = clusters == cluster
        score = x[idx].T @ residual[idx]
        meat += np.outer(score, score)
    g = len(unique_clusters)
    correction = (g / (g - 1)) * ((n - 1) / (n - rank))
    covariance = bread @ meat @ bread * correction
    se = math.sqrt(max(float(covariance[1, 1]), 0.0))
    statistic = float(beta[1]) / se if se else math.inf
    pvalue = 2 * stats.t.sf(abs(statistic), df=g - 1)
    crit = stats.t.ppf(0.975, df=g - 1)
    return {
        "status": "fit",
        "n": n,
        "clusters": g,
        "beta": float(beta[1]),
        "se_cluster": se,
        "p_cluster": float(pvalue),
        "ci95_low": float(beta[1] - crit * se),
        "ci95_high": float(beta[1] + crit * se),
        "standardized_beta": float(beta[1] * np.std(x[:, 1], ddof=1) / np.std(y, ddof=1)),
    }


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index)
    finite = values.dropna().sort_values()
    if not finite.between(0, 1).all():
        raise ValueError("BH requires P values in [0, 1] or missing values")
    m = len(finite)
    if not m:
        return result
    adjusted = np.minimum.accumulate(
        (finite.to_numpy() * m / np.arange(1, m + 1))[::-1]
    )[::-1]
    result.loc[finite.index] = np.minimum(adjusted, 1)
    return result


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(OLD_MATRIX, sep="\t", na_values=["NA"])
    samples = old["sample"].tolist()
    prov = load_provirus()
    exposures = derive_exposures(prov, samples)
    write_tsv(OUTDIR / "person_level_general_exposures.tsv", exposures)
    locus_exposures = derive_locus_exposures(load_catalog_rows(), samples)
    write_tsv(OUTDIR / "person_level_artifact_affected_locus_exposures.tsv", locus_exposures)

    comparison = []
    for feature in GENERAL_FEATURES:
        old_values = old.set_index("sample")[feature].reindex(samples).astype(float)
        new_values = exposures.set_index("sample")[feature].reindex(samples).astype(float)
        delta = new_values - old_values
        comparison.append(
            {
                "exposure_id": feature,
                "changed_people": int((delta.abs() > 1e-12).sum()),
                "max_absolute_delta": float(delta.abs().max()),
                "mean_delta": float(delta.mean()),
                "old_min": float(old_values.min()),
                "old_max": float(old_values.max()),
                "new_min": float(new_values.min()),
                "new_max": float(new_values.max()),
            }
        )
    write_tsv(OUTDIR / "old_vs_artifact_filtered_exposure_audit.tsv", pd.DataFrame(comparison))

    locus_comparison = []
    for feature in LOCUS_FEATURES:
        old_values = old.set_index("sample")[feature].reindex(samples).astype(float)
        new_values = (
            locus_exposures.set_index("sample")[feature].reindex(samples).astype(float)
        )
        changed = ~(
            (old_values.isna() & new_values.isna())
            | (old_values.fillna(-999999) == new_values.fillna(-999999))
        )
        locus_comparison.append(
            {
                "exposure_id": feature,
                "changed_people": int(changed.sum()),
                "changed_people_with_Mandage_outcome": int(
                    (
                        changed
                        & old.set_index("sample")["log2_ebv_load"]
                        .reindex(samples)
                        .notna()
                    ).sum()
                ),
            }
        )
    write_tsv(
        OUTDIR / "artifact_affected_locus_exposure_audit.tsv",
        pd.DataFrame(locus_comparison),
    )

    base = old.drop(
        columns=[
            column
            for column in [*GENERAL_FEATURES, *LOCUS_FEATURES]
            if column in old
        ]
    )
    matrix = base.merge(exposures, on="sample", how="left", validate="one_to_one")
    matrix = matrix.merge(
        locus_exposures, on="sample", how="left", validate="one_to_one"
    )
    outcomes: dict[str, pd.DataFrame] = {}
    mandage = matrix[matrix["log2_ebv_load"].notna()].copy()
    mandage["outcome"] = mandage["log2_ebv_load"]
    outcomes["Mandage2017_EBV_in_silico"] = mandage

    hold = pd.read_csv(HOLD, sep="\t", na_values=["NA"])
    hold["outcome"] = np.log2(hold["relative_ebv_copy_number_qpcr_mean"])
    outcomes["Houldcroft2014_EBV_qPCR"] = matrix.merge(
        hold[["sample", "outcome"]], on="sample", how="inner"
    )
    im = pd.read_csv(IM, sep="\t", na_values=["NA"])
    im["outcome"] = im["intrinsic_growth_rate"] / 10000.0
    outcomes["Im2012_intrinsic_growth"] = matrix.merge(
        im[["sample", "outcome"]], on="sample", how="inner"
    )

    model_rows = []
    for outcome_name, data in outcomes.items():
        for feature in GENERAL_FEATURES:
            if feature.startswith("general::type1_compatible_"):
                adjustment = "general::type1_physical_units"
            elif feature.startswith("general::type2_compatible_"):
                adjustment = "general::type2_physical_units"
            else:
                adjustment = None
            primary = cluster_fit(data, feature, adjustment, "superpopulation")
            sensitivity = cluster_fit(data, feature, adjustment, "population")
            model_rows.append(
                {
                    "outcome": outcome_name,
                    "exposure_id": feature,
                    "adjustment": adjustment or "",
                    "model_status": primary["status"],
                    "n": primary.get("n", 0),
                    "beta": primary.get("beta", np.nan),
                    "se_pedigree_cluster": primary.get("se_cluster", np.nan),
                    "p_pedigree_cluster": primary.get("p_cluster", np.nan),
                    "ci95_low": primary.get("ci95_low", np.nan),
                    "ci95_high": primary.get("ci95_high", np.nan),
                    "standardized_beta": primary.get("standardized_beta", np.nan),
                    "population_model_status": sensitivity["status"],
                    "population_beta": sensitivity.get("beta", np.nan),
                    "population_p_cluster": sensitivity.get("p_cluster", np.nan),
                }
            )
    models = pd.DataFrame(model_rows)
    models["q_bh_51_model_suite"] = bh(models["p_pedigree_cluster"])
    write_tsv(
        OUTDIR / "artifact_filtered_general_burden_model_results.tsv",
        models.sort_values(["q_bh_51_model_suite", "p_pedigree_cluster"]),
    )

    locus_model_rows = []
    for outcome_name, data in outcomes.items():
        # Replace the outcome frames' stale locus columns with the corrected vectors.
        data = data.drop(
            columns=[column for column in LOCUS_FEATURES if column in data],
            errors="ignore",
        ).merge(locus_exposures, on="sample", how="left", validate="one_to_one")
        for feature in LOCUS_FEATURES:
            adjustment = (
                "HML-2_1q22::orf::orf_gag::provirus_dosage_adjustment"
                if feature == "HML-2_1q22::orf::orf_gag"
                and "HML-2_1q22::orf::orf_gag::provirus_dosage_adjustment"
                in data
                else None
            )
            known = data[feature].dropna()
            zero = int((known == 0).sum())
            nonzero = int((known != 0).sum())
            if len(known) < 30:
                primary = {"status": "underpowered_n_lt30", "n": len(known)}
            elif min(zero, nonzero) < 5:
                primary = {"status": "underpowered_arm_lt5", "n": len(known)}
            else:
                primary = cluster_fit(data, feature, adjustment, "superpopulation")
            locus_model_rows.append(
                {
                    "outcome": outcome_name,
                    "exposure_id": feature,
                    "adjustment": adjustment or "",
                    "model_status": primary["status"],
                    "n": primary.get("n", 0),
                    "beta": primary.get("beta", np.nan),
                    "se_pedigree_cluster": primary.get("se_cluster", np.nan),
                    "p_pedigree_cluster": primary.get("p_cluster", np.nan),
                    "ci95_low": primary.get("ci95_low", np.nan),
                    "ci95_high": primary.get("ci95_high", np.nan),
                }
            )
    locus_models = pd.DataFrame(locus_model_rows)
    locus_models["q_bh_12_model_audit"] = bh(
        locus_models["p_pedigree_cluster"]
    )
    write_tsv(
        OUTDIR / "artifact_filtered_affected_locus_model_results.tsv",
        locus_models.sort_values(["q_bh_12_model_audit", "p_pedigree_cluster"]),
    )

    # Reinsert refitted exposures into the 79-feature x three-outcome grid.
    # BH uses only finite P values. Unestimable models retain missing P and q.
    # Unlike the earlier source-specific correction, this combined family
    # retains exact-vector aliases as separate exposure-outcome rows.
    old_mandage = pd.read_csv(OLD_MANDAGE_MODELS, sep="\t", na_values=["NA"])
    old_mandage = old_mandage.assign(
        outcome="Mandage2017_EBV_in_silico",
        old_p_pedigree_cluster=old_mandage["p_pedigree_cluster"],
    )
    old_secondary = pd.read_csv(
        OLD_SECONDARY_MODELS, sep="\t", na_values=["NA"]
    ).assign(old_p_pedigree_cluster=lambda frame: frame["p_pedigree_cluster"])
    full = pd.concat(
        [
            old_mandage[
                ["outcome", "exposure_id", "old_p_pedigree_cluster", "model_status"]
            ],
            old_secondary[
                ["outcome", "exposure_id", "old_p_pedigree_cluster", "model_status"]
            ],
        ],
        ignore_index=True,
    )
    replacements = pd.concat(
        [
            models[
                [
                    "outcome",
                    "exposure_id",
                    "model_status",
                    "p_pedigree_cluster",
                    "beta",
                    "n",
                ]
            ],
            locus_models[
                [
                    "outcome",
                    "exposure_id",
                    "model_status",
                    "p_pedigree_cluster",
                    "beta",
                    "n",
                ]
            ],
        ],
        ignore_index=True,
    ).rename(
        columns={
            "p_pedigree_cluster": "refit_p_pedigree_cluster",
            "model_status": "refit_model_status",
            "beta": "refit_beta",
            "n": "refit_n",
        }
    )
    full = full.merge(
        replacements,
        on=["outcome", "exposure_id"],
        how="left",
        validate="one_to_one",
    )
    full["was_refit"] = full["refit_model_status"].notna()
    full["analysis_p_pedigree_cluster"] = np.where(
        full["was_refit"],
        full["refit_p_pedigree_cluster"],
        full["old_p_pedigree_cluster"],
    )
    full["analysis_model_status"] = full["refit_model_status"].fillna(full["model_status"])
    full = full.drop(columns=["model_status"])
    full["included_in_bh_family"] = np.isfinite(full["analysis_p_pedigree_cluster"])
    full["bh_family_size"] = int(full["included_in_bh_family"].sum())
    full["bh_exclusion_reason"] = np.where(
        full["included_in_bh_family"], "", full["analysis_model_status"]
    )
    full["q_bh_finite_model_suite"] = bh(
        full["analysis_p_pedigree_cluster"]
    )
    write_tsv(
        OUTDIR / "complete_three_outcome_refit_multiplicity.tsv",
        full.sort_values(
            [
                "q_bh_finite_model_suite",
                "analysis_p_pedigree_cluster",
            ]
        ),
    )

    excluded = pd.read_csv(EXCLUSIONS, sep="\t", na_values=["NA"])
    artifact_samples = sorted(set(excluded["ID"]))
    overlap = []
    for outcome_name, data in outcomes.items():
        cohort = set(data["sample"])
        overlap.append(
            {
                "outcome": outcome_name,
                "cohort_n": len(cohort),
                "artifact_sample_overlap_n": len(cohort.intersection(artifact_samples)),
                "artifact_sample_overlap": ";".join(
                    sorted(cohort.intersection(artifact_samples))
                ),
            }
        )
    write_tsv(OUTDIR / "artifact_sample_outcome_overlap.tsv", pd.DataFrame(overlap))

    summary = {
        "schema": "hml2.biologically-filtered-functional-refit.v2",
        "analysis_catalog": str(CATALOG),
        "analysis_rows": 59_656,
        "excluded_rows": {
            "non_HML2_HML11_sequence_identity": 1168,
            "assembly_artifact_not_supported_by_CNV_depth": 35,
            "alias_duplicate_of_8q24.3c": 584,
            "duplicate_catalog_label_for_same_assembled_interval": 78,
        },
        "catalog_rows_used": int(len(prov)),
        "people": len(samples),
        "features_refit": len(GENERAL_FEATURES),
        "outcomes_refit": list(outcomes),
        "models": len(models),
        "artifact_affected_locus_models": len(locus_models),
        "pol_rule": (
            "compatible Pol units require compatible gag, pro, and pol calls on "
            "the same proviral copy; pol is not counted as a standalone product"
        ),
        "changed_exposures": int(
            sum(row["changed_people"] > 0 for row in comparison)
        ),
        "nominal_p_lt_0_05": int(
            (models["p_pedigree_cluster"] < 0.05).fillna(False).sum()
        ),
        "suite_q_lt_0_05": int(
            (models["q_bh_51_model_suite"] < 0.05).fillna(False).sum()
        ),
        "combined_attempted_models": len(full),
        "combined_finite_p_values": int(full["included_in_bh_family"].sum()),
        "combined_missing_p_values": int((~full["included_in_bh_family"]).sum()),
        "combined_family_deduplicates_exposure_aliases": False,
        "combined_finite_model_q_lt_0_05": int(
            (full["q_bh_finite_model_suite"] < 0.05)
            .fillna(False)
            .sum()
        ),
    }
    (OUTDIR / "run_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
