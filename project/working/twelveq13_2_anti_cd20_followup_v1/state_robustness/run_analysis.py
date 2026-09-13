#!/usr/bin/env python3
"""Robust observational follow-up of the 12q13.2 anti-CD20 association."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RESULTS = HERE / "results"

ANTI = ROOT / "project/working/direct_anti_cd20_viability_screen_v1"
OUTCOMES = ANTI / "results/direct_outcomes_qc.tsv"
STRUCTURE = ANTI / "results/twelveq13_identity_semantics.tsv"
EXISTING_FOCAL = ANTI / "results/focal_results.tsv"
EXISTING_STRUCTURAL = ANTI / "results/twelveq13_posthoc_structural_contrasts.tsv"
DIRECT = ROOT / "project/working/locus_marker_expansion_v1/association_ready_matrix.tsv"
PERSON_MARKERS = ROOT / "project/working/locus_marker_expansion_v1/person_functional_markers_long.tsv"
MARKER_CATALOG = ROOT / "project/working/locus_marker_expansion_v1/marker_catalog.tsv"

OUTCOME_LABELS = {
    "Obin_media": "Obinutuzumab + media",
    "Obin_serum": "Obinutuzumab + serum/complement",
    "Ofat_media": "Ofatumumab + media",
    "Ofat_serum": "Ofatumumab + serum/complement",
    "Ritux_serum": "Rituximab + serum/complement",
}

MARKERS = {
    "internal": "12q13.2::internal_presence::all_people",
    "retained": "12q13.2::retained_element_presence::all_people",
    "provirus": "12q13.2::provirus_presence::all_people",
    "dose": "12q13.2::internal_haplotype_dose::all_people",
    "gag_orf": "12q13.2::gag_orf::extended",
    "pro_orf": "12q13.2::pro_orf::strict_phase0",
    "env_orf": "12q13.2::env_orf::extended",
    "ordered_route": "12q13.2::ordered_gag_pro_pol_route::corrected_clean_sequence",
    "np9": "12q13.2::np9_sequence::extended",
    "rec": "12q13.2::rec_sequence::extended",
}

MODEL_IDS = [
    "internal_vs_pooled_noninternal",
    "any_retained_element_vs_true_absence",
    "provirus_vs_true_absence_excluding_retained_noninternal",
    "internal_haplotype_dose_linear",
    "pro_orf_positive_vs_known_negative",
    "structural_internal_vs_retained_noninternal",
    "structural_internal_vs_true_absence",
    "structural_retained_noninternal_vs_true_absence",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def seed_for(text: str) -> int:
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:4], "little")


def adjust_bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    keep = values.notna()
    if not keep.any():
        return result
    p = values.loc[keep].to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    m = len(p)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty(m, dtype=float)
    restored[order] = np.minimum(adjusted, 1.0)
    result.loc[keep] = restored
    return result


def fit_contrast(y: np.ndarray, design: np.ndarray, contrast: np.ndarray) -> dict[str, float | int]:
    n, _ = design.shape
    rank = int(np.linalg.matrix_rank(design))
    df = n - rank
    if df <= 0:
        raise ValueError("nonpositive residual degrees of freedom")
    xtx_inv = np.linalg.pinv(design.T @ design)
    beta_all = xtx_inv @ design.T @ y
    residual = y - design @ beta_all
    estimate = float(contrast @ beta_all)
    sigma2 = float(residual @ residual / df)
    classical_var = float(contrast @ xtx_inv @ contrast)
    classical_se = float(np.sqrt(max(sigma2 * classical_var, 0.0)))
    classical_t = estimate / classical_se if classical_se > 0 else np.nan
    classical_p = float(2 * stats.t.sf(abs(classical_t), df)) if np.isfinite(classical_t) else np.nan

    leverage = np.sum((design @ xtx_inv) * design, axis=1)
    score = design * (residual / np.maximum(1.0 - leverage, 1e-10))[:, None]
    hc3_cov = xtx_inv @ (score.T @ score) @ xtx_inv
    hc3_var = float(contrast @ hc3_cov @ contrast)
    hc3_se = float(np.sqrt(max(hc3_var, 0.0)))
    hc3_t = estimate / hc3_se if hc3_se > 0 else np.nan
    hc3_p = float(2 * stats.t.sf(abs(hc3_t), df)) if np.isfinite(hc3_t) else np.nan
    total = float(np.sum((y - y.mean()) ** 2))
    return {
        "n": n,
        "rank": rank,
        "residual_df": df,
        "condition_number": float(np.linalg.cond(design)),
        "beta": estimate,
        "classical_se": classical_se,
        "classical_t": classical_t,
        "classical_p": classical_p,
        "hc3_se": hc3_se,
        "hc3_t": hc3_t,
        "hc3_p": hc3_p,
        "r_squared": float(1.0 - (residual @ residual) / total) if total > 0 else np.nan,
    }


def model_frame(frame: pd.DataFrame, model_id: str) -> tuple[pd.DataFrame, list[str], np.ndarray]:
    data = frame.copy()
    if model_id == "internal_vs_pooled_noninternal":
        data["x_internal"] = data["internal"]
        return data, ["x_internal"], np.array([1.0])
    if model_id == "any_retained_element_vs_true_absence":
        data["x_retained"] = data["retained"]
        return data, ["x_retained"], np.array([1.0])
    if model_id == "provirus_vs_true_absence_excluding_retained_noninternal":
        data = data[data["provirus"].notna()].copy()
        data["x_provirus"] = data["provirus"]
        return data, ["x_provirus"], np.array([1.0])
    if model_id == "internal_haplotype_dose_linear":
        data["x_dose"] = data["dose"]
        return data, ["x_dose"], np.array([1.0])
    if model_id == "pro_orf_positive_vs_known_negative":
        data = data[data["pro_orf_exposure"].notna()].copy()
        data["x_pro_orf"] = data["pro_orf_exposure"]
        return data, ["x_pro_orf"], np.array([1.0])
    data["x_internal_class"] = (data["structural_class"] == "internal_element_present").astype(float)
    data["x_absent_class"] = (data["structural_class"] == "biological_locus_absent").astype(float)
    if model_id == "structural_internal_vs_retained_noninternal":
        contrast = np.array([1.0, 0.0])
    elif model_id == "structural_internal_vs_true_absence":
        contrast = np.array([1.0, -1.0])
    elif model_id == "structural_retained_noninternal_vs_true_absence":
        contrast = np.array([0.0, -1.0])
    else:
        raise KeyError(model_id)
    return data, ["x_internal_class", "x_absent_class"], contrast


def build_design(
    frame: pd.DataFrame, exposure_columns: list[str], contrast_exposure: np.ndarray, ancestry: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    covariates = pd.DataFrame(index=frame.index)
    covariates["sex_male"] = frame["sex"].map({"female": 0.0, "male": 1.0})
    dummies = pd.get_dummies(frame[ancestry], prefix=ancestry, drop_first=True, dtype=float)
    covariates = pd.concat([covariates, dummies], axis=1)
    exposure = frame[exposure_columns].astype(float).to_numpy()
    design = np.column_stack([np.ones(len(frame)), exposure, covariates.to_numpy(dtype=float)])
    contrast = np.zeros(design.shape[1], dtype=float)
    contrast[1 : 1 + len(exposure_columns)] = contrast_exposure
    reduced = np.column_stack([np.ones(len(frame)), covariates.to_numpy(dtype=float)])
    return design, contrast, reduced, len(exposure_columns)


def arm_status(data: pd.DataFrame, exposure_columns: list[str], minimum_arm: int = 5) -> tuple[str, str]:
    if len(exposure_columns) == 2:
        counts = data["structural_class"].value_counts()
        expected = ["internal_element_present", "retained_noninternal_reference", "biological_locus_absent"]
        if any(counts.get(state, 0) < minimum_arm for state in expected):
            return "blocked_minimum_structural_class", json.dumps(counts.to_dict(), sort_keys=True)
        return "fitted", json.dumps(counts.to_dict(), sort_keys=True)
    x = data[exposure_columns[0]].dropna()
    counts = x.value_counts().sort_index()
    if x.nunique() < 2:
        return "blocked_nonvariable", json.dumps({str(k): int(v) for k, v in counts.items()})
    if set(x.unique()).issubset({0.0, 1.0}) and counts.min() < minimum_arm:
        return "blocked_minimum_binary_arm", json.dumps({str(k): int(v) for k, v in counts.items()})
    if not set(x.unique()).issubset({0.0, 1.0}) and x.nunique() < 3:
        return "blocked_insufficient_dose_support", json.dumps({str(k): int(v) for k, v in counts.items()})
    return "fitted", json.dumps({str(k): int(v) for k, v in counts.items()})


def permute_within_groups(n: int, groups: pd.Series, n_resamples: int, rng: np.random.Generator) -> np.ndarray:
    indices = np.tile(np.arange(n), (n_resamples, 1))
    group_values = groups.to_numpy()
    for label in pd.unique(group_values):
        positions = np.flatnonzero(group_values == label)
        if len(positions) <= 1:
            continue
        random_order = np.argsort(rng.random((n_resamples, len(positions))), axis=1)
        indices[:, positions] = positions[random_order]
    return indices


def freedman_lane_contrast(
    y: np.ndarray,
    full: np.ndarray,
    reduced: np.ndarray,
    contrast: np.ndarray,
    groups: pd.Series,
    seed: int,
    n_resamples: int = 20000,
) -> dict[str, float | int | str]:
    reduced_beta = np.linalg.pinv(reduced) @ y
    fitted = reduced @ reduced_beta
    residual = y - fitted
    permutations = permute_within_groups(len(y), groups, n_resamples, np.random.default_rng(seed))
    y_star = fitted[:, None] + residual[permutations].T

    full_pinv = np.linalg.pinv(full)
    beta = full_pinv @ y_star
    estimates = contrast @ beta
    residual_star = y_star - full @ beta
    df = len(y) - np.linalg.matrix_rank(full)
    base_var = float(contrast @ np.linalg.pinv(full.T @ full) @ contrast)
    se = np.sqrt(np.maximum(np.sum(residual_star * residual_star, axis=0) / df * base_var, 0.0))
    t_star = np.divide(estimates, se, out=np.zeros_like(estimates), where=se > 0)

    observed = fit_contrast(y, full, contrast)["classical_t"]
    extreme = int(np.sum(np.abs(t_star) >= abs(float(observed)) - 1e-14))
    p = (extreme + 1) / (n_resamples + 1)
    return {
        "permutation_method": "Freedman-Lane residual permutation within 1000 Genomes population; classical studentized contrast",
        "permutation_resamples": n_resamples,
        "permutation_seed": seed,
        "permutation_p": p,
        "permutation_mc_se": float(np.sqrt(p * (1 - p) / (n_resamples + 1))),
    }


def freedman_lane_omnibus(
    y: np.ndarray,
    full: np.ndarray,
    reduced: np.ndarray,
    groups: pd.Series,
    seed: int,
    n_resamples: int = 20000,
) -> dict[str, float | int | str]:
    full_rank = int(np.linalg.matrix_rank(full))
    reduced_rank = int(np.linalg.matrix_rank(reduced))
    q = full_rank - reduced_rank
    df = len(y) - full_rank
    full_residual = y - full @ (np.linalg.pinv(full) @ y)
    reduced_residual = y - reduced @ (np.linalg.pinv(reduced) @ y)
    rss_full = float(full_residual @ full_residual)
    rss_reduced = float(reduced_residual @ reduced_residual)
    observed_f = ((rss_reduced - rss_full) / q) / (rss_full / df)
    parametric_p = float(stats.f.sf(observed_f, q, df))

    fitted = y - reduced_residual
    permutations = permute_within_groups(len(y), groups, n_resamples, np.random.default_rng(seed))
    y_star = fitted[:, None] + reduced_residual[permutations].T
    full_beta = np.linalg.pinv(full) @ y_star
    reduced_beta = np.linalg.pinv(reduced) @ y_star
    rss_full_star = np.sum((y_star - full @ full_beta) ** 2, axis=0)
    rss_reduced_star = np.sum((y_star - reduced @ reduced_beta) ** 2, axis=0)
    f_star = ((rss_reduced_star - rss_full_star) / q) / (rss_full_star / df)
    extreme = int(np.sum(f_star >= observed_f - 1e-14))
    p = (extreme + 1) / (n_resamples + 1)
    return {
        "partial_f": observed_f,
        "partial_f_df1": q,
        "partial_f_df2": df,
        "parametric_partial_f_p": parametric_p,
        "permutation_method": "Freedman-Lane residual permutation within 1000 Genomes population; partial F",
        "permutation_resamples": n_resamples,
        "permutation_seed": seed,
        "permutation_p": p,
        "permutation_mc_se": float(np.sqrt(p * (1 - p) / (n_resamples + 1))),
    }


def permutation_mean_difference(
    a: np.ndarray,
    b: np.ndarray,
    seed: int,
    exact_threshold: int = 1_000_000,
    mc_resamples: int = 100_000,
) -> dict[str, float | int | str]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    y = np.concatenate([a, b])
    n_a = len(a)
    observed = float(a.mean() - b.mean())
    total_allocations = math.comb(len(y), n_a)
    tolerance = 1e-14
    if total_allocations <= exact_threshold:
        total_sum = float(y.sum())
        extreme = 0
        for selected in itertools.combinations(range(len(y)), n_a):
            sum_a = float(y[list(selected)].sum())
            diff = sum_a / n_a - (total_sum - sum_a) / (len(y) - n_a)
            extreme += abs(diff) >= abs(observed) - tolerance
        p = extreme / total_allocations
        return {
            "permutation_mean_difference": observed,
            "permutation_p": p,
            "permutation_method": "exact_all_label_allocations",
            "permutation_allocations_or_resamples": total_allocations,
            "permutation_seed": np.nan,
            "permutation_mc_se": 0.0,
        }

    rng = np.random.default_rng(seed)
    extreme = 0
    completed = 0
    batch_size = 5000
    total_sum = float(y.sum())
    while completed < mc_resamples:
        batch = min(batch_size, mc_resamples - completed)
        selected = np.argpartition(rng.random((batch, len(y))), n_a - 1, axis=1)[:, :n_a]
        sum_a = y[selected].sum(axis=1)
        diff = sum_a / n_a - (total_sum - sum_a) / (len(y) - n_a)
        extreme += int(np.sum(np.abs(diff) >= abs(observed) - tolerance))
        completed += batch
    p = (extreme + 1) / (mc_resamples + 1)
    return {
        "permutation_mean_difference": observed,
        "permutation_p": p,
        "permutation_method": "Monte_Carlo_unrestricted_label_permutation",
        "permutation_allocations_or_resamples": mc_resamples,
        "permutation_seed": seed,
        "permutation_mc_se": float(np.sqrt(p * (1 - p) / (mc_resamples + 1))),
    }


def spearman_permutation(
    x: np.ndarray, y: np.ndarray, seed: int, n_resamples: int = 100_000
) -> dict[str, float | int | str]:
    x_rank = stats.rankdata(x).astype(float)
    y_rank = stats.rankdata(y).astype(float)
    x_center = x_rank - x_rank.mean()
    y_center = y_rank - y_rank.mean()
    denominator = float(np.sqrt(np.sum(x_center**2) * np.sum(y_center**2)))
    observed = float(np.dot(x_center, y_center) / denominator)
    rng = np.random.default_rng(seed)
    extreme = 0
    completed = 0
    batch_size = 5000
    while completed < n_resamples:
        batch = min(batch_size, n_resamples - completed)
        order = np.argsort(rng.random((batch, len(y))), axis=1)
        permuted = y_center[order]
        rho = permuted @ x_center / denominator
        extreme += int(np.sum(np.abs(rho) >= abs(observed) - 1e-14))
        completed += batch
    p = (extreme + 1) / (n_resamples + 1)
    return {
        "spearman_rho": observed,
        "permutation_p": p,
        "permutation_method": "Monte_Carlo_unrestricted_rank_permutation",
        "permutation_resamples": n_resamples,
        "permutation_seed": seed,
        "permutation_mc_se": float(np.sqrt(p * (1 - p) / (n_resamples + 1))),
    }


def summarize_numeric(values: pd.Series) -> dict[str, float | int]:
    values = values.dropna().astype(float)
    return {
        "n_outcome_observed": int(len(values)),
        "mean": float(values.mean()) if len(values) else np.nan,
        "sd": float(values.std(ddof=1)) if len(values) > 1 else np.nan,
        "median": float(values.median()) if len(values) else np.nan,
        "q1": float(values.quantile(0.25)) if len(values) else np.nan,
        "q3": float(values.quantile(0.75)) if len(values) else np.nan,
        "minimum": float(values.min()) if len(values) else np.nan,
        "maximum": float(values.max()) if len(values) else np.nan,
    }


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    required = [OUTCOMES, STRUCTURE, EXISTING_FOCAL, EXISTING_STRUCTURAL, DIRECT, PERSON_MARKERS, MARKER_CATALOG]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("missing required sources: " + ", ".join(missing))

    outcomes = pd.read_csv(OUTCOMES, sep="\t").set_index("sample")
    structure = pd.read_csv(STRUCTURE, sep="\t").set_index("sample")
    direct = pd.read_csv(DIRECT, sep="\t").set_index("sample")
    sample_ids = structure.index.tolist()
    if len(sample_ids) != 49 or set(sample_ids) != set(outcomes.index):
        raise ValueError("expected the same 49 individuals in structural and outcome authorities")

    people = structure.loc[sample_ids].copy()
    for short in ["internal", "retained", "provirus", "dose"]:
        authoritative = pd.to_numeric(direct.loc[sample_ids, MARKERS[short]], errors="coerce")
        observed = pd.to_numeric(people[short], errors="coerce")
        if not np.allclose(authoritative.fillna(-999), observed.fillna(-999)):
            raise ValueError(f"12q13.2 {short} differs between direct and structural authority")

    marker_long = pd.read_csv(PERSON_MARKERS, sep="\t")
    marker_long = marker_long[
        (marker_long["locus"] == "12q13.2")
        & marker_long["sample"].isin(sample_ids)
        & marker_long["marker_id"].isin(MARKERS.values())
    ].copy()
    if marker_long.duplicated(["sample", "marker_id"]).any():
        raise ValueError("duplicate accepted person-marker state")
    marker_state = marker_long.pivot(index="sample", columns="marker_id", values="person_state")
    marker_exposure = marker_long.pivot(index="sample", columns="marker_id", values="exposure_value")
    for short in ["gag_orf", "pro_orf", "env_orf", "ordered_route", "np9", "rec"]:
        marker_id = MARKERS[short]
        people[f"{short}_state"] = marker_state.reindex(sample_ids)[marker_id]
        people[f"{short}_exposure"] = pd.to_numeric(marker_exposure.reindex(sample_ids)[marker_id], errors="coerce")

    people["provirus_contrast_inclusion"] = np.where(
        people["provirus"].notna(), "included_provirus_or_true_absence", "excluded_retained_noninternal_structure"
    )
    for short in ["gag_orf", "pro_orf", "env_orf", "ordered_route", "np9", "rec"]:
        people[f"{short}_analysis_inclusion"] = np.where(
            people[f"{short}_exposure"].notna(), "included_callable_internal", "excluded_by_structural_or_callable_state"
        )
    people = people.join(outcomes.loc[sample_ids], how="left")
    for outcome_id in OUTCOME_LABELS:
        people[f"{outcome_id}_missing_reason"] = np.where(
            people[outcome_id].isna(), "source_endpoint_QC_missing", "observed"
        )
    people.index.name = "sample"
    people.reset_index().to_csv(RESULTS / "all_individual_rows.tsv", sep="\t", index=False, na_rep="")

    # Accepted sequence-state catalog, including descriptive states that cannot
    # support a test in this overlap.
    catalog = pd.read_csv(MARKER_CATALOG, sep="\t").set_index("marker_id")
    catalog_rows = []
    for short in ["gag_orf", "pro_orf", "env_orf", "ordered_route", "np9", "rec"]:
        marker_id = MARKERS[short]
        states = people[f"{short}_state"]
        exposure = people[f"{short}_exposure"]
        counts = exposure.dropna().value_counts().sort_index()
        if exposure.notna().sum() == 0:
            testability = "blocked_no_callable_numeric_contrast"
        elif exposure.dropna().nunique() < 2:
            testability = "blocked_nonvariable"
        elif counts.min() < 5:
            testability = "blocked_minimum_binary_arm"
        else:
            testability = "testable"
        info = catalog.loc[marker_id]
        catalog_rows.append(
            {
                "state_short_name": short,
                "marker_id": marker_id,
                "endpoint": info["endpoint"],
                "functional_claim_boundary": info["functional_claim_boundary"],
                "n_total": len(people),
                "n_numeric_callable": int(exposure.notna().sum()),
                "n_structural_or_callable_excluded": int(exposure.isna().sum()),
                "n_unique_numeric": int(exposure.dropna().nunique()),
                "minimum_numeric_arm": int(counts.min()) if len(counts) > 1 else np.nan,
                "person_state_counts_json": json.dumps(states.value_counts(dropna=False).to_dict(), sort_keys=True),
                "exposure_counts_json": json.dumps({str(k): int(v) for k, v in counts.items()}, sort_keys=True),
                "testability": testability,
            }
        )
    pd.DataFrame(catalog_rows).to_csv(RESULTS / "accepted_orf_state_catalog.tsv", sep="\t", index=False, na_rep="")

    # Raw distributions preserve excluded states instead of silently dropping
    # them from the descriptive audit.
    state_variables = {
        "structural_class": people["structural_class"],
        "internal_presence": people["internal"].map({0: "no_internal_element", 1: "internal_element_present"}),
        "retained_element_presence": people["retained"].map({0: "true_biological_absence", 1: "retained_element_present"}),
        "provirus_contrast": pd.Series(
            np.where(
                people["provirus"].isna(),
                "excluded_retained_noninternal_structure",
                people["provirus"].map({0.0: "true_biological_absence", 1.0: "provirus_present"}),
            ),
            index=people.index,
        ),
        "internal_haplotype_dose": people["dose"].map(lambda value: f"dose_{int(value)}"),
        "gag_orf_accepted_state": people["gag_orf_state"],
        "pro_orf_accepted_state": people["pro_orf_state"],
        "env_orf_accepted_state": people["env_orf_state"],
        "ordered_gag_pro_pol_route_state": people["ordered_route_state"],
    }
    distribution_rows = []
    for outcome_id, outcome_label in OUTCOME_LABELS.items():
        for variable, state_series in state_variables.items():
            for state, indices in state_series.groupby(state_series, dropna=False).groups.items():
                values = people.loc[list(indices), outcome_id]
                distribution_rows.append(
                    {
                        "outcome_id": outcome_id,
                        "outcome_label": outcome_label,
                        "state_variable": variable,
                        "state": state,
                        "n_people_in_state": len(indices),
                        "n_outcome_missing": int(values.isna().sum()),
                        **summarize_numeric(values),
                    }
                )
    pd.DataFrame(distribution_rows).to_csv(RESULTS / "raw_state_distributions.tsv", sep="\t", index=False, na_rep="")

    # Unadjusted pairwise tests are distributional sensitivities. Exact label
    # enumeration is used when feasible; otherwise a deterministic Monte Carlo
    # permutation is reported with its simulation SE.
    pair_specs = {
        "internal_vs_pooled_noninternal": (
            people["structural_class"] == "internal_element_present",
            people["structural_class"] != "internal_element_present",
            "internal_element_present",
            "pooled_retained_noninternal_plus_true_absence",
        ),
        "structural_internal_vs_retained_noninternal": (
            people["structural_class"] == "internal_element_present",
            people["structural_class"] == "retained_noninternal_reference",
            "internal_element_present",
            "retained_noninternal_reference",
        ),
        "structural_internal_vs_true_absence": (
            people["structural_class"] == "internal_element_present",
            people["structural_class"] == "biological_locus_absent",
            "internal_element_present",
            "true_biological_absence",
        ),
        "structural_retained_noninternal_vs_true_absence": (
            people["structural_class"] == "retained_noninternal_reference",
            people["structural_class"] == "biological_locus_absent",
            "retained_noninternal_reference",
            "true_biological_absence",
        ),
        "any_retained_element_vs_true_absence": (
            people["retained"] == 1,
            people["retained"] == 0,
            "any_retained_element",
            "true_biological_absence",
        ),
        "pro_orf_positive_vs_known_negative": (
            people["pro_orf_exposure"] == 1,
            people["pro_orf_exposure"] == 0,
            "pro_orf_positive_carrier",
            "pro_orf_known_negative_only",
        ),
        "env_orf_positive_vs_known_negative": (
            people["env_orf_exposure"] == 1,
            people["env_orf_exposure"] == 0,
            "env_orf_positive_carrier",
            "env_orf_known_negative_only",
        ),
    }
    pair_rows = []
    for outcome_id, outcome_label in OUTCOME_LABELS.items():
        for comparison, (mask_a, mask_b, label_a, label_b) in pair_specs.items():
            a = people.loc[mask_a, outcome_id].dropna().to_numpy(dtype=float)
            b = people.loc[mask_b, outcome_id].dropna().to_numpy(dtype=float)
            if len(a) == 0 or len(b) == 0:
                continue
            welch = stats.ttest_ind(a, b, equal_var=False) if len(a) > 1 and len(b) > 1 else None
            mann = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
            permutation = permutation_mean_difference(a, b, seed_for(f"pair|{outcome_id}|{comparison}"))
            pair_rows.append(
                {
                    "outcome_id": outcome_id,
                    "outcome_label": outcome_label,
                    "comparison": comparison,
                    "group_a": label_a,
                    "group_b": label_b,
                    "n_a": len(a),
                    "n_b": len(b),
                    "mean_a": float(a.mean()),
                    "mean_b": float(b.mean()),
                    "median_a": float(np.median(a)),
                    "median_b": float(np.median(b)),
                    "mean_difference_a_minus_b": float(a.mean() - b.mean()),
                    "median_difference_a_minus_b": float(np.median(a) - np.median(b)),
                    "welch_t": float(welch.statistic) if welch is not None else np.nan,
                    "welch_p": float(welch.pvalue) if welch is not None else np.nan,
                    "mann_whitney_u": float(mann.statistic),
                    "mann_whitney_p_asymptotic_tie_corrected": float(mann.pvalue),
                    "rank_biserial_a_minus_b": float(2 * mann.statistic / (len(a) * len(b)) - 1),
                    "support_status": "inferential_sensitivity" if min(len(a), len(b)) >= 5 else "descriptive_low_support_arm_below_5",
                    **permutation,
                }
            )
    pairwise = pd.DataFrame(pair_rows)
    supported = pairwise["support_status"] == "inferential_sensitivity"
    pairwise["permutation_bh_supported_family"] = np.nan
    supported_p = pairwise.loc[supported, "permutation_p"]
    pairwise.loc[supported, "permutation_bh_supported_family"] = adjust_bh(supported_p)
    pairwise.to_csv(RESULTS / "unadjusted_pairwise_tests.tsv", sep="\t", index=False, na_rep="")

    rank_rows = []
    for outcome_id, outcome_label in OUTCOME_LABELS.items():
        data = people[["dose", outcome_id]].dropna()
        spearman = stats.spearmanr(data["dose"], data[outcome_id])
        kendall = stats.kendalltau(data["dose"], data[outcome_id])
        permutation = spearman_permutation(
            data["dose"].to_numpy(dtype=float),
            data[outcome_id].to_numpy(dtype=float),
            seed_for(f"dose_rank|{outcome_id}"),
        )
        rank_rows.append(
            {
                "outcome_id": outcome_id,
                "outcome_label": outcome_label,
                "n": len(data),
                "spearman_parametric_p": float(spearman.pvalue),
                "kendall_tau_b": float(kendall.statistic),
                "kendall_p": float(kendall.pvalue),
                **permutation,
            }
        )
    rank_tests = pd.DataFrame(rank_rows)
    rank_tests["permutation_bh_five_outcomes"] = adjust_bh(rank_tests["permutation_p"])
    rank_tests.to_csv(RESULTS / "dose_rank_trend_tests.tsv", sep="\t", index=False, na_rep="")

    # Covariate-adjusted models plus superpopulation and permutation
    # sensitivities. Positive beta always means greater viability/resistance.
    adjusted_rows = []
    status_rows = []
    fitted_context: dict[tuple[str, str], tuple[pd.DataFrame, list[str], np.ndarray]] = {}
    for outcome_id, outcome_label in OUTCOME_LABELS.items():
        for model_id in MODEL_IDS:
            base = people[people[outcome_id].notna()].copy()
            data, exposure_columns, contrast_exposure = model_frame(base, model_id)
            status, counts_json = arm_status(data, exposure_columns)
            status_row = {
                "outcome_id": outcome_id,
                "model_id": model_id,
                "n": len(data),
                "state_counts_json": counts_json,
                "status": status,
            }
            if status != "fitted":
                status_rows.append(status_row)
                continue
            full, contrast, reduced, _ = build_design(data, exposure_columns, contrast_exposure, "population")
            residual_df = len(data) - np.linalg.matrix_rank(full)
            condition = float(np.linalg.cond(full))
            if residual_df < 10:
                status_row["status"] = "blocked_residual_df_below_10"
                status_rows.append(status_row)
                continue
            if not np.isfinite(condition) or condition > 1e8:
                status_row["status"] = "blocked_ill_conditioned"
                status_rows.append(status_row)
                continue
            status_rows.append(status_row)
            y = data[outcome_id].to_numpy(dtype=float)
            primary = fit_contrast(y, full, contrast)
            super_full, super_contrast, _, _ = build_design(data, exposure_columns, contrast_exposure, "superpopulation")
            sensitivity = fit_contrast(y, super_full, super_contrast)
            permutation = freedman_lane_contrast(
                y,
                full,
                reduced,
                contrast,
                data["population"],
                seed_for(f"adjusted|{outcome_id}|{model_id}"),
            )
            adjusted_rows.append(
                {
                    "outcome_id": outcome_id,
                    "outcome_label": outcome_label,
                    "model_id": model_id,
                    "positive_beta_interpretation": "greater viability/resistance to antibody-mediated killing",
                    "analysis_role": "same_cohort_observational_followup_not_independent_replication",
                    **primary,
                    "superpopulation_fe_beta": sensitivity["beta"],
                    "superpopulation_fe_hc3_se": sensitivity["hc3_se"],
                    "superpopulation_fe_hc3_p": sensitivity["hc3_p"],
                    **permutation,
                }
            )
            fitted_context[(outcome_id, model_id)] = (data, exposure_columns, contrast_exposure)
    adjusted = pd.DataFrame(adjusted_rows)
    adjusted["hc3_bh_all_adjusted_state_models"] = adjust_bh(adjusted["hc3_p"])
    adjusted["permutation_bh_all_adjusted_state_models"] = adjust_bh(adjusted["permutation_p"])
    adjusted.to_csv(RESULTS / "adjusted_state_models.tsv", sep="\t", index=False, na_rep="")
    pd.DataFrame(status_rows).to_csv(RESULTS / "adjusted_model_status.tsv", sep="\t", index=False, na_rep="")

    # Three-class omnibus test.
    omnibus_rows = []
    for outcome_id, outcome_label in OUTCOME_LABELS.items():
        data = people[people[outcome_id].notna()].copy()
        data, exposure_columns, contrast_exposure = model_frame(data, "structural_internal_vs_retained_noninternal")
        full, _, reduced, _ = build_design(data, exposure_columns, contrast_exposure, "population")
        kruskal_groups = [
            group[outcome_id].to_numpy(dtype=float)
            for _, group in data.groupby("structural_class")
            if len(group)
        ]
        kruskal = stats.kruskal(*kruskal_groups)
        permutation = freedman_lane_omnibus(
            data[outcome_id].to_numpy(dtype=float),
            full,
            reduced,
            data["population"],
            seed_for(f"omnibus|{outcome_id}"),
        )
        omnibus_rows.append(
            {
                "outcome_id": outcome_id,
                "outcome_label": outcome_label,
                "n": len(data),
                "kruskal_wallis_h": float(kruskal.statistic),
                "kruskal_wallis_p": float(kruskal.pvalue),
                **permutation,
            }
        )
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["permutation_bh_five_outcomes"] = adjust_bh(omnibus["permutation_p"])
    omnibus.to_csv(RESULTS / "structural_omnibus_tests.tsv", sep="\t", index=False, na_rep="")

    # Individual and population leave-one-out. These are influence audits, not
    # new hypothesis tests; all estimates are retained even if P changes.
    individual_rows = []
    population_rows = []
    for adjusted_row in adjusted.itertuples(index=False):
        outcome_id = adjusted_row.outcome_id
        model_id = adjusted_row.model_id
        data, exposure_columns, contrast_exposure = fitted_context[(outcome_id, model_id)]
        for omitted_sample in data.index:
            retained_data = data.drop(index=omitted_sample)
            loo_status, counts_json = arm_status(retained_data, exposure_columns, minimum_arm=1)
            row = {
                "outcome_id": outcome_id,
                "model_id": model_id,
                "omitted_sample": omitted_sample,
                "omitted_population": data.loc[omitted_sample, "population"],
                "n": len(retained_data),
                "state_counts_json": counts_json,
                "status": loo_status,
            }
            if loo_status == "fitted":
                full, contrast, _, _ = build_design(retained_data, exposure_columns, contrast_exposure, "population")
                if len(retained_data) - np.linalg.matrix_rank(full) > 0 and np.linalg.cond(full) <= 1e8:
                    result = fit_contrast(retained_data[outcome_id].to_numpy(dtype=float), full, contrast)
                    row.update({"beta": result["beta"], "hc3_se": result["hc3_se"], "hc3_p": result["hc3_p"]})
                else:
                    row["status"] = "blocked_rank_or_condition"
            individual_rows.append(row)

        for omitted_population in sorted(data["population"].unique()):
            retained_data = data[data["population"] != omitted_population].copy()
            loo_status, counts_json = arm_status(retained_data, exposure_columns, minimum_arm=1)
            row = {
                "outcome_id": outcome_id,
                "model_id": model_id,
                "omitted_population": omitted_population,
                "n": len(retained_data),
                "state_counts_json": counts_json,
                "status": loo_status,
            }
            if loo_status == "fitted":
                full, contrast, _, _ = build_design(retained_data, exposure_columns, contrast_exposure, "population")
                if len(retained_data) - np.linalg.matrix_rank(full) > 0 and np.linalg.cond(full) <= 1e8:
                    result = fit_contrast(retained_data[outcome_id].to_numpy(dtype=float), full, contrast)
                    row.update({"beta": result["beta"], "hc3_se": result["hc3_se"], "hc3_p": result["hc3_p"]})
                else:
                    row["status"] = "blocked_rank_or_condition"
            population_rows.append(row)

    individual = pd.DataFrame(individual_rows)
    individual.to_csv(RESULTS / "individual_leave_one_out.tsv", sep="\t", index=False, na_rep="")
    population_loo = pd.DataFrame(population_rows)
    population_loo.to_csv(RESULTS / "population_leave_one_out.tsv", sep="\t", index=False, na_rep="")

    summary_rows = []
    for (outcome_id, model_id), group in individual.groupby(["outcome_id", "model_id"]):
        fitted = group[group["status"] == "fitted"].copy()
        full_beta = float(adjusted[(adjusted["outcome_id"] == outcome_id) & (adjusted["model_id"] == model_id)]["beta"].iloc[0])
        if len(fitted):
            fitted["absolute_beta_change"] = (fitted["beta"] - full_beta).abs()
            influential = fitted.sort_values("absolute_beta_change", ascending=False).iloc[0]
            summary_rows.append(
                {
                    "outcome_id": outcome_id,
                    "model_id": model_id,
                    "full_beta": full_beta,
                    "n_successful_leave_one_out": len(fitted),
                    "n_blocked_leave_one_out": len(group) - len(fitted),
                    "beta_min": fitted["beta"].min(),
                    "beta_max": fitted["beta"].max(),
                    "sign_agreement_fraction": float((np.sign(fitted["beta"]) == np.sign(full_beta)).mean()),
                    "hc3_p_min": fitted["hc3_p"].min(),
                    "hc3_p_max": fitted["hc3_p"].max(),
                    "maximum_absolute_beta_change": influential["absolute_beta_change"],
                    "most_influential_omitted_sample": influential["omitted_sample"],
                    "most_influential_omitted_population": influential["omitted_population"],
                }
            )
    pd.DataFrame(summary_rows).to_csv(RESULTS / "individual_leave_one_out_summary.tsv", sep="\t", index=False, na_rep="")

    population_summary_rows = []
    for (outcome_id, model_id), group in population_loo.groupby(["outcome_id", "model_id"]):
        fitted = group[group["status"] == "fitted"]
        full_beta = float(adjusted[(adjusted["outcome_id"] == outcome_id) & (adjusted["model_id"] == model_id)]["beta"].iloc[0])
        population_summary_rows.append(
            {
                "outcome_id": outcome_id,
                "model_id": model_id,
                "full_beta": full_beta,
                "n_successful_population_omissions": len(fitted),
                "n_blocked_population_omissions": len(group) - len(fitted),
                "beta_min": fitted["beta"].min() if len(fitted) else np.nan,
                "beta_max": fitted["beta"].max() if len(fitted) else np.nan,
                "sign_agreement_fraction": float((np.sign(fitted["beta"]) == np.sign(full_beta)).mean()) if len(fitted) else np.nan,
                "hc3_p_min": fitted["hc3_p"].min() if len(fitted) else np.nan,
                "hc3_p_max": fitted["hc3_p"].max() if len(fitted) else np.nan,
            }
        )
    pd.DataFrame(population_summary_rows).to_csv(
        RESULTS / "population_leave_one_out_summary.tsv", sep="\t", index=False, na_rep=""
    )

    source_rows = [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path), "role": role}
        for path, role in [
            (OUTCOMES, "QC-adjusted five anti-CD20 outcomes for the 49 direct overlaps"),
            (STRUCTURE, "authoritative 12q13.2 retained/internal/provirus/dose decomposition"),
            (DIRECT, "complete direct person-locus marker grid"),
            (PERSON_MARKERS, "accepted per-person functional/ORF state authority"),
            (MARKER_CATALOG, "marker definitions and functional claim boundaries"),
            (EXISTING_FOCAL, "cross-check against primary adjusted internal-presence result"),
            (EXISTING_STRUCTURAL, "cross-check against existing three-class structural contrasts"),
        ]
    ]
    pd.DataFrame(source_rows).to_csv(RESULTS / "source_manifest.tsv", sep="\t", index=False)

    summary = {
        "n_people": len(people),
        "structural_counts": people["structural_class"].value_counts().to_dict(),
        "dose_counts": {str(int(k)): int(v) for k, v in people["dose"].value_counts().sort_index().items()},
        "pro_orf_counts": people["pro_orf_state"].value_counts().to_dict(),
        "outcome_nonmissing": {outcome: int(people[outcome].notna().sum()) for outcome in OUTCOME_LABELS},
        "adjusted_models_fitted": len(adjusted),
        "permutation_resamples_adjusted_per_cell": 20000,
        "unadjusted_pairwise_mc_resamples": 100000,
        "interpretation": "exploratory same-cohort observational follow-up; positive beta is greater viability/resistance",
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
