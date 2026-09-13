#!/usr/bin/env python3
"""Read-only copy-state likelihood authority for retained HML-2 CNV summaries.

The module deliberately separates calibration rows from candidate rows.  A
candidate CNV row is reduced to one region-support statistic and enters the
state likelihood exactly once.  No prior is applied: reported state weights
are normalized likelihoods, not posterior probabilities.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence


STATE_ARTIFACT = "ASSEMBLY_ARTIFACT_OR_SPURIOUS_CONTIG"
STATE_TRUE = "TRUE_LATER_DUPLICATION"


def _number(value: object, field: str) -> float:
    if value is None or str(value).strip() in {"", "NA", "N/A", "nan"}:
        raise ValueError(f"{field} is missing")
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{field} must be a finite non-negative number")
    return result


def _quantile(values: Sequence[float], probability: float) -> float:
    """Linear-interpolated quantile, matching the common type-7 definition."""
    if not values:
        raise ValueError("quantile requires at least one value")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * probability
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def region_metrics(
    row: Mapping[str, object], *, body_crossmap_ambiguous: bool
) -> dict[str, object]:
    """Create the single assay statistic used by the state likelihood.

    ``sample_ref_cov`` is the retained sample-level diploid target-flank
    baseline.  Therefore a single-copy body expectation is ref/2, while a
    represented local flank expectation is ref.  A non-ambiguous region must
    support both body and local flank, so its support is their conservative
    minimum.  If body reads can cross-map from an identical copy, body is
    retained for diagnostics but excluded from the likelihood statistic.
    Values above the one-copy expectation saturate at one because this test is
    presence-versus-artifact, not a copy-number magnitude estimator.
    """
    body = _number(row.get("body_median_all"), "body_median_all")
    flank = _number(row.get("flank_median_all"), "flank_median_all")
    reference = _number(row.get("sample_ref_cov"), "sample_ref_cov")
    if reference == 0:
        raise ValueError("sample_ref_cov must be positive")

    expected_haploid = reference / 2.0
    body_ratio = body / expected_haploid
    flank_ratio = flank / reference
    if body_crossmap_ambiguous:
        support_raw = flank_ratio
        likelihood_features = ["local_flank_ratio"]
    else:
        support_raw = min(body_ratio, flank_ratio)
        likelihood_features = ["min(body_ratio,local_flank_ratio)"]
    support = min(1.0, max(0.0, support_raw))

    if body == 0 and flank == 0:
        concordance = None
    elif body_ratio == 0 or flank_ratio == 0:
        concordance = math.inf
    else:
        concordance = abs(math.log2(body_ratio / flank_ratio))

    return {
        "body_median_all": body,
        "local_flank_median_all": flank,
        "sample_diploid_target_flank_baseline": reference,
        "expected_haploid_body_coverage": expected_haploid,
        "body_ratio_to_one_copy": body_ratio,
        "local_flank_ratio_to_represented_contig": flank_ratio,
        "absolute_log2_body_flank_concordance": concordance,
        "body_crossmap_ambiguous": body_crossmap_ambiguous,
        "likelihood_features": likelihood_features,
        "region_support_ratio_uncapped": max(0.0, support_raw),
        "region_support_ratio": support,
    }


@dataclass(frozen=True)
class Calibration:
    n_positive_controls: int
    sigma_narrow: float
    sigma_central: float
    sigma_wide: float
    control_keys: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "n_positive_controls": self.n_positive_controls,
            "state_centers": {STATE_ARTIFACT: 0.0, STATE_TRUE: 1.0},
            "scales": {
                "narrow_q75_lower_residual": self.sigma_narrow,
                "central_rmse_lower_residual": self.sigma_central,
                "wide_q90_lower_residual": self.sigma_wide,
            },
            "control_keys": list(self.control_keys),
            "prior": None,
            "weight_semantics": "normalized_state_likelihood_not_posterior",
        }


def calibrate_positive_controls(
    controls: Iterable[tuple[str, Mapping[str, object]]], *,
    body_crossmap_ambiguous: bool,
) -> Calibration:
    """Calibrate sensitivity scales from independently declared real tips.

    Callers must supply controls whose physical reality is external to this CNV
    assay (authenticated segmental duplications, explicit spreadsheet plus
    shared-LTR arrays, or ordinary singleton tips).  Candidate keys may not be
    included.  The cross-map setting is explicit and shared across controls.
    """
    keys: list[str] = []
    residuals: list[float] = []
    for key, row in controls:
        if key in keys:
            raise ValueError(f"duplicate calibration key: {key}")
        support = float(region_metrics(
            row, body_crossmap_ambiguous=body_crossmap_ambiguous
        )["region_support_ratio_uncapped"])
        keys.append(key)
        # Only undercoverage can mimic a missing/spurious contig.  Coverage
        # above the represented-copy expectation is therefore a zero lower-
        # tail residual, not symmetric error that broadens the artifact model.
        residuals.append(max(0.0, 1.0 - support))
    if len(residuals) < 3:
        raise ValueError("at least three external positive controls are required")

    sigma_narrow = _quantile(residuals, 0.75)
    sigma_central = math.sqrt(sum(value * value for value in residuals) / len(residuals))
    sigma_wide = _quantile(residuals, 0.90)
    positive_scales = [value for value in (sigma_narrow, sigma_central, sigma_wide) if value > 0]
    if not positive_scales:
        raise ValueError("positive controls have zero residual spread; likelihood scale is unidentified")
    # Preserve named estimators while preventing a zero-width likelihood.
    replacement = min(positive_scales)
    sigma_narrow = sigma_narrow or replacement
    sigma_central = sigma_central or replacement
    sigma_wide = sigma_wide or replacement
    return Calibration(
        len(keys), sigma_narrow, sigma_central, sigma_wide, tuple(keys)
    )


def _normalized_likelihood(support: float, sigma: float) -> dict[str, float]:
    if sigma <= 0 or not math.isfinite(sigma):
        raise ValueError("sigma must be finite and positive")
    log_likelihood = {
        STATE_ARTIFACT: -0.5 * ((support - 0.0) / sigma) ** 2,
        STATE_TRUE: -0.5 * ((support - 1.0) / sigma) ** 2,
    }
    anchor = max(log_likelihood.values())
    likelihood = {key: math.exp(value - anchor) for key, value in log_likelihood.items()}
    total = sum(likelihood.values())
    return {key: value / total for key, value in likelihood.items()}


def score_candidate(
    *,
    candidate_key: str,
    row: Mapping[str, object],
    calibration: Calibration,
    always_biological_tips: Sequence[str],
    candidate_tip: str,
    terminal_host: str,
    body_crossmap_ambiguous: bool,
) -> dict[str, object]:
    """Return state likelihoods and state-dependent biological membership."""
    if candidate_key in calibration.control_keys:
        raise ValueError("candidate row must not also be a calibration row")
    if candidate_tip in always_biological_tips:
        raise ValueError("candidate_tip must be state-dependent, not always present")
    metrics = region_metrics(row, body_crossmap_ambiguous=body_crossmap_ambiguous)
    support = float(metrics["region_support_ratio"])
    scales = {
        "narrow_q75_lower_residual": calibration.sigma_narrow,
        "central_rmse_lower_residual": calibration.sigma_central,
        "wide_q90_lower_residual": calibration.sigma_wide,
    }
    likelihoods = {name: _normalized_likelihood(support, sigma) for name, sigma in scales.items()}
    central = likelihoods["central_rmse_lower_residual"]
    favored = max(central, key=central.get)
    return {
        "candidate_key": candidate_key,
        "assay_use_count": 1,
        "metrics": metrics,
        "state_likelihoods": likelihoods,
        "central_state_likelihoods": central,
        "central_likelihood_favors": favored,
        "state_contract": {
            STATE_ARTIFACT: {
                "biological_tips": list(always_biological_tips),
                "excluded_candidate_measurements": [candidate_tip],
                "terminal_hosts": [terminal_host],
            },
            STATE_TRUE: {
                "biological_tips": list(always_biological_tips) + [candidate_tip],
                "excluded_candidate_measurements": [],
                "terminal_hosts": [terminal_host],
            },
        },
        "consumer_rule": (
            "Marginalize state-conditional biological results exactly once with central_state_likelihoods; "
            "do not reuse CNV metrics as a prior, likelihood, filter, or downstream weight."
        ),
    }
