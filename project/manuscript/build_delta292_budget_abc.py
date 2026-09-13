#!/usr/bin/env python3
"""Mechanistic budget ABC for the success of the HML-2 Delta292 state.

This extends the component-gain sensitivity analysis by parameterizing the
candidate mechanisms as measurable budgets:

* RNA redirection: Type-II transcript budget diverted from transmissible
  genomic RNA, and the fraction recovered by Delta292.
* Producer-cost avoidance: limiting producer output lost to deleted products,
  and the fraction of that cost avoided by Delta292.
* Packaging: gain per genomic RNA molecule.

The count likelihood identifies the total gain.  Posterior contribution shares
show how the explicit mechanism priors divide that total and which empirical
measurements can collapse the remaining ridge.  No "unknown" model is used.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript/delta292_budget_abc_v1"
AUTONOMY = ROOT / "manuscript/delta292_autonomy_hypothesis_audit_v1/summary.json"
SOURCE = ROOT / "manuscript/delta292_source_vs_effect_v1/model_comparison.tsv"

SEED = 292_2026_07_27 + 91
N_DRAWS = 500_000
TYPE2_LENGTH = 9472
TYPE1_LENGTH = 9180
G_LENGTH = TYPE2_LENGTH / TYPE1_LENGTH


@dataclass(frozen=True)
class Prior:
    name: str
    moi_low: float
    moi_high: float
    ascertainment_low: float
    ascertainment_high: float
    max_diverted_rna_fraction: float
    max_producer_cost_fraction: float
    max_packaging_gain: float


PRIORS = (
    Prior("conservative", 1.0, 5.0, 0.67, 1.5, 0.35, 0.35, 1.10),
    Prior("moderate", 1.0, 5.0, 0.67, 1.5, 0.50, 0.50, 1.20),
    Prior("broad", 0.5, 10.0, 0.25, 4.0, 0.70, 0.70, 1.35),
)

MODELS = {
    "length_only": (),
    "rna_redirection": ("rna",),
    "producer_cost_avoidance": ("producer",),
    "packaging_per_genomic_rna": ("packaging",),
    "rna_plus_producer": ("rna", "producer"),
    "rna_plus_packaging": ("rna", "packaging"),
    "producer_plus_packaging": ("producer", "packaging"),
    "all_three": ("rna", "producer", "packaging"),
}

OBSERVATIONS = (
    ("orthology_aware_floor", 10, 16),
    ("inclusive_ceiling", 17, 23),
)

PARAMETERS = (
    "moi",
    "ascertainment_odds",
    "rna_diverted_fraction",
    "rna_recovery_fraction",
    "rna_gain",
    "producer_cost_fraction",
    "producer_cost_avoidance_fraction",
    "producer_gain",
    "packaging_gain",
    "effective_gain",
    "latent_fraction",
    "observed_fraction",
    "rna_log_gain_share",
    "producer_log_gain_share",
    "packaging_log_gain_share",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def log_uniform(
    rng: np.random.Generator, low: float, high: float, size: int
) -> np.ndarray:
    return np.exp(rng.uniform(math.log(low), math.log(high), size))


def weighted_quantile(
    values: np.ndarray, weights: np.ndarray, probabilities: tuple[float, ...]
) -> np.ndarray:
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    if cumulative[-1] <= 0:
        return np.full(len(probabilities), np.nan)
    cumulative /= cumulative[-1]
    return np.interp(probabilities, cumulative, values)


def binomial_likelihood(k: int, n: int, p: np.ndarray) -> np.ndarray:
    p = np.clip(p, 1e-12, 1.0 - 1e-12)
    return math.comb(n, k) * p**k * (1.0 - p) ** (n - k)


def equilibrium_fraction(gain: np.ndarray, moi: np.ndarray) -> np.ndarray:
    fraction = np.zeros_like(gain)
    viable = gain > 1.0
    fraction[viable] = 1.0 + np.log1p(-1.0 / gain[viable]) / moi[viable]
    return np.clip(fraction, 0.0, 1.0)


def observed_probability(
    latent: np.ndarray, ascertainment_odds: np.ndarray
) -> np.ndarray:
    denominator = 1.0 - latent + ascertainment_odds * latent
    return np.divide(
        ascertainment_odds * latent,
        denominator,
        out=np.zeros_like(latent),
        where=denominator > 0,
    )


def draw_model(
    prior: Prior, model: str, seed: int
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = N_DRAWS
    active = MODELS[model]
    moi = log_uniform(rng, prior.moi_low, prior.moi_high, n)
    ascertainment = log_uniform(
        rng, prior.ascertainment_low, prior.ascertainment_high, n
    )

    rna_diverted = np.zeros(n)
    rna_recovery = np.zeros(n)
    rna_gain = np.ones(n)
    if "rna" in active:
        # Beta(2, 3) avoids making either zero or the upper bound the default.
        rna_diverted = (
            prior.max_diverted_rna_fraction * rng.beta(2.0, 3.0, n)
        )
        rna_recovery = rng.beta(2.0, 2.0, n)
        rna_gain = (
            1.0
            - rna_diverted
            + rna_recovery * rna_diverted
        ) / (1.0 - rna_diverted)

    producer_cost = np.zeros(n)
    producer_avoidance = np.zeros(n)
    producer_gain = np.ones(n)
    if "producer" in active:
        producer_cost = (
            prior.max_producer_cost_fraction * rng.beta(2.0, 3.0, n)
        )
        producer_avoidance = rng.beta(2.0, 2.0, n)
        producer_gain = 1.0 / (
            1.0 - producer_cost * producer_avoidance
        )

    packaging_gain = np.ones(n)
    if "packaging" in active:
        packaging_gain = log_uniform(rng, 1.0, prior.max_packaging_gain, n)

    effective_gain = G_LENGTH * rna_gain * producer_gain * packaging_gain
    latent = equilibrium_fraction(effective_gain, moi)
    observed = observed_probability(latent, ascertainment)

    nonlength_log = np.log(rna_gain * producer_gain * packaging_gain)
    denominator = np.where(nonlength_log > 0, nonlength_log, 1.0)
    return {
        "moi": moi,
        "ascertainment_odds": ascertainment,
        "rna_diverted_fraction": rna_diverted,
        "rna_recovery_fraction": rna_recovery,
        "rna_gain": rna_gain,
        "producer_cost_fraction": producer_cost,
        "producer_cost_avoidance_fraction": producer_avoidance,
        "producer_gain": producer_gain,
        "packaging_gain": packaging_gain,
        "effective_gain": effective_gain,
        "latent_fraction": latent,
        "observed_fraction": observed,
        "rna_log_gain_share": np.log(rna_gain) / denominator,
        "producer_log_gain_share": np.log(producer_gain) / denominator,
        "packaging_log_gain_share": np.log(packaging_gain) / denominator,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    comparison_rows: list[dict[str, object]] = []
    posterior_rows: list[dict[str, object]] = []
    dominant_rows: list[dict[str, object]] = []
    summary_results: dict[str, object] = {}

    for prior_i, prior in enumerate(PRIORS):
        for obs_i, (observation, k, n) in enumerate(OBSERVATIONS):
            evidence: dict[str, float] = {}
            simulations: dict[str, dict[str, np.ndarray]] = {}
            weights: dict[str, np.ndarray] = {}
            for model_i, model in enumerate(MODELS):
                sim = draw_model(
                    prior,
                    model,
                    SEED + prior_i * 100_000 + obs_i * 10_000 + model_i,
                )
                weight = binomial_likelihood(k, n, sim["observed_fraction"])
                simulations[model] = sim
                weights[model] = weight
                evidence[model] = float(weight.mean())

            total = sum(evidence.values())
            key = f"{prior.name}::{observation}"
            summary_results[key] = {}
            for model in MODELS:
                model_share = evidence[model] / total
                comparison_rows.append(
                    {
                        "prior": prior.name,
                        "observation": observation,
                        "observed_typeI_units": k,
                        "total_units": n,
                        "model": model,
                        "draws": N_DRAWS,
                        "mean_binomial_likelihood": evidence[model],
                        "likelihood_ratio_vs_length_only": (
                            evidence[model] / evidence["length_only"]
                        ),
                        "equal_model_prior_share": model_share,
                        "moi_low": prior.moi_low,
                        "moi_high": prior.moi_high,
                        "max_diverted_rna_fraction": (
                            prior.max_diverted_rna_fraction
                        ),
                        "max_producer_cost_fraction": (
                            prior.max_producer_cost_fraction
                        ),
                        "max_packaging_gain": prior.max_packaging_gain,
                    }
                )
                normalized = weights[model] / weights[model].sum()
                for parameter in PARAMETERS:
                    q025, q25, median, q75, q975 = weighted_quantile(
                        simulations[model][parameter],
                        normalized,
                        (0.025, 0.25, 0.5, 0.75, 0.975),
                    )
                    posterior_rows.append(
                        {
                            "prior": prior.name,
                            "observation": observation,
                            "model": model,
                            "parameter": parameter,
                            "q025": q025,
                            "q25": q25,
                            "median": median,
                            "q75": q75,
                            "q975": q975,
                        }
                    )

                sim = simulations[model]
                for component, field in (
                    ("rna", "rna_log_gain_share"),
                    ("producer", "producer_log_gain_share"),
                    ("packaging", "packaging_log_gain_share"),
                ):
                    dominant_rows.append(
                        {
                            "prior": prior.name,
                            "observation": observation,
                            "model": model,
                            "component": component,
                            "posterior_probability_component_supplies_majority_of_nonlength_log_gain": float(
                                np.sum(
                                    normalized
                                    * (sim[field] > 0.5)
                                )
                            ),
                            "posterior_probability_component_supplies_at_least_quarter_of_nonlength_log_gain": float(
                                np.sum(
                                    normalized
                                    * (sim[field] >= 0.25)
                                )
                            ),
                        }
                    )
                summary_results[key][model] = {
                    "likelihood_ratio_vs_length_only": (
                        evidence[model] / evidence["length_only"]
                    ),
                    "equal_model_prior_share": model_share,
                }

    write_tsv(OUT / "budget_model_comparison.tsv", comparison_rows)
    write_tsv(OUT / "posterior_budget_summary.tsv", posterior_rows)
    write_tsv(OUT / "posterior_component_dominance.tsv", dominant_rows)

    summary = {
        "analysis": "delta292_budget_abc_v1",
        "draws_per_prior_model_observation": N_DRAWS,
        "fixed_length_gain": G_LENGTH,
        "models": {name: list(active) for name, active in MODELS.items()},
        "results": summary_results,
        "empirical_branch_updates": {
            "retained_upstream_autonomy": (
                "not supported by current artifact-filtered locus-level "
                "Gag-Pro-Pol sequence compatibility"
            ),
            "RNA_redirection": (
                "active; locus-resolved full-genome NCCIT remap is the direct "
                "likelihood intended to replace the provisional transcript proxy"
            ),
            "producer_cost_avoidance": (
                "active; posterior cost thresholds are reported rather than "
                "converted into an unsupported protein claim"
            ),
            "packaging_per_genomic_RNA": (
                "active; posterior gain thresholds are reported; known core "
                "packaging intervals are retained outside Delta292"
            ),
        },
        "claim_boundary": (
            "The count likelihood identifies a total effective gain and tests "
            "whether each explicitly parameterized budget can supply it. "
            "Relative molecular contributions remain prior-sensitive until a "
            "component-specific assay is added, but the analysis returns "
            "quantitative posterior targets rather than stopping at that ridge."
        ),
        "input_sha256": {
            str(AUTONOMY.relative_to(ROOT)): sha256(AUTONOMY),
            str(SOURCE.relative_to(ROOT)): sha256(SOURCE),
        },
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    (OUT / "MODEL_SPECIFICATION.md").write_text(
        "# Δ292 mechanistic budget ABC\n\n"
        "For the RNA branch, a Type-II diverted fraction `f` and recovery "
        "fraction `r` produce `G_RNA=(1-f+r*f)/(1-f)`. For the producer branch, "
        "cost fraction `c` and avoided fraction `a` produce "
        "`G_producer=1/(1-c*a)`. Packaging is a direct gain per genomic RNA. "
        "These multiply with the fixed 9472/9180 length term in the same "
        "frequency-dependent helper equilibrium used by the component ABC.\n\n"
        "The same prevalence likelihood cannot by itself tell RNA redirection "
        "from producer-cost avoidance when their priors induce the same gain "
        "distribution. That equality is a testable structural result, not an "
        "endpoint: `posterior_budget_summary.tsv` gives the RNA, producer, and "
        "packaging measurements that would move the posterior. The independent "
        "artifact-filtered Gag-Pro-Pol audit is carried alongside this model and "
        "removes retained extant upstream autonomy as positive support.\n"
    )
    print(OUT)


if __name__ == "__main__":
    main()
