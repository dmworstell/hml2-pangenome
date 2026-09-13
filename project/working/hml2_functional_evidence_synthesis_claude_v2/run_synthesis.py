#!/usr/bin/env python
"""
Standardized-effect synthesis of the donor-disjoint MAGE / GEUVADIS
1p31.1b internal-fragment -> SLC44A5 expression results.

Lane T1-S.  Exclusive write path:
    project/working/hml2_functional_evidence_synthesis_claude_v2/results/synthesis_*

Everything upstream is READ-ONLY.  This module does not refit either cohort:
it verifies the frozen upstream coefficients against the real result files,
harmonizes them onto one standardized scale, and pools them.

Run:
    python run_synthesis.py
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from dataclasses import dataclass, asdict, field

import numpy as np
import pandas as pd
from scipy import stats

ROOT = str(Path(__file__).resolve().parents[3])
HERE = os.path.join(ROOT, "project/working/hml2_functional_evidence_synthesis_claude_v2")
OUT = os.path.join(HERE, "results")

# The one biological endpoint class that may be pooled here.  Two sources may be
# combined only if this string is identical for both.
ENDPOINT_STEADY_STATE_SLC44A5_LCL = "steady_state_SLC44A5_mRNA_abundance_in_EBV_LCL"

TOL = 1e-9


class EndpointClassMismatch(RuntimeError):
    """Raised when a pool is attempted across different biological endpoint classes."""


class SourceVerificationError(RuntimeError):
    """Raised when a pinned upstream value disagrees with the real result file."""


# --------------------------------------------------------------------------
# Source records
# --------------------------------------------------------------------------

@dataclass
class Source:
    cohort: str
    n: int
    endpoint_class: str
    endpoint_detail: str
    exposure_coding: str
    expression_normalization: str
    covariates: str
    se_kind: str
    beta_raw: float
    se_raw: float
    p_raw: float
    resid_df: int
    n_covariates: int          # nuisance covariate columns, excluding intercept and exposure
    sd_x: float
    sd_y: float
    correction_family: str
    source_file: str
    harmonization: str = ""
    beta_std: float = field(default=float("nan"))
    se_std: float = field(default=float("nan"))


# --------------------------------------------------------------------------
# Standardization algebra
# --------------------------------------------------------------------------

def standardize(beta_raw: float, se_raw: float, sd_x: float, sd_y: float):
    """Map a raw regression slope onto the fully standardized scale.

    Convention (identical for every source in this synthesis):

        beta_std = beta_raw * sd_x / sd_y
        se_std   = se_raw   * sd_x / sd_y

    Units of beta_std: SD of expression per SD of exposure.

    sd_x and sd_y are the marginal within-analysis-sample standard deviations of
    the exposure and of the expression outcome.  They are treated as fixed
    constants, so this is a pure linear rescale by k = sd_x / sd_y.  The Wald
    statistic beta/se is therefore invariant, which is what makes the transform
    safe: it moves the point estimate and its uncertainty together and cannot
    manufacture or destroy significance.

    This is exactly the convention already used by the GEUVADIS package
    (`standardized_beta = beta * sd(x) / y_sd` in its run_analysis.R), applied
    unchanged to MAGE.
    """
    if not (sd_x > 0 and sd_y > 0):
        raise ValueError("sd_x and sd_y must be strictly positive")
    k = sd_x / sd_y
    return beta_raw * k, se_raw * k


def unstandardize(beta_std: float, se_std: float, sd_x: float, sd_y: float):
    """Exact inverse of `standardize`.  beta_raw = beta_std * sd_y / sd_x."""
    if not (sd_x > 0 and sd_y > 0):
        raise ValueError("sd_x and sd_y must be strictly positive")
    k = sd_y / sd_x
    return beta_std * k, se_std * k


def partial_correlation(t_stat: float, resid_df: int) -> float:
    """Partial correlation implied by a Wald t and its residual df.

    r = t / sqrt(t^2 + df).  Invariant to any linear rescaling of x or y, so it
    is immune to the expression-normalization asymmetry between the two cohorts.
    """
    return t_stat / math.sqrt(t_stat * t_stat + resid_df)


def fisher_z(r: float) -> float:
    return math.atanh(r)


def fisher_z_se(n: int, n_covariates: int) -> float:
    """SE of Fisher's z for a partial correlation controlling `n_covariates`."""
    denom = n - n_covariates - 3
    if denom <= 0:
        raise ValueError("not enough residual information for a Fisher z SE")
    return 1.0 / math.sqrt(denom)


# --------------------------------------------------------------------------
# Endpoint-class guard
# --------------------------------------------------------------------------

def assert_poolable(sources) -> str:
    """Refuse to pool biologically different endpoints.

    Pooling is only meaningful when every source measures the same biological
    quantity.  This guard is deliberately strict: it raises rather than warns.
    """
    if len(sources) < 2:
        raise ValueError("need at least two sources to pool")
    classes = {s.endpoint_class for s in sources}
    if len(classes) != 1:
        raise EndpointClassMismatch(
            "refusing to pool different biological endpoint classes: "
            + " | ".join(sorted(classes))
        )
    return classes.pop()


# --------------------------------------------------------------------------
# Meta-analysis
# --------------------------------------------------------------------------

def fixed_effect(betas, ses):
    """Inverse-variance fixed-effect pool."""
    b = np.asarray(betas, dtype=float)
    s = np.asarray(ses, dtype=float)
    if np.any(s <= 0):
        raise ValueError("standard errors must be positive")
    w = 1.0 / s**2
    est = float(np.sum(w * b) / np.sum(w))
    se = float(math.sqrt(1.0 / np.sum(w)))
    z = est / se
    p = float(2.0 * stats.norm.sf(abs(z)))
    return {
        "estimate": est,
        "se": se,
        "ci95_low": est - 1.959963984540054 * se,
        "ci95_high": est + 1.959963984540054 * se,
        "z": float(z),
        "p": p,
        "weights": w.tolist(),
        "weight_fraction": (w / np.sum(w)).tolist(),
    }


def cochran_q(betas, ses):
    """Cochran's Q about the fixed-effect mean, with df, P, and I^2."""
    b = np.asarray(betas, dtype=float)
    s = np.asarray(ses, dtype=float)
    w = 1.0 / s**2
    mu = np.sum(w * b) / np.sum(w)
    q = float(np.sum(w * (b - mu) ** 2))
    df = len(b) - 1
    p = float(stats.chi2.sf(q, df)) if df > 0 else float("nan")
    i2 = float(max(0.0, (q - df) / q) * 100.0) if q > 0 else 0.0
    return {"Q": q, "df": df, "p": p, "I2_percent": i2}


def tau2_dersimonian_laird(betas, ses):
    b = np.asarray(betas, dtype=float)
    s = np.asarray(ses, dtype=float)
    w = 1.0 / s**2
    mu = np.sum(w * b) / np.sum(w)
    q = np.sum(w * (b - mu) ** 2)
    df = len(b) - 1
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    if c <= 0:
        return 0.0
    return float(max(0.0, (q - df) / c))


def tau2_paule_mandel(betas, ses, max_iter=500, tol=1e-12):
    """Paule-Mandel tau^2 by bisection on the generalized-Q estimating equation.

    Solve for tau^2 >= 0 such that  sum_i (b_i - mu(tau^2))^2 / (se_i^2 + tau^2) = k - 1.
    The left side is strictly decreasing in tau^2, so bisection is safe.
    """
    b = np.asarray(betas, dtype=float)
    s = np.asarray(ses, dtype=float)
    k = len(b)
    target = k - 1
    if target <= 0:
        return 0.0

    def gen_q(t2):
        w = 1.0 / (s**2 + t2)
        mu = np.sum(w * b) / np.sum(w)
        return np.sum(w * (b - mu) ** 2)

    if gen_q(0.0) <= target:
        return 0.0
    lo, hi = 0.0, max(1.0, float(np.var(b, ddof=0)) * 10.0 + 1.0)
    it = 0
    while gen_q(hi) > target and it < 100:
        hi *= 2.0
        it += 1
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if gen_q(mid) > target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return float(0.5 * (lo + hi))


def random_effects(betas, ses, tau2):
    """Random-effects pool at a supplied tau^2 (DL or PM)."""
    b = np.asarray(betas, dtype=float)
    s = np.asarray(ses, dtype=float)
    w = 1.0 / (s**2 + tau2)
    est = float(np.sum(w * b) / np.sum(w))
    se = float(math.sqrt(1.0 / np.sum(w)))
    z = est / se
    p = float(2.0 * stats.norm.sf(abs(z)))
    return {
        "tau2": float(tau2),
        "tau": float(math.sqrt(tau2)),
        "estimate": est,
        "se": se,
        "ci95_low": est - 1.959963984540054 * se,
        "ci95_high": est + 1.959963984540054 * se,
        "z": float(z),
        "p": p,
        "weights": w.tolist(),
        "weight_fraction": (w / np.sum(w)).tolist(),
    }


def synthesize(sources, betas_attr="beta_std", ses_attr="se_std"):
    """Full synthesis for a poolable source list.  Enforces the endpoint guard."""
    endpoint = assert_poolable(sources)
    betas = [getattr(s, betas_attr) for s in sources]
    ses = [getattr(s, ses_attr) for s in sources]
    fe = fixed_effect(betas, ses)
    het = cochran_q(betas, ses)
    t2_dl = tau2_dersimonian_laird(betas, ses)
    t2_pm = tau2_paule_mandel(betas, ses)
    return {
        "endpoint_class": endpoint,
        "sources": [s.cohort for s in sources],
        "fixed_effect": fe,
        "random_effects_dersimonian_laird": random_effects(betas, ses, t2_dl),
        "random_effects_paule_mandel": random_effects(betas, ses, t2_pm),
        "heterogeneity": het,
    }


# --------------------------------------------------------------------------
# Upstream verification
# --------------------------------------------------------------------------

def _close(a, b, tol=1e-9):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


def independent_disjointness_derivation():
    """Derive the GEUVADIS/MAGE donor intersection from primary ID tables.

    This deliberately does NOT read `mage_overlap_audit.tsv`.  That file is the
    GEUVADIS package's own *assertion* of the overlap; using it to verify the
    overlap would be circular.  Instead we intersect:

      - the MAGE analysis roster, from the MAGE package's sample_alignment.tsv
        (the 39 corrected biological identities actually fitted), and
      - the GEUVADIS analysis roster, from the GEUVADIS package's
        direct33_exposures.tsv (the 33 donors actually fitted).

    The donor-disjointness premise of the primary pool stands or falls here.
    """
    mage_ids = sorted(pd.read_csv(
        os.path.join(ROOT, "project/working/onep31b_slc44a5_followup_v1/robust_stats/results/sample_alignment.tsv"),
        sep="\t")["sample"].astype(str))
    geu_ids = sorted(pd.read_csv(
        os.path.join(ROOT, "project/working/direct_geuvadis_function_screen_v1/data/direct33_exposures.tsv"),
        sep="\t")["sample"].astype(str))

    intersection = sorted(set(mage_ids) & set(geu_ids))
    geu_disjoint = sorted(set(geu_ids) - set(intersection))
    residual = sorted(set(geu_disjoint) & set(mage_ids))

    # Cross-check the package's asserted overlap against our own derivation.
    asserted_file = os.path.join(
        ROOT, "project/working/direct_geuvadis_function_screen_v1/data/mage_overlap_audit.tsv")
    asserted = sorted(pd.read_csv(asserted_file, sep="\t").query("in_prior_MAGE_direct39")["sample"].astype(str))

    return {
        "method": "direct intersection of the two packages' primary analysis rosters; "
                  "mage_overlap_audit.tsv deliberately NOT used as input",
        "n_mage": len(mage_ids),
        "n_geuvadis_total": len(geu_ids),
        "independently_derived_shared_donors": intersection,
        "n_shared": len(intersection),
        "n_geuvadis_after_removal": len(geu_disjoint),
        "intersection_after_removal": residual,
        "disjoint": len(residual) == 0,
        "asserted_shared_donors": asserted,
        "independent_derivation_agrees_with_assertion": intersection == asserted,
        "status": ("VERIFIED_INDEPENDENTLY" if (len(residual) == 0 and intersection == asserted)
                   else "CONFLICT_DO_NOT_POOL"),
    }


def load_and_verify():
    """Read the real upstream result files and build verified Source records.

    Every pinned number below is checked against the file it came from.  A
    disagreement raises rather than silently proceeding.
    """
    disagreements = []

    # ---- GEUVADIS ---------------------------------------------------------
    geu_file = os.path.join(
        ROOT, "project/working/direct_geuvadis_function_screen_v1/results/slc44a5_targeted_results.tsv"
    )
    geu = pd.read_csv(geu_file, sep="\t")
    g28 = geu[(geu.model_id == "MAGE_nonoverlap_replication") & (geu.ancestry_model == "population")]
    g33 = geu[(geu.model_id == "prespecified_MAGE_replication") & (geu.ancestry_model == "population")]
    if len(g28) != 1 or len(g33) != 1:
        raise SourceVerificationError("expected exactly one population-model row per GEUVADIS model_id")
    g28 = g28.iloc[0]
    g33 = g33.iloc[0]

    for label, row, exp_beta, exp_p, exp_n in [
        ("GEUVADIS-28", g28, 0.519734082414696, 0.00122798925532045, 28),
        ("GEUVADIS-33", g33, 0.609216023155312, 1.54608134461494e-05, 33),
    ]:
        if not _close(row.beta, exp_beta):
            disagreements.append(f"{label} beta: file={row.beta} briefed={exp_beta}")
        if not _close(row.p_two_sided_hc3, exp_p, 1e-6):
            disagreements.append(f"{label} P: file={row.p_two_sided_hc3} briefed={exp_p}")
        if int(row.n) != exp_n:
            disagreements.append(f"{label} n: file={row.n} briefed={exp_n}")

    # Recompute the GEUVADIS standardization from the raw materialized data so
    # that sd_x / sd_y are ours, not taken on trust.
    expr = pd.read_csv(
        os.path.join(ROOT, "project/working/direct_geuvadis_function_screen_v1/data/direct33_expression.tsv"),
        sep="\t",
    )
    expo = pd.read_csv(
        os.path.join(ROOT, "project/working/direct_geuvadis_function_screen_v1/data/direct33_exposures.tsv"),
        sep="\t",
    )
    # Donor disjointness is derived independently, not taken from the package's
    # own overlap assertion.  Fail closed on any conflict: the primary pool is
    # only meaningful if the two rosters truly share no donor.
    disjoint = independent_disjointness_derivation()
    if disjoint["status"] != "VERIFIED_INDEPENDENTLY":
        raise SourceVerificationError(
            "donor-disjointness could not be verified independently; refusing to pool. "
            f"derived_shared={disjoint['independently_derived_shared_donors']} "
            f"asserted_shared={disjoint['asserted_shared_donors']} "
            f"residual_intersection={disjoint['intersection_after_removal']}"
        )
    shared = disjoint["independently_derived_shared_donors"]
    ids33 = list(expo["sample"])
    ids28 = [s for s in ids33 if s not in set(shared)]
    if len(shared) != 5 or len(ids28) != 28:
        raise SourceVerificationError(f"donor-disjointness construction wrong: shared={shared} n28={len(ids28)}")

    yrow = expr[expr.gene_symbol == "SLC44A5"]
    if len(yrow) != 1:
        raise SourceVerificationError(f"expected exactly one SLC44A5 row, found {len(yrow)}")
    y = yrow[ids33].iloc[0].astype(float)
    x = expo.set_index("sample")["HML-2_1p31.1b::internal_fragment_present"].astype(float)

    sd_y28 = float(y[ids28].std(ddof=1))
    sd_x28 = float(x[ids28].std(ddof=1))
    sd_y33 = float(y[ids33].std(ddof=1))
    sd_x33 = float(x[ids33].std(ddof=1))

    # our recomputed standardized beta must match the package's published one
    for label, row, sdx, sdy in [
        ("GEUVADIS-28", g28, sd_x28, sd_y28),
        ("GEUVADIS-33", g33, sd_x33, sd_y33),
    ]:
        ours, _ = standardize(float(row.beta), float(row.se_hc3), sdx, sdy)
        if not _close(ours, float(row.standardized_beta), 1e-8):
            disagreements.append(
                f"{label} standardized_beta: recomputed={ours} file={row.standardized_beta}"
            )

    # ---- MAGE -------------------------------------------------------------
    mage_sens_file = os.path.join(
        ROOT, "project/working/onep31b_slc44a5_followup_v1/robust_stats/results/model_sensitivities.tsv"
    )
    ms = pd.read_csv(mage_sens_file, sep="\t")
    prim = ms[ms.model_id == "primary_reproduction"]
    if len(prim) != 1:
        raise SourceVerificationError("expected exactly one MAGE primary_reproduction row")
    prim = prim.iloc[0]

    mage_align = pd.read_csv(
        os.path.join(ROOT, "project/working/onep31b_slc44a5_followup_v1/robust_stats/results/sample_alignment.tsv"),
        sep="\t",
    )
    if len(mage_align) != 39:
        disagreements.append(f"MAGE n: alignment rows={len(mage_align)} expected=39")
    my = mage_align["SLC44A5_inverse_normal_expression"].astype(float)
    mx = mage_align["exposure"].astype(float)
    sd_y39 = float(my.std(ddof=1))
    sd_x39 = float(mx.std(ddof=1))

    # cross-check against the independent targeted screen, which publishes the
    # same coefficient already multiplied by sd(x)
    tgt = pd.read_csv(
        os.path.join(ROOT, "project/working/type1_functional_evidence_synthesis_v1/results/mage_targeted_slc44a5.tsv"),
        sep="\t",
    )
    t14 = tgt[tgt.representative_alias == "1p31.1b::internal_fragment_presence::copy_resolved"]
    if len(t14) != 1:
        raise SourceVerificationError("expected exactly one MAGE targeted 1p31.1b internal-fragment row")
    t14 = t14.iloc[0]
    implied = float(prim.beta) * sd_x39
    if not _close(implied, float(t14.beta_per_exposure_sd), 1e-8):
        disagreements.append(
            f"MAGE beta_per_exposure_sd: beta*sd(x)={implied} targeted-screen={t14.beta_per_exposure_sd}"
        )
    if not _close(float(prim.conventional_p), float(t14.p_value), 1e-8):
        disagreements.append(
            f"MAGE P disagreement between packages: robust_stats={prim.conventional_p} targeted={t14.p_value}"
        )

    # ---- exposure-contrast identity check --------------------------------
    states = pd.read_csv(os.path.join(ROOT, "project/working/locus_marker_expansion_v1/person_locus_states.tsv"), sep="\t")
    s1 = states[states.locus == "1p31.1b"]
    g_states = s1[s1["sample"].isin(ids28)]
    if len(g_states) != 28:
        raise SourceVerificationError(f"1p31.1b state rows for GEUVADIS-28: {len(g_states)}")
    same_contrast = (
        bool((g_states.retained_element_haplotype_dose == 2).all())
        and bool(((g_states.solo_ltr_haplotype_dose + g_states.internal_element_haplotype_dose) == 2).all())
        and bool((mage_align.retained_haplotype_dose == 2).all())
        and bool(((mage_align.solo_ltr_haplotype_dose + mage_align.internal_haplotype_dose) == 2).all())
    )
    if not same_contrast:
        raise SourceVerificationError(
            "the two cohorts do not share the identical two-retained-allele reference contrast"
        )

    geu_dose = g_states.internal_element_haplotype_dose.value_counts().reindex([0, 1, 2]).fillna(0).astype(int)
    mage_dose = mage_align.internal_haplotype_dose.value_counts().reindex([0, 1, 2]).fillna(0).astype(int)

    # ---- build records ----------------------------------------------------
    mage = Source(
        cohort="MAGE_v1_direct39",
        n=39,
        endpoint_class=ENDPOINT_STEADY_STATE_SLC44A5_LCL,
        endpoint_detail="SLC44A5 (ENSG00000137968) steady-state mRNA, 1000 Genomes EBV-LCL RNA-seq, MAGE v1",
        exposure_coding=(
            "binary: >=1 haplotype carrying a 1p31.1b Type-I internal fragment (13) vs two solo-LTR "
            "alleles (26); all 39 donors carry two retained 1p31.1b alleles; no array carriers"
        ),
        expression_normalization="inverse-normal transformed TMM (inverse_normal_TMM.filtered.TSS.MAGE.v1.0), no PEER pre-residualization",
        covariates="sex + genotype PC1-5 + PEER1-5 (in model)",
        se_kind="HC3",
        beta_raw=float(prim.beta),
        se_raw=float(prim.hc3_se),
        p_raw=float(prim.hc3_p),
        resid_df=int(prim.residual_df),
        n_covariates=11,
        sd_x=sd_x39,
        sd_y=sd_y39,
        correction_family=(
            "prespecified cis lead; also charged against a broadened all-gene targeted family in "
            "type1_functional_evidence_synthesis_v1 (targeted-family BH=0.127, no global hit)"
        ),
        source_file=mage_sens_file,
        harmonization="beta*sd(x)/sd(y); HC3 SE adopted to match the GEUVADIS SE convention",
    )
    mage.beta_std, mage.se_std = standardize(mage.beta_raw, mage.se_raw, mage.sd_x, mage.sd_y)

    geuv28 = Source(
        cohort="GEUVADIS_28_MAGE_disjoint",
        n=28,
        endpoint_class=ENDPOINT_STEADY_STATE_SLC44A5_LCL,
        endpoint_detail="SLC44A5 (ENSG00000137968) steady-state mRNA, 1000 Genomes EBV-LCL RNA-seq, E-GEUV-1",
        exposure_coding=(
            "binary: >=1 haplotype carrying a 1p31.1b Type-I internal fragment (14) vs two solo-LTR "
            "alleles (14); all 28 donors carry two retained 1p31.1b alleles; 1 array carrier (HG00329, CN=5)"
        ),
        expression_normalization="RPKM, quantile-normalized, PEER k=10 residuals (GEUVADIS_gene_RPKM_50FN_resk10)",
        covariates="population + sex (PEER factors already removed from the outcome)",
        se_kind="HC3",
        beta_raw=float(g28.beta),
        se_raw=float(g28.se_hc3),
        p_raw=float(g28.p_two_sided_hc3),
        resid_df=int(g28.residual_df),
        n_covariates=4,
        sd_x=sd_x28,
        sd_y=sd_y28,
        correction_family="independent_donor_sensitivity; q_bh_secondary_family=0.0184198388298067",
        source_file=geu_file,
        harmonization="beta*sd(x)/sd(y); reproduced the package's own standardized_beta exactly",
    )
    geuv28.beta_std, geuv28.se_std = standardize(geuv28.beta_raw, geuv28.se_raw, geuv28.sd_x, geuv28.sd_y)

    geuv33 = Source(
        cohort="GEUVADIS_33_NOT_disjoint",
        n=33,
        endpoint_class=ENDPOINT_STEADY_STATE_SLC44A5_LCL,
        endpoint_detail="SLC44A5 steady-state mRNA, E-GEUV-1; SHARES 5 DONORS WITH MAGE - secondary only",
        exposure_coding="binary internal-fragment presence (17) vs two solo-LTR alleles (16)",
        expression_normalization="RPKM, quantile-normalized, PEER k=10 residuals",
        covariates="population + sex",
        se_kind="HC3",
        beta_raw=float(g33.beta),
        se_raw=float(g33.se_hc3),
        p_raw=float(g33.p_two_sided_hc3),
        resid_df=int(g33.residual_df),
        n_covariates=4,
        sd_x=sd_x33,
        sd_y=sd_y33,
        correction_family="prespecified_replication",
        source_file=geu_file,
        harmonization="NOT DONOR-DISJOINT FROM MAGE (5 shared donors) - never used in the primary pool",
    )
    geuv33.beta_std, geuv33.se_std = standardize(geuv33.beta_raw, geuv33.se_raw, geuv33.sd_x, geuv33.sd_y)

    # MAGE with its own published conventional-OLS SE, for the SE-convention sensitivity
    mage_conv = Source(**{**asdict(mage), "se_kind": "conventional_OLS",
                          "se_raw": float(prim.conventional_se), "p_raw": float(prim.conventional_p)})
    mage_conv.beta_std, mage_conv.se_std = standardize(
        mage_conv.beta_raw, mage_conv.se_raw, mage_conv.sd_x, mage_conv.sd_y
    )

    facts = {
        "disjointness": disjoint,
        "shared_donors": shared,
        "geuvadis_dose_counts": {str(k): int(v) for k, v in geu_dose.items()},
        "mage_dose_counts": {str(k): int(v) for k, v in mage_dose.items()},
        "identical_reference_contrast": same_contrast,
        "disagreements": disagreements,
    }
    return {"mage": mage, "mage_conv": mage_conv, "geuv28": geuv28, "geuv33": geuv33}, facts


# --------------------------------------------------------------------------
# Leave-source-out
# --------------------------------------------------------------------------

def leave_source_out(sources):
    """With k sources, drop each in turn.

    With exactly two sources every 'jackknife' replicate is simply the other
    single source, unpooled.  That is reported literally, not as evidence of
    stability.
    """
    rows = []
    for i, s in enumerate(sources):
        kept = [t for j, t in enumerate(sources) if j != i]
        if len(kept) >= 2:
            res = synthesize(kept)["fixed_effect"]
            est, se, z, p = res["estimate"], res["se"], res["z"], res["p"]
            kind = "pooled_remaining_sources"
        else:
            k = kept[0]
            est, se = k.beta_std, k.se_std
            z = est / se
            p = float(2.0 * stats.norm.sf(abs(z)))
            kind = "single_remaining_source_NOT_a_pool"
        rows.append({
            "dropped_source": s.cohort,
            "remaining_sources": ";".join(t.cohort for t in kept),
            "n_remaining_sources": len(kept),
            "result_kind": kind,
            "estimate_std": est,
            "se_std": se,
            "ci95_low": est - 1.959963984540054 * se,
            "ci95_high": est + 1.959963984540054 * se,
            "z": z,
            "p": p,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    os.makedirs(OUT, exist_ok=True)
    src, facts = load_and_verify()
    mage, mage_conv, geuv28, geuv33 = src["mage"], src["mage_conv"], src["geuv28"], src["geuv33"]

    if facts["disagreements"]:
        print("!! UPSTREAM DISAGREEMENTS:")
        for d in facts["disagreements"]:
            print("   -", d)
    else:
        print("upstream verification: all pinned values match the real result files")

    primary_pair = [mage, geuv28]

    # ---- 1. inputs table --------------------------------------------------
    dj = facts["disjointness"]["status"]
    rows = []
    for s, role, djs in [
        (mage, "PRIMARY_disjoint_pair", dj),
        (geuv28, "PRIMARY_disjoint_pair", dj),
        (geuv33, "SECONDARY_not_disjoint", "NOT_DISJOINT_shares_5_donors_with_MAGE"),
        (mage_conv, "SENSITIVITY_se_convention", dj),
    ]:
        d = asdict(s)
        d["role"] = role
        d["donor_disjointness_status"] = djs
        d["standardization_factor_k_sdx_over_sdy"] = s.sd_x / s.sd_y
        d["z"] = s.beta_std / s.se_std
        d["partial_r"] = partial_correlation(s.beta_raw / s.se_raw, s.resid_df)
        rows.append(d)
    inputs = pd.DataFrame(rows)[[
        "role", "cohort", "n", "donor_disjointness_status", "endpoint_class", "endpoint_detail",
        "exposure_coding", "expression_normalization", "covariates", "se_kind", "beta_raw", "se_raw",
        "p_raw", "resid_df", "n_covariates", "sd_x", "sd_y", "standardization_factor_k_sdx_over_sdy",
        "beta_std", "se_std", "z", "partial_r", "correction_family", "source_file", "harmonization",
    ]]
    inputs.to_csv(os.path.join(OUT, "synthesis_slc44a5_inputs.tsv"), sep="\t", index=False)
    print("\n--- harmonized inputs (standardized scale: SD expression per SD exposure) ---")
    for _, r in inputs.iterrows():
        print(f"  {r['role']:28s} {r['cohort']:26s} n={r['n']:3d}  beta_std={r['beta_std']:+.6f}  "
              f"se_std={r['se_std']:.6f}  P={r['p_raw']:.4g}")

    # ---- 2. pooled --------------------------------------------------------
    primary = synthesize(primary_pair)
    secondary_nondisjoint = synthesize([mage, geuv33])
    sens_se = synthesize([mage_conv, geuv28])

    # scale-free robustness: partial correlation on Fisher's z
    zs, zses = [], []
    for s in primary_pair:
        r = partial_correlation(s.beta_raw / s.se_raw, s.resid_df)
        zs.append(fisher_z(r))
        zses.append(fisher_z_se(s.n, s.n_covariates))
    fe_z = fixed_effect(zs, zses)
    het_z = cochran_q(zs, zses)
    fisher_block = {
        "convention": "partial correlation r = t/sqrt(t^2+df), Fisher z = atanh(r), SE = 1/sqrt(n-k-3)",
        "why": "invariant to expression normalization and exposure scaling; removes the PEER asymmetry "
               "between the two cohorts",
        "per_source": [
            {"cohort": s.cohort, "partial_r": partial_correlation(s.beta_raw / s.se_raw, s.resid_df),
             "fisher_z": z, "fisher_z_se": zse}
            for s, z, zse in zip(primary_pair, zs, zses)
        ],
        "fixed_effect_z_scale": fe_z,
        "pooled_r": float(math.tanh(fe_z["estimate"])),
        "pooled_r_ci95": [float(math.tanh(fe_z["ci95_low"])), float(math.tanh(fe_z["ci95_high"]))],
        "heterogeneity_z_scale": het_z,
    }

    directions = {s.cohort: ("positive" if s.beta_std > 0 else "negative") for s in primary_pair}
    concordant = len(set(directions.values())) == 1

    pooled = {
        "analysis": "standardized-effect synthesis, 1p31.1b Type-I internal fragment -> SLC44A5 expression",
        "generated_by": "project/working/hml2_functional_evidence_synthesis_claude_v2/run_synthesis.py",
        "standardization_convention": {
            "scale": "SD of expression per SD of exposure",
            "algebra": "beta_std = beta_raw * sd(x)/sd(y); se_std = se_raw * sd(x)/sd(y); z invariant",
            "sd_x": "marginal SD of the binary exposure within the analysis sample",
            "sd_y": "marginal SD of the expression outcome within the analysis sample",
            "exposure_coding": ">=1 internal-fragment haplotype vs two solo-LTR alleles (identical in both cohorts)",
            "se_convention": "HC3 in both sources for the primary pool",
        },
        "endpoint_class_check": {
            "class": primary["endpoint_class"],
            "both_sources_same_class": True,
            "verdict": "POOLING PERMITTED - both are steady-state SLC44A5 mRNA in EBV-transformed LCLs",
        },
        "donor_disjointness_assertion_used": {
            "assertion": "the 28 GEUVADIS donors and the 39 MAGE donors share no individual",
            "status": facts["disjointness"]["status"],
            "verified_how": "INDEPENDENTLY DERIVED by this lane: the MAGE analysis roster "
                            "(onep31b_slc44a5_followup_v1/robust_stats/results/sample_alignment.tsv, 39 fitted "
                            "identities) was intersected directly with the GEUVADIS analysis roster "
                            "(direct_geuvadis_function_screen_v1/data/direct33_exposures.tsv, 33 fitted donors). "
                            "mage_overlap_audit.tsv was NOT used as an input to this derivation - it is the "
                            "package's own assertion and using it would be circular. It was compared afterwards "
                            "and agrees exactly.",
            "independent_derivation": facts["disjointness"],
            "shared_donors_removed": facts["shared_donors"],
            "intersection_after_removal": facts["disjointness"]["intersection_after_removal"],
            "n_geuvadis_after_removal": 28,
            "n_mage": 39,
            "contingency": "If any further independent derivation (e.g. lane T1-L) returns anything other "
                           "than exactly these 5 shared donors / 28 remaining, the pooled estimate must be "
                           "recomputed. run_synthesis.py fails closed rather than pooling on a conflict.",
            "cross_cohort_first_degree_relatives": 0,
            "relatedness_source": "HML2_ProjectResources/data/ref/1kGP.3202_samples.pedigree_info.txt "
                                  "(all 28 and all 39 donors present in the pedigree)",
            "residual_caveat": "first-degree pedigree only; cryptic/distant relatedness and shared "
                               "population structure are not excluded",
        },
        "exposure_contrast_identity": {
            "identical_reference_contrast": facts["identical_reference_contrast"],
            "detail": "every donor in both cohorts carries two retained 1p31.1b alleles; reference = two "
                      "solo-LTR alleles; this is NOT insertion presence vs absence",
            "geuvadis28_internal_dose_counts": facts["geuvadis_dose_counts"],
            "mage39_internal_dose_counts": facts["mage_dose_counts"],
            "composition_difference": "GEUVADIS-28 carriers are 10 het + 4 hom incl. 1 array carrier "
                                      "(HG00329, CN=5); MAGE-39 carriers are 12 het + 1 hom with no array "
                                      "carriers. The binary contrast is identical; the dose tail is not.",
        },
        "direction_concordance": {"per_source": directions, "concordant": concordant},
        "PRIMARY_donor_disjoint": primary,
        "SECONDARY_not_donor_disjoint_33": {
            "WARNING": "GEUVADIS-33 shares 5 donors with MAGE; this pool double-counts them and is NOT "
                       "independent replication. Reported for completeness only.",
            **secondary_nondisjoint,
        },
        "SENSITIVITY_mage_conventional_se": {
            "note": "MAGE's own published SE convention (conventional OLS) instead of HC3; GEUVADIS stays HC3",
            **sens_se,
        },
        "ROBUSTNESS_partial_correlation_scale_free": fisher_block,
        "heterogeneity_power_caveat": (
            "With k=2 sources Cochran's Q has df=1 and heterogeneity estimation is very low-power. "
            "I^2 and tau^2 here are near-uninformative: a null Q does NOT establish homogeneity, and "
            "I^2=0 must not be read as evidence that the two effects are equal. Do not over-interpret."
        ),
        "upstream_disagreements": facts["disagreements"],
    }
    with open(os.path.join(OUT, "synthesis_slc44a5_pooled.json"), "w") as fh:
        json.dump(pooled, fh, indent=2)

    fe = primary["fixed_effect"]
    dl = primary["random_effects_dersimonian_laird"]
    pm = primary["random_effects_paule_mandel"]
    het = primary["heterogeneity"]
    print("\n--- PRIMARY donor-disjoint pool (MAGE-39 + GEUVADIS-28) ---")
    print(f"  fixed effect  : {fe['estimate']:+.6f}  SE={fe['se']:.6f}  "
          f"95% CI [{fe['ci95_low']:+.6f}, {fe['ci95_high']:+.6f}]  z={fe['z']:.4f}  P={fe['p']:.4g}")
    print(f"  RE (DL)       : {dl['estimate']:+.6f}  SE={dl['se']:.6f}  "
          f"95% CI [{dl['ci95_low']:+.6f}, {dl['ci95_high']:+.6f}]  z={dl['z']:.4f}  P={dl['p']:.4g}  tau2={dl['tau2']:.6g}")
    print(f"  RE (PM)       : {pm['estimate']:+.6f}  SE={pm['se']:.6f}  z={pm['z']:.4f}  P={pm['p']:.4g}  tau2={pm['tau2']:.6g}")
    print(f"  heterogeneity : Q={het['Q']:.6f}  df={het['df']}  P_Q={het['p']:.4f}  I2={het['I2_percent']:.2f}%")
    print(f"  partial-r pool: r={fisher_block['pooled_r']:.4f} "
          f"CI [{fisher_block['pooled_r_ci95'][0]:.4f}, {fisher_block['pooled_r_ci95'][1]:.4f}] "
          f"P={fe_z['p']:.4g}")

    # ---- 3. leave-source-out ---------------------------------------------
    lso = leave_source_out(primary_pair)
    lso.to_csv(os.path.join(OUT, "synthesis_slc44a5_leave_source_out.tsv"), sep="\t", index=False)
    print("\n--- leave-source-out (k=2: each row is one source alone, not a jackknife) ---")
    for _, r in lso.iterrows():
        print(f"  drop {r['dropped_source']:26s} -> {r['remaining_sources']:26s} "
              f"est={r['estimate_std']:+.6f} P={r['p']:.4g}  [{r['result_kind']}]")

    print("\nwrote:", OUT)
    return pooled


if __name__ == "__main__":
    main()
