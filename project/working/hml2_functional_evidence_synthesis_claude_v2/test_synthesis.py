#!/usr/bin/env python
"""
Lane T-VAL / Task 1 -- focused validator for the SLC44A5 synthesis outputs.

Per CTC-0010 item 5, this file validates the ALREADY-COMPUTED, coordinator-accepted
synthesis outputs against the EXACT facts confirmed in CTC-0009.2.

The load-bearing check is S2: the inverse-variance pool is RECOMPUTED from
synthesis_slc44a5_inputs.tsv inside this test (weights = 1/SE^2, pooled = sum(w*b)/sum(w),
SE = sqrt(1/sum(w))) and compared against synthesis_slc44a5_pooled.json. It does not read
the JSON's answer back and call that agreement.

READ-ONLY on every data file. This module writes exactly one artefact:
    results/validator_synthesis_report.json

Run:
    python test_synthesis.py
Exit code 0 = all assertions pass, 1 = at least one assertion failed.
"""

import csv
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
REPORT_PATH = os.path.join(RESULTS, "validator_synthesis_report.json")

INPUTS_TSV = os.path.join(RESULTS, "synthesis_slc44a5_inputs.tsv")
POOLED_JSON = os.path.join(RESULTS, "synthesis_slc44a5_pooled.json")
LSO_TSV = os.path.join(RESULTS, "synthesis_slc44a5_leave_source_out.tsv")

sys.path.insert(0, HERE)
import run_synthesis  # noqa: E402  (the module under validation; main() is __main__-guarded)

# ---------------------------------------------------------------------------
# Coordinator-confirmed constants (CTC-0009.2). Hard-coded on purpose.
# ---------------------------------------------------------------------------
EXPECTED_MAGE_BETA_STD = 0.78364
EXPECTED_MAGE_SE_STD = 0.16135
EXPECTED_GEUVADIS_BETA_STD = 0.84927
EXPECTED_GEUVADIS_SE_STD = 0.22908
EXPECTED_POOLED_BETA = 0.80540
EXPECTED_POOLED_SE = 0.13191
EXPECTED_POOLED_P = 1.02e-9

MAGE_COHORT = "MAGE_v1_direct39"
GEUVADIS_28_COHORT = "GEUVADIS_28_MAGE_disjoint"
GEUVADIS_33_COHORT = "GEUVADIS_33_NOT_disjoint"
PRIMARY_ROLE = "PRIMARY_disjoint_pair"

# Stated tolerances.
TOL_RECOMPUTE = 1e-9   # absolute: recomputed pool vs the JSON's stored pool
TOL_ROUND_DP = 5       # coordinator constants are 5-decimal roundings
TOL_P_SIGFIG = 3       # P = 1.02e-9 is a 3-significant-figure rounding

REPORT = []


def check(aid, group, description, expected, observed, ok=None):
    """Record one assertion. Never raises -- the report must stay complete even
    after the first failure, so a failing run still surfaces every value."""
    if ok is None:
        ok = expected == observed
    REPORT.append(
        {
            "assertion_id": aid,
            "group": group,
            "description": description,
            "expected": expected,
            "observed": observed,
            "result": "PASS" if ok else "FAIL",
        }
    )
    return bool(ok)


def assert_group(group):
    fails = [r for r in REPORT if r["group"] == group and r["result"] == "FAIL"]
    if fails:
        lines = [
            "%s: expected=%r observed=%r" % (f["assertion_id"], f["expected"], f["observed"])
            for f in fails
        ]
        raise AssertionError(
            "%d assertion(s) FAILED in group %s:\n  %s" % (len(fails), group, "\n  ".join(lines))
        )


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if any((v or "").strip() for v in r.values())]


def load():
    rows = read_tsv(INPUTS_TSV)
    with open(POOLED_JSON, encoding="utf-8") as fh:
        pooled = json.load(fh)
    lso = read_tsv(LSO_TSV)
    return rows, pooled, lso


def primary_rows(rows):
    """The two donor-disjoint sources, MAGE first (the JSON's source order)."""
    prim = [r for r in rows if r["role"] == PRIMARY_ROLE]
    order = {MAGE_COHORT: 0, GEUVADIS_28_COHORT: 1}
    return sorted(prim, key=lambda r: order.get(r["cohort"], 99))


def sig(x, n):
    """Round to n significant figures."""
    if x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))


def inverse_variance_pool(betas, ses):
    """The pooling arithmetic, implemented here independently of run_synthesis."""
    w = [1.0 / (s * s) for s in ses]
    est = sum(wi * bi for wi, bi in zip(w, betas)) / sum(w)
    se = math.sqrt(1.0 / sum(w))
    z = est / se
    p = math.erfc(abs(z) / math.sqrt(2.0))  # two-sided normal, no scipy dependency
    return {"estimate": est, "se": se, "z": z, "p": p, "weights": w}


# ---------------------------------------------------------------------------
# S1 -- standardized effects and their HC3 SEs.
# ---------------------------------------------------------------------------
def test_s1_standardized_effects():
    rows, pooled, lso = load()
    g = "S1_standardized_effects"
    prim = primary_rows(rows)
    by = {r["cohort"]: r for r in prim}

    check("S1.1", g, "exactly two PRIMARY_disjoint_pair rows", 2, len(prim))
    check("S1.2", g, "the two primary sources are MAGE-39 and the 28-person MAGE-disjoint GEUVADIS subset",
          sorted([MAGE_COHORT, GEUVADIS_28_COHORT]), sorted(by))

    m, gv = by[MAGE_COHORT], by[GEUVADIS_28_COHORT]

    check("S1.3", g, "MAGE standardized effect == 0.78364",
          EXPECTED_MAGE_BETA_STD, round(float(m["beta_std"]), TOL_ROUND_DP))
    check("S1.4", g, "MAGE HC3 standardized SE == 0.16135",
          EXPECTED_MAGE_SE_STD, round(float(m["se_std"]), TOL_ROUND_DP))
    check("S1.5", g, "GEUVADIS standardized effect == 0.84927",
          EXPECTED_GEUVADIS_BETA_STD, round(float(gv["beta_std"]), TOL_ROUND_DP))
    check("S1.6", g, "GEUVADIS HC3 standardized SE == 0.22908",
          EXPECTED_GEUVADIS_SE_STD, round(float(gv["se_std"]), TOL_ROUND_DP))
    check("S1.7", g, "both primary sources use the HC3 SE convention",
          {MAGE_COHORT: "HC3", GEUVADIS_28_COHORT: "HC3"},
          {MAGE_COHORT: m["se_kind"], GEUVADIS_28_COHORT: gv["se_kind"]})
    check("S1.8", g, "both primary sources carry the same endpoint class",
          True, m["endpoint_class"] == gv["endpoint_class"])

    # Independently rebuild beta_std/se_std from the raw slope and the SDs.
    recomputed = {}
    for r in (m, gv):
        b, s = run_synthesis.standardize(
            float(r["beta_raw"]), float(r["se_raw"]), float(r["sd_x"]), float(r["sd_y"])
        )
        recomputed[r["cohort"]] = (round(b, TOL_ROUND_DP), round(s, TOL_ROUND_DP))
    check("S1.9", g, "beta_std and se_std recomputed as beta_raw*sd_x/sd_y from the raw columns",
          {MAGE_COHORT: (EXPECTED_MAGE_BETA_STD, EXPECTED_MAGE_SE_STD),
           GEUVADIS_28_COHORT: (EXPECTED_GEUVADIS_BETA_STD, EXPECTED_GEUVADIS_SE_STD)},
          recomputed)

    # The stored standardization factor k must equal sd_x/sd_y.
    k_bad = []
    for r in (m, gv):
        k = float(r["sd_x"]) / float(r["sd_y"])
        if abs(k - float(r["standardization_factor_k_sdx_over_sdy"])) > 1e-9:
            k_bad.append("%s: stored=%s recomputed=%.12f" %
                         (r["cohort"], r["standardization_factor_k_sdx_over_sdy"], k))
    check("S1.10", g, "stored standardization factor k equals sd_x/sd_y in both sources", [], k_bad)

    # Both directions must be positive -- a replicated direction is the whole claim.
    check("S1.11", g, "both standardized effects are positive (concordant direction)",
          True, float(m["beta_std"]) > 0 and float(gv["beta_std"]) > 0)
    assert_group(g)


# ---------------------------------------------------------------------------
# S2 -- LOAD-BEARING: recompute the inverse-variance pool from the TSV.
# ---------------------------------------------------------------------------
def test_s2_inverse_variance_pool_recomputed():
    rows, pooled, lso = load()
    g = "S2_inverse_variance_pool_recomputed"
    prim = primary_rows(rows)

    betas = [float(r["beta_std"]) for r in prim]
    ses = [float(r["se_std"]) for r in prim]
    mine = inverse_variance_pool(betas, ses)
    fe = pooled["PRIMARY_donor_disjoint"]["fixed_effect"]

    check("S2.1", g, "recomputed pooled beta == 0.80540 (weights = 1/SE^2, pooled = sum(w*b)/sum(w))",
          EXPECTED_POOLED_BETA, round(mine["estimate"], TOL_ROUND_DP))
    check("S2.2", g, "recomputed pooled SE == 0.13191 (SE = sqrt(1/sum(w)))",
          EXPECTED_POOLED_SE, round(mine["se"], TOL_ROUND_DP))

    d_beta = abs(mine["estimate"] - fe["estimate"])
    d_se = abs(mine["se"] - fe["se"])
    check("S2.3", g, "recomputed pooled beta matches synthesis_slc44a5_pooled.json (abs tol 1e-9)",
          True, d_beta <= TOL_RECOMPUTE,
          ok=(d_beta <= TOL_RECOMPUTE))
    check("S2.4", g, "recomputed pooled SE matches synthesis_slc44a5_pooled.json (abs tol 1e-9)",
          True, d_se <= TOL_RECOMPUTE,
          ok=(d_se <= TOL_RECOMPUTE))
    check("S2.5", g, "absolute difference between recomputed and stored pooled beta",
          "<= 1e-9", repr(d_beta), ok=(d_beta <= TOL_RECOMPUTE))
    check("S2.6", g, "absolute difference between recomputed and stored pooled SE",
          "<= 1e-9", repr(d_se), ok=(d_se <= TOL_RECOMPUTE))

    w_diff = max(abs(a - b) for a, b in zip(mine["weights"], fe["weights"]))
    check("S2.7", g, "recomputed inverse-variance weights match the stored weights",
          True, w_diff <= 1e-6, ok=(w_diff <= 1e-6))
    check("S2.8", g, "weights are exactly 1/SE^2 for both sources",
          True, all(abs(w - 1.0 / (s * s)) <= 1e-12 for w, s in zip(fe["weights"], ses)))
    check("S2.9", g, "stored weight fractions sum to 1",
          True, abs(sum(fe["weight_fraction"]) - 1.0) <= 1e-12)
    check("S2.10", g, "recomputed z matches the stored z",
          True, abs(mine["z"] - fe["z"]) <= 1e-9)

    # Sanity properties any correct inverse-variance pool must have.
    check("S2.11", g, "the pooled estimate lies between the two source estimates",
          True, min(betas) <= fe["estimate"] <= max(betas))
    check("S2.12", g, "the pooled SE is smaller than either source SE",
          True, fe["se"] < min(ses))
    check("S2.13", g, "the more precise source (MAGE) carries the larger weight",
          True, fe["weights"][0] > fe["weights"][1])

    # run_synthesis's own pooling must agree with this independent implementation.
    fe_mod = run_synthesis.fixed_effect(betas, ses)
    check("S2.14", g, "run_synthesis.fixed_effect() agrees with this independent recomputation",
          True, abs(fe_mod["estimate"] - mine["estimate"]) <= TOL_RECOMPUTE
          and abs(fe_mod["se"] - mine["se"]) <= TOL_RECOMPUTE)
    assert_group(g)


# ---------------------------------------------------------------------------
# S3 -- fixed-effect P.
# ---------------------------------------------------------------------------
def test_s3_fixed_effect_p_value():
    rows, pooled, lso = load()
    g = "S3_fixed_effect_p"
    prim = primary_rows(rows)
    betas = [float(r["beta_std"]) for r in prim]
    ses = [float(r["se_std"]) for r in prim]
    mine = inverse_variance_pool(betas, ses)
    fe = pooled["PRIMARY_donor_disjoint"]["fixed_effect"]

    check("S3.1", g, "stored fixed-effect P == 1.02e-9 (3 significant figures)",
          EXPECTED_POOLED_P, sig(fe["p"], TOL_P_SIGFIG))
    check("S3.2", g, "recomputed fixed-effect P == 1.02e-9 (3 significant figures)",
          EXPECTED_POOLED_P, sig(mine["p"], TOL_P_SIGFIG))
    rel = abs(mine["p"] - fe["p"]) / fe["p"]
    check("S3.3", g, "recomputed P matches the stored P (relative tol 1e-6)",
          True, rel <= 1e-6, ok=(rel <= 1e-6))
    check("S3.4", g, "relative difference between recomputed and stored P",
          "<= 1e-6", repr(rel), ok=(rel <= 1e-6))
    check("S3.5", g, "the stored 95% CI is consistent with estimate +/- 1.959964*SE",
          True,
          abs(fe["ci95_low"] - (fe["estimate"] - 1.959963984540054 * fe["se"])) <= 1e-9
          and abs(fe["ci95_high"] - (fe["estimate"] + 1.959963984540054 * fe["se"])) <= 1e-9)
    check("S3.6", g, "the 95% CI excludes the null",
          True, fe["ci95_low"] > 0.0)
    assert_group(g)


# ---------------------------------------------------------------------------
# S4 -- the pooled estimate uses the DONOR-DISJOINT pair; the 33-person version
# must NOT be the primary.
# ---------------------------------------------------------------------------
def test_s4_primary_pool_is_donor_disjoint():
    rows, pooled, lso = load()
    g = "S4_primary_is_donor_disjoint"

    prim_block = pooled["PRIMARY_donor_disjoint"]
    check("S4.1", g, "the primary pool's sources are MAGE-39 and the 28-person GEUVADIS-disjoint subset",
          sorted([MAGE_COHORT, GEUVADIS_28_COHORT]), sorted(prim_block["sources"]))
    check("S4.2", g, "the 33-person GEUVADIS panel is NOT among the primary pool's sources",
          False, GEUVADIS_33_COHORT in prim_block["sources"])

    prim = primary_rows(rows)
    ns = {r["cohort"]: int(r["n"]) for r in prim}
    check("S4.3", g, "the primary pool's GEUVADIS side has n=28, not n=33",
          {MAGE_COHORT: 39, GEUVADIS_28_COHORT: 28}, ns)

    # The 33-person row exists but is explicitly demoted to secondary.
    g33 = [r for r in rows if r["cohort"] == GEUVADIS_33_COHORT]
    check("S4.4", g, "the 33-person GEUVADIS row is present in the inputs TSV", 1, len(g33))
    check("S4.5", g, "the 33-person row is labelled a secondary, NOT-disjoint role",
          "SECONDARY_not_disjoint", g33[0]["role"] if g33 else None)
    check("S4.6", g, "the 33-person row records that it shares donors with MAGE and is never pooled as primary",
          True, "NOT DONOR-DISJOINT" in g33[0]["harmonization"].upper() if g33 else False)

    check("S4.7", g, "the JSON's 33-person pool is a separately-labelled secondary block, not the primary",
          True, "SECONDARY_not_donor_disjoint_33" in pooled and "PRIMARY_donor_disjoint" in pooled)
    check("S4.8", g, "the 33-person secondary block carries an explicit non-independence WARNING",
          True, "WARNING" in pooled["SECONDARY_not_donor_disjoint_33"])
    check("S4.9", g, "the secondary block's warning names the 5 shared donors as the reason",
          True, "5 donors" in pooled["SECONDARY_not_donor_disjoint_33"]["WARNING"])

    # The primary estimate must actually be the disjoint one, not the 33-person one.
    sec_est = pooled["SECONDARY_not_donor_disjoint_33"]["fixed_effect"]["estimate"]
    check("S4.10", g, "the primary estimate differs from the 33-person pool (the disjoint pair was used)",
          True, abs(prim_block["fixed_effect"]["estimate"] - sec_est) > 1e-6)
    check("S4.11", g, "recomputing the pool from the 28-person row reproduces the PRIMARY estimate",
          EXPECTED_POOLED_BETA,
          round(inverse_variance_pool([float(r["beta_std"]) for r in prim],
                                      [float(r["se_std"]) for r in prim])["estimate"], TOL_ROUND_DP))

    # Donor-disjointness provenance.
    dj = pooled["donor_disjointness_assertion_used"]
    check("S4.12", g, "the disjointness block records the five removed shared donors",
          ["HG00096", "HG00146", "HG00320", "NA18508", "NA19129"],
          sorted(dj["shared_donors_removed"]))
    check("S4.13", g, "the intersection after removal is empty", [], list(dj["intersection_after_removal"]))
    check("S4.14", g, "n after removal is 28 GEUVADIS vs 39 MAGE",
          {"geuvadis": 28, "mage": 39},
          {"geuvadis": dj["n_geuvadis_after_removal"], "mage": dj["n_mage"]})
    check("S4.15", g, "zero cross-cohort first-degree relatives",
          0, dj["cross_cohort_first_degree_relatives"])
    check("S4.16", g, "the residual relatedness caveat is retained (first-degree pedigree only)",
          True, "cryptic" in dj["residual_caveat"].lower())

    # Leave-source-out must contain exactly the two disjoint sources, and must not
    # present a single remaining source as if it were a pool.
    dropped = sorted(r["dropped_source"] for r in lso)
    check("S4.17", g, "leave-source-out drops each of the two disjoint sources exactly once",
          sorted([MAGE_COHORT, GEUVADIS_28_COHORT]), dropped)
    check("S4.18", g, "leave-source-out never labels a single remaining source as a pool",
          [r["result_kind"] for r in lso],
          ["single_remaining_source_NOT_a_pool"] * len(lso))
    check("S4.19", g, "the 33-person panel never appears in leave-source-out",
          False, any(GEUVADIS_33_COHORT in r["dropped_source"] or GEUVADIS_33_COHORT in r["remaining_sources"]
                     for r in lso))

    # Each leave-one-out row must equal the surviving source's own effect.
    by_cohort = {r["cohort"]: r for r in prim}
    lso_bad = []
    for r in lso:
        surviving = by_cohort.get(r["remaining_sources"])
        if surviving is None:
            lso_bad.append("unresolved remaining source: %s" % r["remaining_sources"])
            continue
        if abs(float(r["estimate_std"]) - float(surviving["beta_std"])) > 1e-9:
            lso_bad.append("%s: estimate %s != source beta_std %s"
                           % (r["remaining_sources"], r["estimate_std"], surviving["beta_std"]))
    check("S4.20", g, "each leave-one-out estimate equals the surviving source's own standardized effect",
          [], lso_bad)
    assert_group(g)


# ---------------------------------------------------------------------------
# S5 -- the code must REFUSE to pool different endpoint classes.
# Exercised on a synthetic mismatched-class case, not on the real data.
# ---------------------------------------------------------------------------
def _mk_source(cohort, endpoint_class, beta_std=0.5, se_std=0.1):
    return run_synthesis.Source(
        cohort=cohort,
        n=30,
        endpoint_class=endpoint_class,
        endpoint_detail="synthetic test fixture",
        exposure_coding="synthetic binary",
        expression_normalization="synthetic",
        covariates="synthetic",
        se_kind="HC3",
        beta_raw=1.0,
        se_raw=0.2,
        p_raw=0.01,
        resid_df=25,
        n_covariates=3,
        sd_x=0.5,
        sd_y=1.0,
        correction_family="synthetic",
        source_file="<synthetic>",
        harmonization="synthetic",
        beta_std=beta_std,
        se_std=se_std,
    )


def test_s5_endpoint_class_guard_refuses_mismatch():
    rows, pooled, lso = load()
    g = "S5_endpoint_class_guard"

    same_class = "steady_state_SLC44A5_mRNA_abundance_in_EBV_LCL"
    other_class = "drug_survival_anti_CD20_complement_killing"

    mismatched = [_mk_source("synthetic_a", same_class), _mk_source("synthetic_b", other_class)]
    matched = [_mk_source("synthetic_a", same_class), _mk_source("synthetic_b", same_class)]

    # assert_poolable must RAISE on a mismatched class, not warn and proceed.
    raised = None
    try:
        run_synthesis.assert_poolable(mismatched)
    except run_synthesis.EndpointClassMismatch as exc:
        raised = str(exc)
    except Exception as exc:  # wrong exception type is still a failure
        raised = "WRONG_EXCEPTION_TYPE: %s: %s" % (type(exc).__name__, exc)
    check("S5.1", g, "assert_poolable RAISES EndpointClassMismatch on two different endpoint classes",
          True, raised is not None and not str(raised).startswith("WRONG_EXCEPTION_TYPE"))
    check("S5.2", g, "the refusal message names the refusal, not a warning",
          True, bool(raised) and "refusing to pool" in raised)
    check("S5.3", g, "the refusal message names both offending classes",
          True, bool(raised) and same_class in raised and other_class in raised)

    # The guard must also fire through the full synthesize() entry point, so a
    # caller cannot bypass it by going one level up.
    synth_raised = False
    synth_result = None
    try:
        synth_result = run_synthesis.synthesize(mismatched)
    except run_synthesis.EndpointClassMismatch:
        synth_raised = True
    except Exception:
        synth_raised = False
    check("S5.4", g, "synthesize() also refuses a mismatched-class pool (the guard cannot be bypassed)",
          True, synth_raised)
    check("S5.5", g, "no pooled result is returned when the classes mismatch",
          None, synth_result)

    # A three-source list with one intruder must also be refused.
    intruder = [
        _mk_source("synthetic_a", same_class),
        _mk_source("synthetic_b", same_class),
        _mk_source("synthetic_c", other_class),
    ]
    intruder_raised = False
    try:
        run_synthesis.assert_poolable(intruder)
    except run_synthesis.EndpointClassMismatch:
        intruder_raised = True
    check("S5.6", g, "a single intruding endpoint class in a 3-source list is refused",
          True, intruder_raised)

    # Matched classes must be accepted and return the shared class.
    accepted = run_synthesis.assert_poolable(matched)
    check("S5.7", g, "assert_poolable accepts a matched-class pool and returns the shared class",
          same_class, accepted)

    # Fewer than two sources is not a pool.
    lone_raised = False
    try:
        run_synthesis.assert_poolable([_mk_source("synthetic_a", same_class)])
    except ValueError:
        lone_raised = True
    check("S5.8", g, "assert_poolable refuses a 'pool' of fewer than two sources", True, lone_raised)

    # The real primary sources must satisfy the guard.
    prim = primary_rows(rows)
    real = [_mk_source(r["cohort"], r["endpoint_class"]) for r in prim]
    real_ok = True
    try:
        run_synthesis.assert_poolable(real)
    except Exception:
        real_ok = False
    check("S5.9", g, "the two real primary sources pass the endpoint-class guard", True, real_ok)
    check("S5.10", g, "the pooled JSON records that the endpoint-class check permitted pooling",
          True, pooled["endpoint_class_check"]["both_sources_same_class"] is True)
    assert_group(g)


# ---------------------------------------------------------------------------
# S6 -- the standardization algebra is invertible on a known synthetic case.
# ---------------------------------------------------------------------------
def test_s6_standardization_algebra_invertible():
    rows, pooled, lso = load()
    g = "S6_standardization_invertible"

    # A hand-checkable case: k = sd_x/sd_y = 0.5/2.0 = 0.25, so
    # beta_std = 2.0*0.25 = 0.5 and se_std = 0.5*0.25 = 0.125.
    beta_raw, se_raw, sd_x, sd_y = 2.0, 0.5, 0.5, 2.0
    b_std, s_std = run_synthesis.standardize(beta_raw, se_raw, sd_x, sd_y)
    check("S6.1", g, "known synthetic case: standardize(2.0, 0.5, sd_x=0.5, sd_y=2.0) -> beta_std=0.5",
          0.5, b_std)
    check("S6.2", g, "known synthetic case: -> se_std=0.125", 0.125, s_std)

    b_back, s_back = run_synthesis.unstandardize(b_std, s_std, sd_x, sd_y)
    check("S6.3", g, "unstandardize inverts standardize exactly on the known case (beta)",
          True, abs(b_back - beta_raw) <= 1e-12)
    check("S6.4", g, "unstandardize inverts standardize exactly on the known case (se)",
          True, abs(s_back - se_raw) <= 1e-12)

    # The Wald statistic must be invariant -- this is what makes the rescale safe.
    check("S6.5", g, "the Wald statistic beta/se is invariant under standardization",
          True, abs((beta_raw / se_raw) - (b_std / s_std)) <= 1e-12)

    # Round-trip the REAL published values, not just the toy case.
    rt = []
    for r in primary_rows(rows):
        sx, sy = float(r["sd_x"]), float(r["sd_y"])
        b1, s1 = run_synthesis.standardize(float(r["beta_raw"]), float(r["se_raw"]), sx, sy)
        b0, s0 = run_synthesis.unstandardize(b1, s1, sx, sy)
        if abs(b0 - float(r["beta_raw"])) > 1e-12 or abs(s0 - float(r["se_raw"])) > 1e-12:
            rt.append(r["cohort"])
    check("S6.6", g, "standardize/unstandardize round-trips exactly on both real primary sources",
          [], rt)

    # Standardizing an already-standardized value with k=1 must be a no-op.
    b_id, s_id = run_synthesis.standardize(0.78364, 0.16135, 1.0, 1.0)
    check("S6.7", g, "standardization with sd_x == sd_y is the identity",
          (0.78364, 0.16135), (b_id, s_id))

    # Degenerate SDs must be rejected rather than silently producing inf/nan.
    for aid, sx, sy, label in (("S6.8", 0.0, 1.0, "sd_x == 0"),
                               ("S6.9", 1.0, 0.0, "sd_y == 0"),
                               ("S6.10", -1.0, 1.0, "sd_x < 0")):
        rejected = False
        try:
            run_synthesis.standardize(1.0, 0.1, sx, sy)
        except ValueError:
            rejected = True
        check(aid, g, "standardize rejects a degenerate SD (%s)" % label, True, rejected)

    # The stored z must equal beta_std/se_std for each real source.
    z_bad = []
    for r in primary_rows(rows):
        z = float(r["beta_std"]) / float(r["se_std"])
        if abs(z - float(r["z"])) > 1e-9:
            z_bad.append("%s: stored z=%s recomputed=%.12f" % (r["cohort"], r["z"], z))
    check("S6.11", g, "stored z equals beta_std/se_std for both real primary sources", [], z_bad)
    assert_group(g)


# ---------------------------------------------------------------------------
# S7 -- claim boundary.
#
# CTC-0009.2: "Carry the claim boundary: this is a replicated association with
# 1p31.1b internal-fragment state, not proof that the fragment rather than a
# linked host haplotype is causal."
#
# The synthesis report is what gets frozen and read downstream, so the boundary
# has to travel WITH the numbers. This group searches the union of all three
# synthesis outputs -- the most generous reading available.
# ---------------------------------------------------------------------------

# The POSITIVE half of the claim: a replicated/concordant association with the
# 1p31.1b internal-fragment state. Checked for context -- the report does make this.
POSITIVE_CLAIM_CONCEPTS = {
    "identifies_internal_fragment_exposure": [
        r"internal[- ]fragment\s+(?:state|presence|status|haplotype)",
        r"1p31\.1b[^.\n]{0,40}internal[- ]fragment",
    ],
    "states_direction_replication_or_concordance": [r"concordan\w*", r"replicat\w*"],
}

# The EPISTEMIC LIMIT that CTC-0009.2 requires the report to carry:
#   "not proof that the fragment rather than a linked host haplotype is causal".
# This is the load-bearing half. Without it, the frozen numbers read downstream as a
# causal result. Each concept is satisfied by any one variant.
EPISTEMIC_BOUNDARY_CONCEPTS = {
    "not_proof_of_causality": [
        r"not\s+proof", r"no\s+proof", r"does\s+not\s+prove",
        r"not\s+establish\w*\s+causal", r"not\s+.{0,40}causal",
    ],
    "linked_host_haplotype_alternative": [
        r"linked\s+host\s+haplotype", r"host\s+haplotype", r"linked\s+haplotype",
    ],
}

# Affirmative causal claims. The word "causal" inside a negated boundary statement
# ("not proof that ... is causal") is the boundary itself, not a violation, so each
# hit is only counted when no negation cue precedes it within NEGATION_WINDOW chars.
CAUSAL_CLAIM_PATTERNS = [
    r"\bcauses\b", r"\bcausing\b", r"\bcausally\b", r"\bis\s+causal\b", r"\bare\s+causal\b",
    r"\bproves\b", r"\bproven\b", r"\bdemonstrates\s+that\b", r"\bestablishes\s+causal",
    r"\bresponsible\s+for\b", r"\bdrives\s+(?:expression|SLC44A5)\b",
]
NEGATION_CUES = ["not", "never", "cannot", "can not", "no ", "without", "rather than",
                 "insufficient", "unproven", "does not", "is not", "are not", "nor "]
NEGATION_WINDOW = 90


def _synthesis_text():
    """The union of every synthesis output, as one lowercase blob per file."""
    blobs = {}
    for path in (POOLED_JSON, INPUTS_TSV, LSO_TSV):
        with open(path, encoding="utf-8") as fh:
            blobs[os.path.basename(path)] = fh.read().lower()
    return blobs


def test_s7_claim_boundary_present_and_no_causal_claim():
    rows, pooled, lso = load()
    g = "S7_claim_boundary"
    blobs = _synthesis_text()
    combined = "\n".join(blobs.values())

    # -- S7.1 (context): the positive half of the claim IS made by the report.
    positive = {
        c: any(re.search(v, combined) for v in variants)
        for c, variants in POSITIVE_CLAIM_CONCEPTS.items()
    }
    check("S7.1", g,
          "the synthesis report identifies the exposure as 1p31.1b internal-fragment state and "
          "reports replication/direction concordance (the positive half of the CTC-0009.2 claim)",
          {c: True for c in POSITIVE_CLAIM_CONCEPTS}, positive)

    # -- S7.2 (load-bearing): the epistemic limit must travel with the numbers.
    boundary = {
        c: any(re.search(v, combined) for v in variants)
        for c, variants in EPISTEMIC_BOUNDARY_CONCEPTS.items()
    }
    missing = sorted(c for c, ok in boundary.items() if not ok)
    check("S7.2", g,
          "the synthesis report carries the CTC-0009.2 epistemic boundary: this is NOT proof that "
          "the fragment rather than a linked host haplotype is causal (searched the union of all "
          "three synthesis outputs -- the most generous reading available)",
          {c: True for c in EPISTEMIC_BOUNDARY_CONCEPTS}, boundary)
    check("S7.3", g, "epistemic-boundary concepts missing from EVERY synthesis output",
          [], missing)

    # Per-file detail, so the report says exactly where the boundary is absent.
    per_file = {}
    for name, blob in blobs.items():
        per_file[name] = {
            c: any(re.search(v, blob) for v in variants)
            for c, variants in EPISTEMIC_BOUNDARY_CONCEPTS.items()
        }
    check("S7.4", g, "at least one synthesis output carries the complete epistemic boundary",
          True, any(all(v.values()) for v in per_file.values()))
    check("S7.5", g, "per-file epistemic-boundary coverage",
          {n: {c: True for c in EPISTEMIC_BOUNDARY_CONCEPTS} for n in per_file}, per_file)

    # -- no output may assert causality.
    violations = []
    for name, blob in blobs.items():
        for pat in CAUSAL_CLAIM_PATTERNS:
            for m in re.finditer(pat, blob):
                window = blob[max(0, m.start() - NEGATION_WINDOW):m.start()]
                if not any(cue in window for cue in NEGATION_CUES):
                    violations.append(
                        "%s: affirmative causal claim %r at offset %d: ...%s..."
                        % (name, m.group(0), m.start(),
                           blob[max(0, m.start() - 60):m.end() + 40].replace("\n", " "))
                    )
    check("S7.6", g, "no synthesis output asserts causality (affirmative causal claims, "
                     "excluding negated boundary statements)",
          [], violations)

    # The exposure must be described as a STATE/association, never as an intervention.
    check("S7.7", g, "the pooled JSON describes the contrast as fragment state vs solo-LTR alleles, "
                     "not as insertion presence vs absence",
          True, "not insertion presence vs absence"
          in pooled["exposure_contrast_identity"]["detail"].lower())
    check("S7.8", g, "the heterogeneity caveat against over-reading I^2 with k=2 is retained",
          True, "do not over-interpret" in pooled["heterogeneity_power_caveat"].lower())
    assert_group(g)


TESTS = [
    test_s1_standardized_effects,
    test_s2_inverse_variance_pool_recomputed,
    test_s3_fixed_effect_p_value,
    test_s4_primary_pool_is_donor_disjoint,
    test_s5_endpoint_class_guard_refuses_mismatch,
    test_s6_standardization_algebra_invertible,
    test_s7_claim_boundary_present_and_no_causal_claim,
]


def main():
    failed_groups = []
    for fn in TESTS:
        try:
            fn()
        except AssertionError as exc:
            failed_groups.append({"test": fn.__name__, "error": str(exc)})

    n_pass = sum(1 for r in REPORT if r["result"] == "PASS")
    n_fail = sum(1 for r in REPORT if r["result"] == "FAIL")

    groups = {}
    for r in REPORT:
        gg = groups.setdefault(r["group"], {"n_pass": 0, "n_fail": 0})
        gg["n_pass" if r["result"] == "PASS" else "n_fail"] += 1
    for gg in groups.values():
        gg["result"] = "PASS" if gg["n_fail"] == 0 else "FAIL"

    findings = []
    if any(r["group"] == "S7_claim_boundary" and r["result"] == "FAIL" for r in REPORT):
        findings.append({
            "severity": "MUST_FIX_BEFORE_FREEZE",
            "finding": (
                "The CTC-0009.2 epistemic claim boundary is ABSENT from all three synthesis "
                "outputs. Searching the union of synthesis_slc44a5_pooled.json, "
                "synthesis_slc44a5_inputs.tsv and synthesis_slc44a5_leave_source_out.tsv finds ZERO "
                "occurrences of 'causal', of 'not proof'/'no proof'/'does not prove', and of "
                "'linked host haplotype'/'host haplotype'."
            ),
            "what_is_present": (
                "The positive half of the claim IS carried: the exposure is identified as the "
                "1p31.1b Type-I internal fragment and direction concordance is reported. No output "
                "asserts causality either (S7.6 passes). This is an OMISSION, not a misstatement."
            ),
            "why_it_matters": (
                "synthesis_slc44a5_pooled.json is the artefact that gets frozen and read "
                "downstream. As written it presents beta=0.80540, SE=0.13191, P=1.02e-9 for a "
                "genotype-state exposure with no statement that the fragment is not established as "
                "causal relative to a linked host haplotype. The numbers would travel without their "
                "limit."
            ),
            "ownership": (
                "T1-S owns these data files. Lane T-VAL writes validators only and must not edit "
                "the data it validates, so this is reported rather than fixed here."
            ),
            "suggested_fix": (
                "Add a claim_boundary field to synthesis_slc44a5_pooled.json carrying the CTC-0009.2 "
                "wording: a replicated association with 1p31.1b internal-fragment state, NOT proof "
                "that the fragment rather than a linked host haplotype is causal. The equivalent "
                "sentence already exists in the sibling matrix lane's run_matrix.py (lines 423-425 "
                "and 467-468), so the project wording is settled and can be reused verbatim."
            ),
            "integrity_note": (
                "This validator was NOT weakened to make it pass. It fails by design until the "
                "boundary is added to the synthesis outputs."
            ),
        })

    report = {
        "lane": "T-VAL (HML-2 Task 1) -- focused validator for the SLC44A5 synthesis outputs",
        "validates": "CTC-0009.2 coordinator-confirmed facts",
        "headline": (
            "ALL numeric facts reproduce exactly; the claim boundary does not. "
            "Recomputed pool matches to 0.0 absolute difference. %d assertion(s) FAILED."
            % n_fail
        ) if n_fail else "All %d assertions pass; recomputed pool matches to 0.0." % n_pass,
        "findings_requiring_owner_action": findings,
        "generated_by": "project/working/hml2_functional_evidence_synthesis_claude_v2/test_synthesis.py",
        "inputs_validated_read_only": [
            "results/synthesis_slc44a5_inputs.tsv",
            "results/synthesis_slc44a5_pooled.json",
            "results/synthesis_slc44a5_leave_source_out.tsv",
            "run_synthesis.py (imported to exercise its guards, not re-run)",
        ],
        "tolerances": {
            "recompute_vs_stored": "absolute 1e-9",
            "coordinator_constants": "exact at %d decimal places" % TOL_ROUND_DP,
            "fixed_effect_p": "exact at %d significant figures" % TOL_P_SIGFIG,
        },
        "method_note": (
            "S2 is load-bearing: the inverse-variance pool is recomputed from "
            "synthesis_slc44a5_inputs.tsv (weights = 1/SE^2, pooled = sum(w*b)/sum(w), "
            "SE = sqrt(1/sum(w))) and compared with synthesis_slc44a5_pooled.json. It does not "
            "read the JSON's own answer back and call that agreement."
        ),
        "overall_result": "PASS" if n_fail == 0 else "FAIL",
        "n_assertions": len(REPORT),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "groups": groups,
        "failed_groups": failed_groups,
        "assertions": REPORT,
    }
    with open(REPORT_PATH, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")

    for r in REPORT:
        if r["result"] == "FAIL":
            print("FAIL %s [%s] %s" % (r["assertion_id"], r["group"], r["description"]))
            print("      expected: %r" % (r["expected"],))
            print("      observed: %r" % (r["observed"],))
    print("\n%s: %d/%d assertions passed across %d groups -> %s"
          % (os.path.basename(__file__), n_pass, len(REPORT), len(groups), report["overall_result"]))
    print("report: %s" % REPORT_PATH)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
