#!/usr/bin/env python3
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "project/working/direct_ebv_fitness_screen_v1/results"
SEARCH = ROOT / "project/working/data_search_direct_cellular_phenotypes_v1/derived"


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(value):
    return float(value) if value not in {"", "NA"} else math.nan


summary = {r["outcome"]: r for r in rows(RESULTS / "secondary_outcome_summary.tsv")}
assert {k: tuple(int(float(v[f])) for f in ("people", "exposures", "unique_vectors",
                                                    "fitted_representatives", "nominal_p_lt_0.05",
                                                    "outcome_bh_lt_0.05"))
        for k, v in summary.items()} == {
    "Houldcroft2014_EBV_qPCR": (49, 79, 68, 49, 5, 1),
    "Im2012_intrinsic_growth": (34, 79, 65, 42, 5, 2),
}

h_source = rows(SEARCH / "Houldcroft2014_EBV_qPCR_unique_direct292_join.tsv")
im_source = rows(SEARCH / "Im2012_LCL_intrinsic_growth_direct292_join.tsv")
assert len(h_source) == 49 and len({r["sample"] for r in h_source}) == 49
assert len(im_source) == 34 and len({r["sample"] for r in im_source}) == 34
duplicates = {r["sample"]: r for r in h_source if int(r["n_source_rows"]) > 1}
assert set(duplicates) == {"NA19238", "NA19682", "NA19776"}
for row in duplicates.values():
    source_values = [float(x) for x in row["source_values"].split("|")]
    assert math.isclose(float(row["relative_ebv_copy_number_qpcr_mean"]),
                        sum(source_values) / len(source_values), rel_tol=1e-8)

h_matrix = rows(RESULTS / "Houldcroft2014_EBV_qPCR/person_level_outcome_matrix.tsv")
im_matrix = rows(RESULTS / "Im2012_intrinsic_growth/person_level_outcome_matrix.tsv")
assert len(h_matrix) == 49 and len(im_matrix) == 34
assert all(math.isclose(float(r["log2_relative_ebv_copy_number_qpcr"]),
                        math.log2(float(r["relative_ebv_copy_number_qpcr_mean"])), rel_tol=1e-10)
           for r in h_matrix)
assert all(math.isclose(float(r["intrinsic_growth_rate_per_10000"]),
                        float(r["intrinsic_growth_rate"]) / 10000, rel_tol=1e-10)
           for r in im_matrix)

secondary = rows(RESULTS / "secondary_outcome_model_results.tsv")
assert len(secondary) == 158
by_outcome = defaultdict(list)
for row in secondary:
    by_outcome[row["outcome"]].append(row)
for outcome_rows in by_outcome.values():
    groups = defaultdict(list)
    for row in outcome_rows:
        groups[row["observed_vector_group"]].append(row)
    assert all(sum(r["multiplicity_representative"] == "TRUE" for r in group) == 1
               for group in groups.values())

h_np9 = next(r for r in secondary if r["outcome"] == "Houldcroft2014_EBV_qPCR"
             and r["exposure_id"] == "general::type1_compatible_np9_units")
assert math.isclose(numeric(h_np9["beta"]), 0.0432632088, rel_tol=1e-7)
assert numeric(h_np9["p_pedigree_cluster"]) < 0.001
assert numeric(h_np9["population_p_cluster"]) < 0.015
assert numeric(h_np9["q_bh_outcome_global"]) < 0.05
assert numeric(h_np9["q_bh_secondary_suite"]) < 0.05

all_three = rows(RESULTS / "all_three_outcome_multiplicity.tsv")
assert len(all_three) == 156
three_key = {(r["outcome"], r["exposure_id"]): r for r in all_three}
h_np9_three = three_key[("Houldcroft2014_EBV_qPCR", "general::type1_compatible_np9_units")]
im_six_three = three_key[("Im2012_intrinsic_growth",
                          "HML-2_6q14.1::structural::solo_ltr_vs_any_provirus")]
assert 0.05 < numeric(h_np9_three["q_bh_three_outcome_suite"]) < 0.07
assert numeric(im_six_three["q_bh_three_outcome_suite"]) < 0.05
assert numeric(im_six_three["population_p_cluster"]) < 0.001

np9_loo = rows(RESULTS / "Houldcroft2014_EBV_qPCR/type1_np9_leave_one_person_out.tsv")
assert len(np9_loo) == 47
assert all(numeric(r["beta_superpopulation"]) > 0 for r in np9_loo)
assert max(numeric(r["p_superpopulation"]) for r in np9_loo) < 0.004
assert max(numeric(r["p_population"]) for r in np9_loo) < 0.05

np9_decomp = rows(RESULTS / "Houldcroft2014_EBV_qPCR/type1_np9_locus_decomposition.tsv")
assert len(np9_decomp) == 68
leave_primary = {r["label"]: r for r in np9_decomp
                 if r["mode"] == "leave_locus_out" and r["ancestry_model"] == "superpopulation"}
assert len(leave_primary) == 17
assert numeric(leave_primary["HML-2_5q33.3"]["p_pedigree_cluster"]) > 0.1
assert numeric(leave_primary["HML-2_19p12c"]["p_pedigree_cluster"]) < 0.04

six_loo = rows(RESULTS / "Im2012_intrinsic_growth/growth_6q14_leave_one_person_out.tsv")
assert len(six_loo) == 34
assert sum(float(r["dropped_exposure"]) > 0 for r in six_loo) == 7
assert all(numeric(r["beta_superpopulation"]) < 0 for r in six_loo)
assert max(numeric(r["p_superpopulation"]) for r in six_loo) < 0.002
assert max(numeric(r["p_population"]) for r in six_loo) < 0.001

groups = rows(RESULTS / "Im2012_intrinsic_growth/growth_6q14_group_summary.tsv")
carrier_groups = [r for r in groups if float(r["exposure"]) > 0]
assert sum(int(float(r["n"])) for r in carrier_groups) == 7
assert {r["population"] for r in carrier_groups} == {"ASW", "CHB", "JPT", "YRI"}

burden = rows(RESULTS / "Im2012_intrinsic_growth/type2_burden_6q14_exclusion_audit.tsv")
without = [r for r in burden if r["label"] == "type2_physical_units_without_6q14"]
assert len(without) == 2
assert all(numeric(r["beta"]) > 0 for r in without)
assert all(numeric(r["p_pedigree_cluster"]) < 0.012 for r in without)

print("direct_ebv_fitness_screen_v1 secondary validation PASS")
