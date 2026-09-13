#!/usr/bin/env python3

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
WORK = ROOT / "project/working/direct_geuvadis_function_screen_v1"
DATA = WORK / "data"
RESULTS = WORK / "results"


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(value):
    return float(value) if value not in {"", "NA"} else math.nan


materialization = json.loads((DATA / "materialization_summary.json").read_text())
assert materialization == {
    "expression_people": 462,
    "sdrf_people": 462,
    "direct_people": 292,
    "exact_direct_overlap": 33,
    "expression_genes": 23722,
    "gencode_mapped_genes": 22416,
    "populations": ["FIN", "GBR", "TSI", "YRI"],
    "superpopulations": ["AFR", "EUR"],
    "females": 22,
    "males": 11,
    "repeated_pedigree_components": 0,
    "prior_MAGE_overlap_people": 5,
    "GEUVADIS_only_people": 28,
    "frozen_exposures": 79,
    "fixed_type1_loci": 17,
}

identity = rows(DATA / "identity_audit.tsv")
assert len(identity) == 33 and len({r["sample"] for r in identity}) == 33
assert all(r["expression_header_exact"] == r["sdrf_source_exact"] == r["sdrf_individual_exact"] == "TRUE"
           for r in identity)
assert all(r["official_sex"] == r["direct_metadata_sex"] for r in identity)
assert Counter(r["population"] for r in identity) == {"FIN": 11, "GBR": 10, "TSI": 7, "YRI": 5}
assert Counter(r["official_sex"] for r in identity) == {"female": 22, "male": 11}
assert len({r["pedigree_component"] for r in identity}) == 33

exposures = rows(DATA / "direct33_exposures.tsv")
assert len(exposures) == 33 and len({r["sample"] for r in exposures}) == 33
assert Counter(r["HML-2_1p31.1b::internal_fragment_present"] for r in exposures) == {"0": 16, "1": 17}
assert Counter(r["HML-2_1p31.1b::total_internal_fragment_cn"] for r in exposures) == {
    "0": 16, "1": 13, "2": 3, "5": 1
}
assert Counter(r["HML-2_7p22.1::total_physical_cn"] for r in exposures) == {"2": 16, "3": 14, "4": 3}
for field in ("general::type1_physical_units", "general::type2_physical_units"):
    assert all(r[field] not in {"", "NA"} for r in exposures)

type1 = rows(DATA / "type1_locus_units.tsv")
assert len(type1) == 33 and len(type1[0]) == 18
exposure_by_id = {r["sample"]: r for r in exposures}
for row in type1:
    assert sum(int(v) for k, v in row.items() if k != "sample") == int(
        exposure_by_id[row["sample"]]["general::type1_physical_units"]
    )

expression_path = DATA / "direct33_expression.tsv"
slc = None
gene_count = 0
with expression_path.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        gene_count += 1
        if row["gene_id_base"] == "ENSG00000137968":
            slc = row
assert gene_count == 23722 and slc is not None
assert slc["gene_symbol"] == "SLC44A5" and slc["chromosome_grch38"] == "chr1"
assert all(math.isfinite(float(slc[r["sample"]])) for r in identity)

run_summary = {r["metric"]: int(float(r["value"])) for r in rows(RESULTS / "run_summary.tsv")}
assert run_summary == {
    "official_expression_people": 462,
    "exact_direct_people": 33,
    "genes": 23722,
    "gencode_mapped_genes": 22416,
    "prior_MAGE_overlap_people": 5,
    "GEUVADIS_only_people": 28,
    "frozen_exposures": 79,
    "unique_direct33_vectors": 64,
    "fitted_vector_groups": 35,
    "broad_unique_tests": 403666,
    "broad_nominal_p_lt_0.05": 12780,
    "broad_global_bh_lt_0.05": 0,
    "broad_global_holm_lt_0.05": 0,
}

alias = rows(RESULTS / "exposure_alias_map.tsv")
assert len(alias) == 79
groups = defaultdict(list)
for row in alias:
    groups[row["observed_direct33_vector_group"]].append(row)
assert len(groups) == 64
assert all(sum(r["direct33_multiplicity_representative"] == "TRUE" for r in group) == 1
           for group in groups.values())

targeted = rows(RESULTS / "slc44a5_targeted_results.tsv")
def target(model_id, ancestry="population"):
    matches = [r for r in targeted if r["model_id"] == model_id and r["ancestry_model"] == ancestry]
    assert len(matches) == 1
    return matches[0]

replication = target("prespecified_MAGE_replication")
assert math.isclose(numeric(replication["beta"]), 0.609216023, rel_tol=1e-7)
assert numeric(replication["p_two_sided_hc3"]) < 1.6e-5
assert numeric(replication["p_one_sided_positive"]) < 8e-6
assert numeric(target("prespecified_MAGE_replication", "superpopulation")["p_two_sided_hc3"]) < 1e-6
independent = target("MAGE_nonoverlap_replication")
assert int(independent["n"]) == 28
assert numeric(independent["beta"]) > 0.51 and numeric(independent["p_two_sided_hc3"]) < 0.0013
assert numeric(target("general_type1_physical_burden")["p_two_sided_hc3"]) > 0.98
local_adjusted = target("local_state_adjusted_type1_burden")
assert numeric(local_adjusted["beta"]) > 0.60 and numeric(local_adjusted["p_two_sided_hc3"]) < 2.1e-5
assert numeric(target("type1_burden_adjusted_local_state")["p_two_sided_hc3"]) > 0.8
assert 0.04 < numeric(target("onep31b_feature::HML-2_1p31.1b::total_internal_fragment_cn")["p_two_sided_hc3"]) < 0.06

locus_models = [r for r in targeted if r["family"] == "type1_locus_heterogeneity"
                and r["ancestry_model"] == "population" and r["status"] == "fit"]
assert len(locus_models) == 2
assert min(numeric(r["p_two_sided_hc3"]) for r in locus_models) > 0.4

loo = rows(RESULTS / "slc44a5_leave_one_person_out.tsv")
assert len(loo) == 33
assert min(numeric(r["beta_population"]) for r in loo) > 0.57
assert max(numeric(r["beta_population"]) for r in loo) < 0.70
assert max(numeric(r["p_population_hc3"]) for r in loo) < 1.4e-4

permutation = rows(RESULTS / "slc44a5_population_blocked_permutation.tsv")
assert len(permutation) == 1 and permutation[0]["permutations"] == "50000"
assert numeric(permutation[0]["empirical_p_two_sided"]) < 5e-5

inverse = rows(RESULTS / "slc44a5_inverse_normal_sensitivity.tsv")
assert len(inverse) == 2 and all(numeric(r["beta"]) > 1.4 for r in inverse)
assert all(numeric(r["p_two_sided_hc3"]) < 2.1e-5 for r in inverse)

# Stream the large screen once: global multiplicity, focal cis results, and the
# exact prespecified row all remain linked to the same materialized statistics.
broad_count = nominal = bh_hits = holm_hits = 0
slc_broad = iqgap3 = dpm3 = None
sevenp_q = []
with (RESULTS / "broad_model_results.tsv").open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        broad_count += 1
        nominal += numeric(row["p_population_sex_hc3"]) < 0.05
        bh_hits += numeric(row["q_bh_global"]) < 0.05
        holm_hits += numeric(row["p_holm_global"]) < 0.05
        if row["exposure_id"] == "HML-2_1p31.1b::internal_fragment_present" and row["gene_id"] == "ENSG00000137968":
            slc_broad = row
        if row["exposure_id"] == "HML-2_1q22::orf::orf_gag" and row["gene_symbol"] == "IQGAP3":
            iqgap3 = row
        if row["exposure_id"] == "HML-2_1q22::orf::orf_gag" and row["gene_symbol"] == "DPM3":
            dpm3 = row
        if row["exposure_id"].startswith("HML-2_7p22.1::"):
            sevenp_q.append(numeric(row["q_bh_within_exposure"]))
assert (broad_count, nominal, bh_hits, holm_hits) == (403666, 12780, 0, 0)
assert slc_broad is not None and numeric(slc_broad["q_bh_global"]) > 0.49
assert numeric(slc_broad["q_bh_lane"]) < 0.007
assert numeric(slc_broad["q_bh_within_exposure"]) < 0.0002
assert iqgap3 is not None and numeric(iqgap3["p_population_sex_hc3"]) > 0.8
assert dpm3 is not None and numeric(dpm3["p_population_sex_hc3"]) < 0.002
assert numeric(dpm3["q_bh_within_exposure"]) > 0.1
assert sevenp_q and min(sevenp_q) > 0.6

print("direct_geuvadis_function_screen_v1 validation PASS")
