#!/usr/bin/env python3
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "project/working/direct_ebv_fitness_screen_v1/results"


def rows(name):
    with (RESULTS / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(value):
    return float(value) if value not in {"", "NA"} else math.nan


summary = {r["metric"]: int(float(r["value"])) for r in rows("run_summary.tsv")}
assert summary == {
    "source_phenotype_people": 1753,
    "direct_sequence_people": 292,
    "direct_phenotype_overlap": 116,
    "overlap_superpopulations": 5,
    "overlap_populations": 16,
    "related_pedigree_components": 2,
    "exposures_total": 79,
    "unique_direct_truth_vectors": 70,
    "models_fit": 65,
    "global_bh_lt_0.05": 0,
    "global_holm_lt_0.05": 0,
}

matrix = rows("person_level_direct_matrix.tsv")
assert len(matrix) == 292
assert len({r["sample"] for r in matrix}) == 292
overlap = [r for r in matrix if r["ebv_load"] != "NA"]
assert len(overlap) == 116
assert all(float(r["ebv_load"]) > 0 for r in overlap)
assert all(abs(float(r["log2_ebv_load"]) - math.log2(float(r["ebv_load"]))) < 1e-10 for r in overlap)
assert set(r["superpopulation"] for r in overlap) == {"AFR", "AMR", "EAS", "EUR", "SAS"}
component_counts = Counter(r["pedigree_component"] for r in overlap)
assert sum(v > 1 for v in component_counts.values()) == 2

# Catalog-wide no-row semantics were implemented as a zero contribution, never NA.
for field in ("type1_units", "type2_units", "type1_loci_present", "type2_loci_present"):
    assert all(r[field] not in {"", "NA"} for r in matrix)

authority = rows("type_authority_audit.tsv")
assert len(authority) == 54
assert Counter(r["type_class"] for r in authority) == {"fixed_type1": 17, "fixed_type2": 37}
assert all(";" not in r["observed_types"] for r in authority)
absence = {r["rule"]: r["value"] for r in rows("absence_semantics_audit.tsv")}
assert absence["expected_locus_haplotype_cells"] == "57232"
assert absence["observed_catalog_cells"] == "57017"
assert absence["no_row_cells"] == "215"
assert absence["catalog_wide_burden_no_row_value"] == "0"
assert absence["internal_orf_no_row_value"] == "not_applicable"

models = rows("model_results.tsv")
assert len(models) == 79
fit = [r for r in models if r["model_status"] == "fit" and r["multiplicity_representative"] == "TRUE"]
assert len(fit) == 65
assert all(math.isfinite(numeric(r["p_pedigree_cluster"])) for r in fit)
assert all(0 <= numeric(r["q_bh_global"]) <= 1 for r in fit)
assert not any(numeric(r["q_bh_global"]) < 0.05 for r in fit)
assert all(int(float(r["pedigree_clusters"])) <= int(float(r["n"])) for r in fit)

aliases = rows("truth_vector_alias_map.tsv")
groups = defaultdict(list)
for row in aliases:
    groups[row["observed_vector_group"]].append(row)
assert len(groups) == 70
assert all(sum(r["multiplicity_representative"] == "TRUE" for r in group) == 1 for group in groups.values())

by_id = {r["exposure_id"]: r for r in models}
assert math.isclose(numeric(by_id["general::type1_physical_units"]["p_pedigree_cluster"]), 0.65620509, rel_tol=1e-6)
assert math.isclose(numeric(by_id["HML-2_1q22::orf::orf_gag"]["p_pedigree_cluster"]), 0.89821401, rel_tol=1e-6)
assert numeric(by_id["HML-2_3q27.2::structural::solo_ltr_vs_any_provirus"]["p_pedigree_cluster"]) < 0.05
assert numeric(by_id["HML-2_3q27.2::structural::solo_ltr_vs_any_provirus"]["population_p_cluster"]) > 0.2
assert numeric(by_id["HML-2_3q27.2::structural::solo_ltr_vs_any_provirus"]["q_bh_global"]) > 0.9
assert by_id["HML-2_7p22.1::any_haplotype_cn_ge3"]["model_status"] == "underpowered_arm_lt5"

decomp = rows("type_burden_leave_one_locus_out.tsv")
assert len(decomp) == 54
assert Counter(r["locus_type"] for r in decomp) == {"type1": 17, "type2": 37}

print("direct_ebv_fitness_screen_v1 validation PASS")
