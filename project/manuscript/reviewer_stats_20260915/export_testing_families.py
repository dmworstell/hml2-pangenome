#!/usr/bin/env python3
"""Export the full testing families behind the manuscript's reported q values.

Reads retained analysis outputs without refitting models or replacing missing
P values. Outputs machine-readable tables for the supplementary-data package.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "exports"


def bh(p):
    p = pd.to_numeric(p, errors="raise")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    finite = p[np.isfinite(p)].sort_values(kind="stable")
    if len(finite):
        assert finite.between(0, 1).all()
        values = finite.to_numpy() * len(finite) / np.arange(1, len(finite) + 1)
        out.loc[finite.index] = np.minimum(1, np.minimum.accumulate(values[::-1])[::-1])
    return out


def verify(p, recorded):
    calculated = bh(p)
    assert calculated.isna().equals(recorded.isna())
    delta = (calculated - recorded).abs().max()
    assert delta < 1e-11, delta
    return float(delta)


def read(relative):
    return pd.read_csv(PROJECT / relative, sep="\t", na_values=["NA"])


def save(frame, filename):
    frame.to_csv(OUT / filename, sep="\t", index=False, na_rep="NA",
                 compression="infer")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    refit = read("manuscript/artifact_filtered_functional_refit/complete_three_outcome_refit_multiplicity.tsv")
    old = pd.concat([
        read("working/direct_ebv_fitness_screen_v1/results/model_results.tsv").assign(outcome="Mandage2017_EBV_in_silico"),
        read("working/direct_ebv_fitness_screen_v1/results/secondary_outcome_model_results.tsv"),
    ], ignore_index=True)
    original_q = "q_bh_conservative_237_model_suite"
    if original_q not in refit:
        original_q = "q_bh_finite_model_suite"
    delta = verify(refit.analysis_p_pedigree_cluster, refit[original_q])
    original_cols = ["outcome", "exposure_id", "model_status", "n", "adjustment"]
    refit = refit.drop(columns=["analysis_model_status", "analysis_n", "included_in_bh_family",
                               "bh_family_size", "bh_exclusion_reason"], errors="ignore")
    source = old[original_cols].rename(columns={"model_status": "original_model_status", "n": "original_n"})
    final = refit.merge(source, on=["outcome", "exposure_id"], validate="one_to_one")
    final["analysis_model_status"] = final.refit_model_status.fillna(final.original_model_status)
    final["analysis_n"] = np.where(final.was_refit, final.refit_n, final.original_n).astype(int)
    final["included_in_bh_family"] = np.isfinite(final.analysis_p_pedigree_cluster)
    final["bh_family_size"] = int(final.included_in_bh_family.sum())
    final["bh_exclusion_reason"] = np.where(final.included_in_bh_family, "", final.analysis_model_status)
    final = final.rename(columns={original_q: "q_bh_finite_model_suite"})
    final = final.drop(columns=["original_model_status", "original_n"])
    save(final, "Table_S10a_functional_237_model_results.tsv")
    burden = read("manuscript/artifact_filtered_functional_refit/artifact_filtered_general_burden_model_results.tsv")
    verify(burden.p_pedigree_cluster, burden.q_bh_51_model_suite)
    save(burden, "Table_S10b_functional_51_refitted_models.tsv")

    mage_base = "working/type1_functional_evidence_synthesis_v1/results/"
    mage_cols = ["truth_family_id", "representative_alias", "locus", "endpoint",
                 "gene_id", "gene_symbol", "n", "beta_per_exposure_sd", "se",
                 "statistic", "p_value", "residual_df", "bh_global_targeted_family"]
    mage = pd.read_csv(PROJECT / (mage_base + "mage_targeted_locus_endpoint_results.tsv.gz"),
                       sep="\t", usecols=mage_cols)
    mage_delta = verify(mage.p_value, mage.bh_global_targeted_family)
    assert not mage.duplicated(["truth_family_id", "gene_id"]).any()
    save(mage, "Table_S10e_MAGE_complete_discovery_family.tsv.gz")
    catalog = read(mage_base + "mage_targeted_locus_endpoint_catalog.tsv")
    save(catalog, "Table_S10f_MAGE_candidate_eligibility_and_aliases.tsv")

    geu = read("working/direct_geuvadis_function_screen_v1/results/slc44a5_targeted_results.tsv")
    geu["included_in_bh_family"] = ((geu.status == "fit") &
        (geu.ancestry_model == "population") & (geu.family != "prespecified_replication"))
    mask = geu.included_in_bh_family
    geu_delta = verify(geu.loc[mask, "p_two_sided_hc3"], geu.loc[mask, "q_bh_secondary_family"])
    assert geu.loc[~mask, "q_bh_secondary_family"].isna().all()
    geu["bh_family_size"] = int(mask.sum())
    geu["bh_exclusion_reason"] = np.select(
        [geu.ancestry_model != "population", geu.family == "prespecified_replication", geu.status != "fit"],
        ["superpopulation_sensitivity_not_in_population_BH_family", "primary_replication_not_in_secondary_BH_family", geu.status],
        default="")
    save(geu, "Table_S10g_GEUVADIS_complete_SLC44A5_followup_family.tsv")
    mage_hc3 = read("working/onep31b_slc44a5_followup_v1/robust_stats/results/model_sensitivities.tsv")
    save(mage_hc3, "Table_S10h_MAGE_SLC44A5_HC3_sensitivity_models.tsv")

    status_counts = final.groupby(["outcome", "analysis_model_status"]).size()
    summary = {
        "functional": {
            "attempted": len(final), "finite_p": int(final.included_in_bh_family.sum()),
            "missing_p": int((~final.included_in_bh_family).sum()),
            "status_counts": {" | ".join(key): int(value) for key, value in status_counts.items()},
            "max_q_difference": delta,
            "significant_q_rows": int((final.q_bh_finite_model_suite < 0.05).sum()),
            "deduplication": "No exact-vector deduplication in final combined correction. All finite exposure-outcome rows retained.",
        },
        "MAGE_discovery": {
            "candidate_aliases": len(catalog), "eligible_aliases": int((catalog.eligibility_status == "eligible").sum()),
            "unique_exposure_adjustment_vectors": int(mage.truth_family_id.nunique()),
            "genes": int(mage.gene_id.nunique()), "tests": len(mage),
            "max_q_difference": mage_delta,
            "significant_q_rows": int((mage.bh_global_targeted_family < 0.05).sum()),
            "SLC44A5_lead": mage.loc[mage.p_value.idxmin()].to_dict(),
        },
        "GEUVADIS_followup": {
            "rows": len(geu), "models": int(geu.model_id.nunique()),
            "corrected_population_models": int(mask.sum()), "max_q_difference": geu_delta,
            "significant_q_rows": int((geu.q_bh_secondary_family < 0.05).sum()),
            "disjoint_replication": geu.loc[(geu.model_id == "MAGE_nonoverlap_replication") &
                (geu.ancestry_model == "population")].iloc[0].to_dict(),
        },
    }
    (OUT.parent / "verification.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
