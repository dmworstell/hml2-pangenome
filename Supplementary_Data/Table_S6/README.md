# Table S3 and main Table 2

The reconstruction comparison includes 282 matched donors and 83 shared loci after catalog filtering. It is distinct from the published-VCF structural comparison in Figure S1, whose source remains `../Table_S6_short_read_vcf_locus_diagnostics.tsv`.

`results/summary_metrics.tsv` supplies coding and structural recovery, callability and confidence intervals. `results/per_locus_coding_recovery.tsv` supplies the locus-level results. `derived/` records the catalog universe and original source provenance.

The two compressed input tables retain every field consumed by the comparison. Long-read rows carry their original `analysis_include` flag. The loader excludes rows unless that flag is 1. The short-read projection retains only public donor IDs and loci represented in the retained long-read catalog. No sequence or private collaborator data are added. `input_projection.json` identifies the original catalogs and their exact projections.

From the repository root, run

```sh
python manuscript_figures/python/analysis/longread_shortread_comparison/run_analysis.py \
  --long-read Supplementary_Data/Table_S6/inputs/long_read.tsv.gz \
  --short-read Supplementary_Data/Table_S6/inputs/short_read.tsv.gz \
  --output-dir outputs/short_read_comparison --bootstrap 2000 --skip-plots
```

The results count one donor per locus and feature. Linked Gag–Pro–Pol annotations must occur on the same provirus. Missing data remain unknown. Env denotes the screened Env-region ORF, including the truncated Type-I frame. These are results for the tested short-read reconstruction workflow.

Main Table 2 reports the seven `overall_workflow_recovery` rows as n/N (%). N includes all long-read-positive calls. The optional plot contains aggregate recovery counts only. No binary locus heat map is used in the manuscript.
