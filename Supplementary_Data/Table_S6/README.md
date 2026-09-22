# Table S3 and main Table 2

Main Table 2 and Table S3 compare coding and structural annotations derived from the phased high-coverage 1000 Genomes VCF panel with long-read calls in 282 matched donors at 83 shared loci after catalog filtering. The historical upstream workflow applied each donor's phased variants to hg38 with `bcftools consensus`, then annotated the resulting haplotype sequences against KCON. Those sequences were VCF-derived intermediates, not a separate assembled short-read dataset.

The source panel was the `20220422_3202` release, with chromosome files named `1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel*.vcf.gz`. The comparison below consumes the retained ORF/structure catalog from that workflow. It does not reread raw sequencing data or regenerate haplotype FASTAs.

`results/summary_metrics.tsv` supplies coding and structural recovery, callability and confidence intervals. `results/per_locus_coding_recovery.tsv` supplies the locus-level results. `derived/` records the catalog universe and original source provenance.

The two compressed input tables retain every field consumed by the comparison. Long-read rows carry their original `analysis_include` flag. The loader excludes rows unless that flag is 1. The short-read projection retains only public donor IDs and loci represented in the retained long-read catalog. No sequence or private collaborator data are added. `input_projection.json` identifies the original catalogs, their exact projections and the short-read source callset.

The original short-read input was the retained project-data file `1kgp_comparison/orf_analysis/combined_1kgp_hml2_orf_analysis.tsv`, 1,030,110,706 bytes, with SHA256 `4cf6c2f579356be88f34d27e4aaaff797b2d97b108bffc0289afa0c334922872`. A later same-named catalog in `HML2_ProjectResources/data/catalog/` is not this input. Use the bundled `inputs/short_read.tsv.gz` to reproduce these results. The full original VCFs and intermediate haplotype FASTAs are not bundled here.

From the repository root, run

```sh
python manuscript_figures/python/analysis/longread_shortread_comparison/run_analysis.py \
  --long-read Supplementary_Data/Table_S6/inputs/long_read.tsv.gz \
  --short-read Supplementary_Data/Table_S6/inputs/short_read.tsv.gz \
  --output-dir outputs/short_read_comparison --bootstrap 2000 --skip-plots
```

The results count one donor per locus and feature. Linked Gag–Pro–Pol annotations must occur on the same provirus. Missing data remain unknown. Env denotes the screened Env-region ORF, including the truncated Type-I frame. Structural recovery refers to the catalog annotations, not direct recovery of every symbolic structural-variant record in the VCF.

Main Table 2 reports the seven `overall_workflow_recovery` rows as n/N (%). N includes all long-read-positive calls. The optional plot contains aggregate recovery counts only. No binary locus heat map is used in the manuscript.

Figure S1 uses `../Table_S6_short_read_vcf_locus_diagnostics.tsv`. Its plotted input is the Illumina ensemble callset `1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz`, not a pool of the four VCF resources examined in the broader work. The figure builder reads the retained direction-corrected E2/E5 diagnostics from `project/working/short_read_direction_corrected_v4/results/bio5_direction_corrected_v4/short_read_numeric_authority_v1.json`. These diagnostics retain records compatible with the direction and reference state of the proposed structural change. Errors are averaged within each donor–locus combination before calculating the locus-level fraction. Each panel shows the 15 loci with the highest disagreement among those with at least ten matched combinations.
