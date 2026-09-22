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

The results count one donor per locus and feature. Linked Gag–Pro–Pol annotations must occur on the same provirus. Missing data remain unknown. The original pooled Env results combine Type-II Env and the theoretical N-terminally truncated Type-I Env annotation. Translation of that Type-I product has not been demonstrated. Structural recovery refers to the catalog annotations, not direct recovery of every symbolic structural-variant record in the VCF.

Main Table 2 now reports separate intact and combined ORF screens, and separates Type-II Env from theoretical Type-I Env. Its structural rows are unchanged. The results are in `Table_2_reanalysis/Table_2_ORF_recovery.tsv`; per-locus counts and the contributing donor–locus observations are supplied alongside it. N includes all long-read-positive calls under the stated criterion, and n counts those also positive in the VCF-derived catalog under the same criterion. The seven original aggregate totals are reproduced exactly; the original pooled Env results remain available in the data.

The intact screen requires the caller's `Intact` classification. The combined screen also retains its stop-free altered-frame candidates. Each gene was scored in its annotated frame; the canonical programmed gag/pro/pol frame transitions were not counted as disruptive sequence mutations. Type-I Env intactness refers only to the residual reading frame, not a full-length Env protein.

Run the new table calculation from any directory after extracting the archive:

```sh
python /path/to/Supplementary_Data/Table_S6/Table_2_reanalysis/recompute_table2.py \
  --data /path/to/Supplementary_Data \
  --out /path/to/reproduced_table_2
```

This script uses only the Python standard library. It reads the two frozen comparison inputs and the full retained catalog for the locus type. Every eligible locus has a single retained type; a mixed-type locus causes the script to stop rather than silently pool types. `provenance.json` records input hashes and the validation results. This Table 2 extension is included in the v0.4.2 code and derived-data archives. The frozen comparison inputs and the previous aggregate results are unchanged.

Figure S1 uses `../Table_S6_short_read_vcf_locus_diagnostics.tsv`. Its plotted input is the Illumina ensemble callset `1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz`, not a pool of the four VCF resources examined in the broader work. The figure builder reads the retained direction-corrected E2/E5 diagnostics from `project/working/short_read_direction_corrected_v4/results/bio5_direction_corrected_v4/short_read_numeric_authority_v1.json`. These diagnostics retain records compatible with the direction and reference state of the proposed structural change. Errors are averaged within each donor–locus combination before calculating the locus-level fraction. Each panel shows the 15 loci with the highest disagreement among those with at least ten matched combinations.
