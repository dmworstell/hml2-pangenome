# Manuscript figure and table map

This map uses the final manuscript order from 22 September 2026. It applies to
release 0.4.1 on GitHub and its matching Zenodo archive. Paths are relative to
the repository root. A source filename is not a manuscript figure number.
Historical names are retained to preserve links to the original inputs.

For the revised analyses, use
[`analysis/september2026_revision/README.md`](analysis/september2026_revision/README.md).
The older scripts in `project/manuscript/` and `manuscript_figures/` retain their
original figure names and layouts. They are analysis provenance, and running
all of them does not assemble the final manuscript. In particular, the older
Figure 4, Figure 5C and TSD builders do not reproduce the September 21 revision.
BioRender artwork and collaborator Fiber-seq source data are not redistributed.

## Main figures

| Final label | Content | Source data and analysis |
|---|---|---|
| Figure 1 | Structural variation | Panel A is licensed BioRender artwork. Panel B uses `Supplementary_Data/Figure_1_structural_source_data.tsv`, `Figure_1_structural_observation_cells.tsv` and `Figure_1_sex_chromosome_eligibility.tsv`; `project/manuscript/build_narrative_main_figures.py`, `build_figure_1`. |
| Figure 2 | Locus copy counts and within-array coding variation | `Supplementary_Data/Table_S7_array_copy_number_distributions.tsv`, `Table_S7_haplotype_copy_number_calls.tsv`, `Figure_2_within_array_variation.tsv`; `project/manuscript/build_narrative_main_figures.py`, `build_figure_2`. Nucleotide differences are separately shown in Figure S5. |
| Figure 3 | Regional phylogenies, sequence sharing and chromosome locations | Current results and rendering source are in `Supplementary_Data/Table_S14/phylogeny_review/`; host-flank evidence is in `Supplementary_Data/Table_S14/4q_acrocentric_relationship/`. The current panel B reports mean pairwise nucleotide differences, not counts of shared sequences. |
| Figure 4 | ORF annotations and person-level coding potential | `Supplementary_Data/Table_S8_artifact_filtered_orf_coding_potential.tsv`, `Figure_4_person_orf_burden.tsv`, `Figure_4_locus_orf_superpopulation_frequencies.tsv`; `project/manuscript/build_narrative_main_figures.py`, `build_figure_4`. |
| Figure 5 | Type-I cassette conservation and propagation models | Revised panel C uses `Supplementary_Data/Table_S17/window_summary.tsv` and `window_pairwise_differences.tsv`, including the LTR5Hs-only Type-II control. Panels A, B and D retain their existing inputs in `project/manuscript/build_narrative_main_figures.py`; model evidence is in `Supplementary_Data/Figure_5_source_effect_model_comparison.tsv`. |
| Figure 6 | Structural states and cellular phenotypes | `Supplementary_Data/Table_S10a_functional_237_model_results.tsv`, `Table_S10c_SLC44A5_disjoint_cohort_estimates.tsv`, `Table_S10d_anti_CD20_adjusted_state_models.tsv`, and `Functional_refit/`; `project/manuscript/build_narrative_main_figures.py`, `build_figure_6`. |
| Figure 7 | 8q11.23 provirus structure and solo-LTR haplotypes | `data/8q11_23_provirus_annotation.json`, `Supplementary_Data/Figure_7_solo_LTR_haplotype_counts.tsv`; `project/manuscript/retained_panels.py`, `draw_eightq_network`. Replot with `scripts/reproduce_compact_panels.py`. |

## Main tables

| Final label | Content | Source |
|---|---|---|
| Table 1 | Cohorts and denominators | Retained donor roster and the listed comparison cohorts. |
| Table 2 | Coding sequences and structural states recovered by short reads | `Supplementary_Data/Table_S6/results/summary_metrics.tsv`, rows for `overall_workflow_recovery`. Counts include all long-read-positive calls, not only callable short-read records. |

## Supplementary figures

| Final label | Content | Source data and analysis |
|---|---|---|
| Figure S1 | Matched short-read and long-read structural calls | `Supplementary_Data/Table_S6_short_read_vcf_locus_diagnostics.tsv`. This is the published-VCF comparison, distinct from the reconstruction benchmark in main Table 2. Historical plotting function `build_short_read_figure` in `project/manuscript/build_artifact_filtered_catalog_figures.py` writes a filename beginning `Figure_2`. |
| Figure S2 | Telomeric Type-II LTR and host-flank trees | `Supplementary_Data/Table_S14/telomeric_phylogeny/`, including `reference_LTRs.aln.fa.nwk`, `B_flank_tree.nwk`, `C_flank_tree.nwk`, bootstrap tables and `precision_figure.py`. |
| Figure S3 | Review of apparent extra copies | `Supplementary_Data/Table_S4_CNV_assembly_artifact_qc.tsv` and `Table_S4_NucFreq_region_summary.tsv`; `build_apparent_copy_review_figure` in `project/manuscript/build_artifact_filtered_catalog_figures.py` writes `Figure_S9_apparent_extra_copy_review`. |
| Figure S4 | Positional ORF retention and array sizes | `Supplementary_Data/Panel_data/tandem_resolved/`; `manuscript_figures/R/figures/Duplication_Analysis.R`. Replot the retained summaries with `Rscript scripts/reproduce_s4.R outputs/Figure_S4`. The original plotting blocks are internally labelled Figure 2 and Figure 3. |
| Figure S5 | Nucleotide differences within tandem arrays | `Supplementary_Data/Table_S15/`, especially pairwise differences, site frequencies and variable-site tables. Historical `figure_arrays.py` uses `Figure_S20_array_nucleotide_variation`. Use the September revision workflow for current reproduction. |
| Figure S6 | Candidate paired-TSD differences and length ambiguity | `Supplementary_Data/Figure_S7_boundary_evidence/`. These are the corrected terminal-junction measurements. `Supplementary_Data/Figure_S10_TSD_per_locus.tsv` and the old TSD builder are historical and do not supply the current panel. The unsupported 3q12.3 paired-TSD result is withdrawn. |
| Figure S7 | Gag and Pro phylogenies | `Supplementary_Data/Phylogeny/hml2_pan_orf_gag_tree.nwk`, `hml2_pan_orf_pro_tree.nwk`, `selected_sequence_clusters.tsv` and `selected_sequence_clusters.fasta`. |
| Figure S8 | LTR, Env and Pol phylogenies | `Supplementary_Data/Phylogeny/hml2_pan_ltr_expanded_tree.nwk`, `hml2_pan_orf_env_tree.nwk`, `hml2_pan_orf_pol_tree.nwk`, and their selected-cluster tables. The directory named `Supplementary_Data/Figure_S9/` instead belongs to current Figure S10. |
| Figure S9 | Protein substitutions, frameshifts and Pol Y195C | `Supplementary_Data/Panel_data/variant_positions_resolved/`; `manuscript_figures/R/figures/lollipop_plots.R`. |
| Figure S10 | ORF combinations and independent longest-frame scan | `Supplementary_Data/Figure_S9/`, including `highest_annotation_tiers.tsv`, `primary_longest_frame_scan.tsv`, `gag_region_longest_frame_witnesses.tsv` and `figure_counts_and_provenance.json`; current Python renderer `manuscript_figures/python/figures/atypical_fusion_figure.py`. `python scripts/reproduce_s10.py --output outputs/Figure_S10` checks the retained witness-table counts and replots them. It does not rescan sequences or replace the accepted annotation tiers. |
| Figure S11 | Direct Type-I cassette calls | `Supplementary_Data/Table_S11_direct_TypeI_locus_calls.tsv` and `Table_S11/`; `project/manuscript/retained_panels.py`, `draw_type_state_counts`; `scripts/reproduce_compact_panels.py`. |
| Figure S12 | Regional nearest-neighbor changes, linked cassette sites and deletion groups | `manuscript_figures/python/analysis/type1_cassette_mechanism/`, `project/manuscript/delta292_ape_counterfactuals_v1/`, and `project/manuscript/build_delta292_ape_counterfactuals.py`. The 83-cluster primate alignment is also in Table S11. This figure uses zero-based KCON positions 6331/6492 and deletion interval [6501,6793). |
| Figure S13 | Human cassette frequency logos and similarity ranking | `Supplementary_Data/Table_S17/human_cassette_site_frequencies.tsv`, `candidate_similarity_to_TypeI.tsv`, `human_representatives.tsv`; September revision workflow. This figure uses one-based KCON sites 6332/6493 and deletion interval 6502–6793. |
| Figure S14 | Alternative models for Δ292 enrichment | `project/manuscript/build_delta292_mechanism_discrimination_figure.py`; current output label supplied by `scripts/reproduce_compact_panels.py`. Main source-effect comparison remains Figure 5D. |
| Figure S15 | Helper dependence and neutral-founder drift | `project/manuscript/build_delta292_helper_abc_figure.py`; current output label supplied by `scripts/reproduce_compact_panels.py`. |
| Figure S16 | Gag dosage, cell growth and EBV abundance | `Supplementary_Data/Table_S9_1q22_Gag_artifact_corrected_truth.tsv`, `Table_S10a_functional_237_model_results.tsv`, `Table_S10b_functional_51_refitted_models.tsv`, and `Functional_refit/`; `project/manuscript/build_current_results_figure.py`. |
| Figure S17 | Fiber-seq chromatin accessibility | `manuscript_figures/R/analysis/Fiberseq_haplotype_reanalysis.R` and `Supplementary_Data/Figure_S17_locus_crosswalk.tsv`. The latter's S17 label is historical. Collaborator source data and artwork are not included in this public release. |
| Figure S18 | Solo-LTR diversity and illustrative substitution expectations | `Supplementary_Data/Table_S16/`; retained `make_Figure_S21.py` has a historical filename. Use the September revision workflow. This comparison does not estimate an insertion or coalescence date. |

## Supplementary tables and alignment

| Final label | Content | Public source |
|---|---|---|
| Table S1 | Structural-analysis catalog | `Supplementary_Data/Table_S1_structural_observations_corrected.tsv`; complete corrected catalog is restored from `data/` as described in the root README. |
| Table S2 | Candidate additional-copy review and read sources | `Supplementary_Data/Table_S4_CNV_assembly_artifact_qc.tsv`, `Table_S4_NucFreq_region_summary.tsv`, `Table_S4_raw_read_source_index.tsv`. |
| Table S3 | Short-read reconstruction and published-VCF comparisons | `Supplementary_Data/Table_S6/` supplies main Table 2. `Table_S6_short_read_vcf_locus_diagnostics.tsv` supplies Figure S1. |
| Table S4 | Solo-LTR diversity and expected differences | `Supplementary_Data/Table_S16/`. Supports Figure S18. |
| Table S5 | Structural-state source table | `Supplementary_Data/Table_S5_artifact_filtered_structural_spectrum.tsv`; `Figure_1_structural_observation_cells.tsv` provides the underlying observation summary. |
| Table S6 | Array copy-number distributions and calls | `Supplementary_Data/Table_S7_array_copy_number_distributions.tsv`, `Table_S7_haplotype_copy_number_calls.tsv`. |
| Table S7 | Tandem-array and solo-LTR allele counts | Table embedded in the manuscript supplement. Sources are `Supplementary_Data/Panel_data/tandem_resolved/`, `Table_S5_artifact_filtered_structural_spectrum.tsv` and `Table_S7_haplotype_copy_number_calls.tsv`. There is no separate file named Table_S3 in this release. |
| Table S8 | Nucleotide variation within arrays | `Supplementary_Data/Table_S15/`. Supports Figure S5, including conditional substitution-clock sensitivity and sequence-based closest-pair ties. |
| Table S9 | Phylogenetic and host-flank relationships | `Supplementary_Data/Table_S14/4q_acrocentric_relationship/`, `telomeric_phylogeny/`, `phylogeny_review/`. Supports Figure 3 and Figure S2. |
| Table S10 | Combined ORF-screen counts and 7p22.1 donor categories | `Supplementary_Data/Table_S8_artifact_filtered_orf_coding_potential.tsv` and `Table_S8_7p22_donor_categories.tsv`. |
| Table S11 | Type-I cassette comparisons | `Supplementary_Data/Table_S17/`. Supports Figure 5C, Figures S12–S13 and Supplementary Alignment 1. |
| Table S12 | Direct type-state audit | `Supplementary_Data/Table_S11_direct_TypeI_locus_calls.tsv` and `Table_S11/`. |
| Table S13 | Recurrent conversion and within-locus Type-I/II polymorphism | `Supplementary_Data/Table_S13/`, including frozen input, 54-scenario results, numerical checks and `analysis_code/`. |
| Table S14 | Functional association results and full testing families | `Supplementary_Data/Table_S10a_functional_237_model_results.tsv`, `Table_S10b_functional_51_refitted_models.tsv`, `Table_S10c_SLC44A5_disjoint_cohort_estimates.tsv`, `Table_S10d_anti_CD20_adjusted_state_models.tsv`, `Table_S10e_MAGE_complete_discovery_family.tsv.gz`, `Table_S10f_MAGE_candidate_eligibility_and_aliases.tsv`, `Table_S10g_GEUVADIS_complete_SLC44A5_followup_family.tsv`, `Table_S10h_MAGE_SLC44A5_HC3_sensitivity_models.tsv`. |
| Table S15 | Corrected 1q22 Gag state | `Supplementary_Data/Table_S9_1q22_Gag_artifact_corrected_truth.tsv`. |
| Table S16 | Archaic 8q11.23 junction reads and ape empty-site evidence | `Supplementary_Data/Table_S12/`. |
| Table S17 | GRCh38-disrupted annotations that pass in other alleles | Table embedded in the manuscript supplement. Source is the corrected record-level catalog, using the combined Intact/Intact_FS_End screen and 292 donors. There is no separate file named Table_S2 in this release. |
| Supplementary Alignment 1 | Complete cassette logos and aligned human/primate sequences | `Supplementary_Data/Table_S17/human_full_cassette_6000_7293.aligned.fa`, `primate_full_cassette_6000_7293.aligned.fa`, `human_1001bp_cassette_flanks.aligned.fa`, manifests and base-frequency table. `build_visuals.py` records the original 24-page PDF rendering. Follow the September revision README for release reproduction. |

## Reproduction labels

`scripts/reproduce_compact_panels.py` uses final labels Figure 7, S11, S14
and S15. Its additional `Type_I_lesion_context` plot is supporting context, not
a numbered manuscript figure. The September revision workflow covers the
revised Figure 3, Figure 5C, S5, S6, S7, S8, S13, S18 and Alignment 1. No ORF-break
insertion-dating analysis is included.

Four former displays (S3, S10, S11 and S12 in version 0.4.0) are omitted from the shortened supplement. Their data remain available. Telomeric assignment counts remain in `Figure_S5_assignment_resolution.tsv`, locus ORF frequencies in Table S10, and 7p22.1 donor categories in `Table_S8_7p22_donor_categories.tsv`. Historical data and script filenames are not final figure numbers.

Supplementary tables are numbered by first citation in the final manuscript. Historical source filenames remain unchanged. `SUPPLEMENTARY_TABLE_NUMBERING.tsv` gives each final label and its retained source label.
