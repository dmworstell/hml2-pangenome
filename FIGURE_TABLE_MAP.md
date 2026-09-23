# Manuscript figure and table map

This map follows the manuscript and supplement revised on 22 September 2026:
seven main figures and Figures S1–S17. It describes the corrected code/data revision prepared for v0.4.3.
The successor archive is not published and has no reserved DOI. The published
v0.4.2 DOI archive remains an immutable historical snapshot.

Paths are relative to the repository root. Retained source filenames often use
older numbering. A source or plotting function listed here identifies provenance,
not proof that running it recreates the current figure. In particular, older
whole-ORF recovery, assembly-only 8q11.23, lower-count Type-I and helper-population
analyses are superseded. Current numerical sources take precedence.

Public figure exports are in `Figures/`, including the current `Figure_7.png`
and `Figure_7.pdf`. The main manuscript, BioRender artwork and private
collaborator Fiber-seq data are not redistributed. The table below identifies
source data even when a complete public figure export is unavailable.

## Main figures

| Final label | Content | Source data and analysis |
|---|---|---|
| Figure 1 | Structural variation | Panel A is licensed BioRender artwork. Panel B uses `Supplementary_Data/Figure_1_structural_source_data.tsv`, `Figure_1_structural_observation_cells.tsv` and `Figure_1_sex_chromosome_eligibility.tsv`; `project/manuscript/build_narrative_main_figures.py`, `build_figure_1`. |
| Figure 2 | Locus copy counts and within-array coding variation | `Supplementary_Data/Table_S7_array_copy_number_distributions.tsv`, `Table_S7_haplotype_copy_number_calls.tsv`, `Figure_2_within_array_variation.tsv`; `project/manuscript/build_narrative_main_figures.py`, `build_figure_2`. Nucleotide differences are separately shown in Figure S5. |
| Figure 3 | Regional phylogenies, sequence sharing and chromosome locations | Current results and rendering source are in `Supplementary_Data/Table_S14/phylogeny_review/`; host-flank evidence is in `Supplementary_Data/Table_S14/4q_acrocentric_relationship/`. The current panel B reports mean pairwise nucleotide differences, not counts of shared sequences. |
| Figure 4 | ORF annotations and person-level coding potential | `Supplementary_Data/Table_S8_artifact_filtered_orf_coding_potential.tsv`, `Figure_4_person_orf_burden.tsv`, `Figure_4_locus_orf_superpopulation_frequencies.tsv`; `project/manuscript/build_narrative_main_figures.py`, `build_figure_4`. |
| Figure 5 | Type-I cassette conservation and deletion-class propagation models | Panel C uses `Supplementary_Data/Table_S17/window_summary.tsv` and `window_pairwise_differences.tsv`, including the LTR5Hs-only Type-II control. Current panel D uses `Supplementary_Data/TypeI_deletion_models/`: 16 resolved Δ292-bearing insertions plus one possible additional insertion, compared with six other deletion classes. The current conservative and inclusive count vectors are `[16,1,1,1,1,1,1]` and `[17,1,1,1,1,1,1]`. The earlier lower-count source-effect table is superseded. |
| Figure 6 | 8q11.23 provirus structure and solo-LTR haplotypes | `data/8q11_23_provirus_annotation.json`, `Supplementary_Data/Figure_7_solo_LTR_haplotype_counts.tsv` and `Supplementary_Data/Table_S16/`. The corrected sequence set contains all 583 solo-LTR donor-haplotypes, 14 sequence haplotypes and 547 observations of the dominant sequence. The `Figure_7` source filename is historical. The compact runner labels its network output `Figure_6`. |
| Figure 7 | Structural states and cellular phenotypes | `Supplementary_Data/Table_S10a_functional_237_model_results.tsv`, `Table_S10c_SLC44A5_disjoint_cohort_estimates.tsv`, `Table_S10d_anti_CD20_adjusted_state_models.tsv`, and `Functional_refit/`. Current public exports are `Figures/Figure_7.png` and `Figures/Figure_7.pdf`. Panel C reports adjusted differences in relative viability for 12q13.2 provirus carriers versus solo-LTR carriers without a provirus. Historical source function `build_figure_6` in `project/manuscript/build_narrative_main_figures.py` predates the current numbering and labeling. |

## Main tables

| Final label | Content | Source |
|---|---|---|
| Table 1 | Cohorts and denominators | Retained donor roster and the listed comparison cohorts. |
| Table 2 | Short-read genotype support for ORF-associated variants observed in long-read assemblies | `Supplementary_Data/Table_S3/ORF_variant_comparison/Data/variant_comparison.summary.tsv` and `variant_comparison.rows.tsv`. The selected panel contains 564 normalized alleles at 36 loci and 20,500 donor–locus–variant observations. Native unphased genotypes are classified as supported alternate, supported reference, other alternate or unresolved. These calls do not establish whole-ORF sequence, phase or all-base callability. |

## Supplementary figures

| Final label | Content | Source data and analysis |
|---|---|---|
| Figure S1 | Illumina ensemble structural-variant calls compared with matched long-read assemblies | `Supplementary_Data/Table_S6_short_read_vcf_locus_diagnostics.tsv`. The plotted calls come from the Illumina ensemble callset alone. Main Table 2 instead uses native short-read genotypes at selected ORF-associated variants. Historical plotting function `build_short_read_figure` in `project/manuscript/build_artifact_filtered_catalog_figures.py` writes a filename beginning `Figure_2`. |
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
| Figure S12 | Regional nearest-neighbor changes, linked cassette sites and deletion-class insertion counts | `manuscript_figures/python/analysis/type1_cassette_mechanism/`, `project/manuscript/delta292_ape_counterfactuals_v1/` and `Supplementary_Data/Table_S17/` retain sequence-comparison sources. Panel C uses the orthology-collapsed inputs in `Supplementary_Data/TypeI_deletion_models/inputs/`: 16 resolved Δ292-bearing insertions and one possible additional insertion. The 83 primate clusters are separate sequence observations before orthology collapse. This figure uses zero-based KCON positions 6331/6492 and deletion interval [6501,6793). |
| Figure S13 | Human cassette frequency logos and similarity ranking | `Supplementary_Data/Table_S17/human_cassette_site_frequencies.tsv`, `candidate_similarity_to_TypeI.tsv`, `human_representatives.tsv`; September revision workflow. This figure uses one-based KCON sites 6332/6493 and deletion interval 6502–6793. |
| Figure S14 | Deletion-class contributions and Δ292 enrichment | `Supplementary_Data/TypeI_deletion_models/`, including the current orthology-collapsed inputs and `reproduce_source_model.py`. Panel A compares conservative 16-plus-six and inclusive 17-plus-six counts under alternative detection odds. Panel B gives the median relative contribution of the Δ292 class without an added deletion-specific effect. The primary conservative-count, narrow-odds curve appears in Figure 5D. Historical mechanism-ranking and helper-population outputs do not supply this figure. |
| Figure S15 | Solo-LTR diversity and illustrative substitution expectations | `Supplementary_Data/Table_S16/`, including `solo_LTR_diversity_comparison.tsv` and `verify_8q11_23_diversity.py`. The corrected 8q11.23 alignment contains 583 sequences. The other twelve loci retain their validated sequence sets. This comparison does not estimate an insertion or coalescence date. The former assembly-only renderer and `Figure_S21` source names are historical. |
| Figure S16 | Gag dosage, cell growth and EBV abundance | `Supplementary_Data/Table_S9_1q22_Gag_artifact_corrected_truth.tsv`, `Table_S10a_functional_237_model_results.tsv`, `Table_S10b_functional_51_refitted_models.tsv`, and `Functional_refit/`; `project/manuscript/build_current_results_figure.py`. |
| Figure S17 | Fiber-seq chromatin accessibility | Retained analysis source `manuscript_figures/R/analysis/Fiberseq_haplotype_reanalysis.R` and public `Supplementary_Data/Figure_S17_locus_crosswalk.tsv`. Private collaborator source data and figure artwork are not included in this public collection. |

## Supplementary tables and alignment data

| Final label | Content | Public source |
|---|---|---|
| Table S1 | Structural-analysis catalog | `Supplementary_Data/Table_S1_structural_observations_corrected.tsv`; complete corrected catalog is restored from `data/` as described in the root README. |
| Table S2 | Candidate additional-copy review and read sources | `Supplementary_Data/Table_S4_CNV_assembly_artifact_qc.tsv`, `Table_S4_NucFreq_region_summary.tsv`, `Table_S4_raw_read_source_index.tsv`. |
| Table S3 | ORF-associated variant comparison and structural-variant diagnostics | `Supplementary_Data/Table_S3/ORF_variant_comparison/` supplies native-genotype comparisons for main Table 2, call-level quality evidence and intact-copy tested/untested inventories. The separate `Table_S6_short_read_vcf_locus_diagnostics.tsv` supplies Figure S1. The retired `Table_S6/` consensus-comparison subtree is not the current source. |
| Table S4 | Solo-LTR diversity and expected differences | `Supplementary_Data/Table_S16/`. Supports Figure 6B and Figure S15 after correction to all 583 8q11.23 solo-LTR donor-haplotypes. |
| Table S5 | Structural-state source table | `Supplementary_Data/Table_S5_artifact_filtered_structural_spectrum.tsv`; `Figure_1_structural_observation_cells.tsv` provides the underlying observation summary. |
| Table S6 | Array copy-number distributions and calls | `Supplementary_Data/Table_S7_array_copy_number_distributions.tsv`, `Table_S7_haplotype_copy_number_calls.tsv`. |
| Table S7 | Tandem-array and solo-LTR allele counts | Table embedded in the manuscript supplement. Sources are `Supplementary_Data/Panel_data/tandem_resolved/`, `Table_S5_artifact_filtered_structural_spectrum.tsv` and `Table_S7_haplotype_copy_number_calls.tsv`. Its historical source label is Table S3; the current `Table_S3/` directory instead contains the native variant comparison for final Table S3. |
| Table S8 | Nucleotide variation within arrays | `Supplementary_Data/Table_S15/`. Supports Figure S5, including conditional substitution-clock sensitivity and sequence-based closest-pair ties. |
| Table S9 | Phylogenetic and host-flank relationships | `Supplementary_Data/Table_S14/4q_acrocentric_relationship/`, `telomeric_phylogeny/`, `phylogeny_review/`. Supports Figure 3 and Figure S2. |
| Table S10 | Combined ORF-screen counts and 7p22.1 donor categories | `Supplementary_Data/Table_S8_artifact_filtered_orf_coding_potential.tsv` and `Table_S8_7p22_donor_categories.tsv`. |
| Table S11 | Type-I cassette comparisons and deletion-class models | `Supplementary_Data/Table_S17/` supplies complete human/primate aligned FASTA files, sequence manifests and nucleotide-frequency tables for Figure 5C and Figures S12–S13. `Supplementary_Data/TypeI_deletion_models/` supplies the corrected orthology inputs and model results for Figure 5D and Figure S14. |
| Table S12 | Direct type-state audit | `Supplementary_Data/Table_S11_direct_TypeI_locus_calls.tsv` and `Table_S11/`. |
| Table S13 | Recurrent conversion and within-locus Type-I/II polymorphism | `Supplementary_Data/Table_S13/`, including frozen input, 54-scenario results, numerical checks and `analysis_code/`. |
| Table S14 | Archaic 8q11.23 junction reads and ape empty-site evidence | `Supplementary_Data/Table_S12/`. |
| Table S15 | Functional association results and full testing families | `Supplementary_Data/Table_S10a_functional_237_model_results.tsv`, `Table_S10b_functional_51_refitted_models.tsv`, `Table_S10c_SLC44A5_disjoint_cohort_estimates.tsv`, `Table_S10d_anti_CD20_adjusted_state_models.tsv`, `Table_S10e_MAGE_complete_discovery_family.tsv.gz`, `Table_S10f_MAGE_candidate_eligibility_and_aliases.tsv`, `Table_S10g_GEUVADIS_complete_SLC44A5_followup_family.tsv`, `Table_S10h_MAGE_SLC44A5_HC3_sensitivity_models.tsv`. |
| Table S16 | Corrected 1q22 Gag state | `Supplementary_Data/Table_S9_1q22_Gag_artifact_corrected_truth.tsv`. |
| Table S17 | GRCh38-disrupted annotations that pass in other alleles | Table embedded in the manuscript supplement. Source is the corrected record-level catalog, using the combined Intact/Intact_FS_End screen and 292 donors. There is no separate file named Table_S2 in this release. |
| Table S11 alignment data | Complete cassette logos and aligned human/primate sequences | `Supplementary_Data/Table_S17/human_full_cassette_6000_7293.aligned.fa`, `primate_full_cassette_6000_7293.aligned.fa`, `human_1001bp_cassette_flanks.aligned.fa`, manifests and base-frequency table. The separate 24-page alignment PDF is omitted from the submission. `build_visuals.py` retains its historical rendering, and all sequence and frequency data remain available. Follow the September revision README for release reproduction. |

## Reproduction labels and limits

The current supplement contains captions S1–S17. The former helper-dependence
figure is absent, and there is no current Figure S18 assignment. Historical
references to eighteen supplementary figures predate this revision.

`scripts/reproduce_compact_panels.py` labels the 8q11.23 network `Figure_6` and
retains the direct cassette-call panel as `Figure_S11`. Its historical mechanism
and helper outputs are labeled `Historical_mechanism_ranking` and
`Historical_helper_equilibrium`; they are not Figures S14 or S15. Supporting
context plots are not numbered manuscript figures.

Use `Supplementary_Data/TypeI_deletion_models/reproduce_source_model.py` for
current Figure 5D/S14 numerical results, and the verifier in
`Supplementary_Data/Table_S16/` for the corrected 8q11.23 data. The earlier
September workflow retains useful phylogeny, array, cassette and
terminal-junction analyses, but its historical numbering and older solo-LTR
inputs do not reproduce every current panel. Retained source renderers are
not a verified end-to-end assembly command for the current manuscript.

`Supplementary_Data/SUPPLEMENTARY_TABLE_NUMBERING.tsv` gives each current
supplementary table label and its retained source label.
`Supplementary_Data/FIGURE_SOURCE_NUMBERING.tsv` identifies the changed
8q11.23, solo-LTR-diversity and Fiber-seq source mappings. Biological identities
and sequence coordinates are unchanged by renumbering. No ORF-break
insertion-dating analysis is included here.

## Current Table 2 and Table S3

Use `Supplementary_Data/Table_S3/ORF_variant_comparison/README.md` for frozen-input
reproduction. Missing records and inadequate genotypes remain unresolved.
The previous consensus-derived whole-ORF recovery analysis is superseded;
its recovery percentages and zero multi-copy/21% solo-LTR statements are not
current results. The selected variant panel cannot establish a complete ORF
or phase across one provirus.
