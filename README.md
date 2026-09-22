# Structural polymorphism and population-variable coding capacity of HERV-K(HML-2) in human pangenomes

Source code and result tables for the HML-2 pangenome manuscript. The current
release 0.4.1 includes the manuscript revision of 22 September 2026. It adds
array nucleotide comparisons, Type I cassette controls, solo-LTR diversity,
phylogenetic bootstrap support and corrected terminal-junction evidence.

The analysis covers structural variation, tandem-copy validation, coding
potential, phylogenetic relationships, Type-I cassette variation, target-site
duplications, short-read comparisons, and exploratory functional associations.

## Reproduce the compact results

With Python 3.12, install `requirements.txt`, then run:

```sh
python scripts/check_release.py
python scripts/check_review_corrections.py
python scripts/reproduce_compact_panels.py --output outputs/retained_panels
python Supplementary_Data/Table_S13/analysis_code/type1_observed_loci_uniformity_20260914.py
```

The release check verifies the bundled tables, chromosome denominators,
structural-state tests, recovered NucFreq method, retained model formulas,
deterministic model evidence, and Python syntax. The review-correction check
reconstructs the reported multiple-testing families and checks the corrected
short-ORF and missing-haplotype records. The compact-panel command
regenerates retained panels for Figure 7 and Supplementary Figures S11, S14 and S15, plus supporting model plots. The final command reproduces the 54 recurrent-conversion and drift scenarios
in Table S13. These commands do not download data or run a cluster job.

For the September 21 analyses and revised Figure 3, follow
[`analysis/september2026_revision/README.md`](analysis/september2026_revision/README.md).
That workflow uses retained sequence inputs and reproduces the revised tables.
The older Figure 4 and TSD builders preserve their historical analyses.
The current TSD measurements are in `Supplementary_Data/Figure_S7_boundary_evidence/`.

See [FIGURE_TABLE_MAP.md](FIGURE_TABLE_MAP.md) for final manuscript numbering
and the relationship to historical source filenames.

## Code map

| Location | Analysis |
|---|---|
| `project/manuscript/` | Retained figure builders, artifact filtering, functional refits and mechanism models. Revised Figure 3, Figure 5C and Figures S5, S6, S13 and S18 use the September workflow below. |
| `project/working/biological_orf_annotation_20260802/` | Host-flank-supported locus assignment and biological ORF annotation |
| `project/working/cnv_copy_state_reinterpretation_v1/` | Copy-state calibration and depth evidence |
| `project/cluster_workflows/cnv_*/` | Sources for the completed targeted ONT copy-validation workflows |
| `project/working/direct_*` and `onep31b_slc44a5_followup_v1/` | EBV/LCL outcomes and expression follow-up |
| `project/working/hml2_functional_evidence_synthesis_claude_v2/` | Donor-disjoint SLC44A5 effect synthesis |
| `manuscript_figures/python/analysis/` | Type-I cassette/backbone comparison and matched long-read/short-read comparison |
| `manuscript_figures/R/` | Source for mutation, duplication, fusion-ORF and Fiber-seq panels |
| `HML2_ProjectResources/cluster_pipeline_source/shared/` | Upstream ORF and flank-anchoring source |
| `Supplementary_Data/` | Current manuscript source tables, alignments, trees, analysis code and checksums |
| `analysis/september2026_revision/` | Reproduction of the revised array, solo-LTR, Type I, phylogeny and terminal-junction analyses |
| `data/` | Additional compact model and copy-validation evidence |

The current plot builders use 292 donor IDs.

## Full analysis inputs

The [v0.4.1 archive](https://doi.org/10.5281/zenodo.22888353),
`HML2_derived_data_v0.4.1.zip`, contains the corrected catalog, extracted HML-2
sequences, supplementary tables and retained analysis inputs. Its
`file_manifest.tsv` records every archived file identity.

For full-input reproduction, extract the matching code and derived-data ZIPs
into the same directory, preserving their relative paths. If using a later
GitHub revision, do not overwrite its files with an older archive.
Run `python scripts/restore_corrected_catalog.py` to extract the two catalogs
bundled in `data/` for the Rec replay. The current figure and functional
refit builders read that corrected catalog. Complete discovery/follow-up
testing families, including all 1,289,856 MAGE discovery tests, are in
`Supplementary_Data/Table_S10*`.

The current catalog is `REC_CORRECTED`, which includes the complete Rec stop
codon. The older `SHORT_ORF_CORRECTED` catalog is retained only as the input to
the documented correction. After extracting both catalogs, run
`python scripts/check_rec_boundary.py` to reproduce all 36,073 Rec-bearing
records from the retained alignment slices and compare the resulting catalog
byte-for-byte with the current distribution. This replay also checks paired
call files for three alignments regenerated after the older catalog snapshot.

Three historical source-provenance fields in the distributed catalog use
neutral `historical_source/` prefixes instead of local or cluster paths. These
labels identify provenance, not additional bundled files. The catalog
verification record retains both the original and distributed checksums.
Biological calls and numerical fields are unchanged by this normalization.
Historical paths in supplementary table provenance are likewise shortened to
relative or neutral historical-source labels. The supplementary manifest records
the original package checksum and the distributed checksum for each file.

See `DATA_REQUIREMENTS.md` for workflow inputs and `DATA_SOURCES.md` for
the source datasets and accessions.

Additional Python packages for full-input analyses are in
`requirements-analysis.txt`. The R sources declare their own package imports.
External executables such as MAFFT, minimap2, samtools, bedtools and odgi are
required by the relevant sequence workflows, not by the compact reproduction.

The targeted CNV workflows use input manifests and a configured compute
allocation. Model assumptions and reproduction checks are in `MODEL_NOTES.md`.

## Provenance

The figure builders label the short-arm HML-2 groups as "Telomeric Type I"
and "Telomeric Type II". Historical filenames may retain earlier figure numbers.
The September revision README maps each analysis to the final manuscript figures.
Previous archived versions remain available through
[Zenodo](https://doi.org/10.5281/zenodo.22759510).

`source_inventory.tsv` records the source-file identity used for this snapshot.
`CHANGES.md` lists the publication-specific corrections. The supplementary
manifest records both the original source checksum and the distributed checksum
where line endings or local-path metadata were normalized. Numeric values were
not changed by that normalization.

Original code is MIT-licensed and original derived data are CC BY 4.0.
Third-party sources retain their existing terms.

## Matched VCF-derived and long-read calls

Main Table 2 and Table S3 compare annotations derived from the phased high-coverage 1000 Genomes VCF panel with long-read calls in 282 shared donors at 83 loci after artifact filtering. The phased variants were applied to hg38 and the resulting haplotype sequences were annotated for coding and structural states. The compact catalog inputs, complete results and reproduction command are in [Table S3](Supplementary_Data/Table_S6/README.md). The comparison retains missing data as unknown and requires linked coding genes on the same provirus.

Figure S1 compares Illumina ensemble structural-variant calls with matched long-read assemblies. It uses that callset alone. Four VCF resources were examined in the broader work, including the phased panel used for the coding comparison. They are not four additional inputs to Table 2 or a pooled source for Figure S1.

Four redundant supplemental displays have been removed from the manuscript. Their data remain in this archive. See the figure map for the current seven main and eighteen supplemental figures.

### Retained supplementary panels S4 and S10

`Rscript scripts/reproduce_s4.R outputs/Figure_S4` uses the two retained tandem-array summary tables and the current `Duplication_Analysis.R` plotting blocks. It requires R with ggplot2, scales, dplyr, readr and stringr. Optional ragg and systemfonts packages control the graphics device and font selection.

`python scripts/reproduce_s10.py --output outputs/Figure_S10` uses the root Python environment. It checks the annotation-tier and longest-frame counts against the retained record-level witness tables before rendering the current S10 palette. This is a replot of the retained result, not a new sequence scan or a reclassification with another catalog. The historical source directory remains `Supplementary_Data/Figure_S9/`.
