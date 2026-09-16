# Structural polymorphism and population-variable coding capacity of HERV-K(HML-2) in human pangenomes

Source code and result tables for the HML-2 pangenome manuscript. The current
GitHub revision includes the manuscript-review corrections of 15 September 2026.
The archived release 0.2.0 remains the earlier 14 September snapshot.

The analysis covers structural variation, tandem-copy validation, coding
potential, phylogenetic relationships, Type-I cassette variation, target-site
duplications, short-read comparisons, and exploratory functional associations.

## Reproduce the compact results

With Python 3.12, install `requirements.txt`, then run:

```sh
python scripts/check_release.py
python scripts/check_review_corrections.py
python scripts/reproduce_compact_panels.py --output outputs/retained_panels
python project/manuscript/restore_figure3_20260914.py
python Supplementary_Data/Table_S13/analysis_code/type1_observed_loci_uniformity_20260914.py
```

The release check verifies the bundled tables, chromosome denominators,
structural-state tests, recovered NucFreq method, retained model formulas,
deterministic model evidence, and Python syntax. The review-correction check
reconstructs the reported multiple-testing families and checks the corrected
short-ORF and missing-haplotype records. The compact-panel command
regenerates Figure 7, Supplementary Figures S6 and S13, and supporting model
plots. The Figure 3 command reproduces the current figure, including all-locus
trees, exact nucleotide sharing and the chromosome-4 connection. Use this
command for Figure 3, not the historical general-purpose compositor.
The final command reproduces the 54 recurrent-conversion and drift scenarios
in Table S13. These commands do not download data or run a cluster job.

## Code map

| Location | Analysis |
|---|---|
| `project/manuscript/` | Current figure builders, artifact filtering, functional refits, TSD analysis and conditional mechanism models |
| `project/working/biological_orf_annotation_20260802/` | Host-flank-supported locus assignment and biological ORF annotation |
| `project/working/cnv_copy_state_reinterpretation_v1/` | Copy-state calibration and depth evidence |
| `project/cluster_workflows/cnv_*/` | Sources for the completed targeted ONT copy-validation workflows |
| `project/working/direct_*` and `onep31b_slc44a5_followup_v1/` | EBV/LCL outcomes and expression follow-up |
| `project/working/hml2_functional_evidence_synthesis_claude_v2/` | Donor-disjoint SLC44A5 effect synthesis |
| `manuscript_figures/python/analysis/` | Type-I cassette/backbone comparison and matched long-read/short-read comparison |
| `manuscript_figures/R/` | Source for mutation, duplication, fusion-ORF and Fiber-seq panels |
| `HML2_ProjectResources/cluster_pipeline_source/shared/` | Upstream ORF and flank-anchoring source |
| `Supplementary_Data/` | Current manuscript source tables, alignments, trees, Table S13 analysis code and checksums |
| `data/` | Additional compact model and copy-validation evidence |

The current plot builders use 292 donor IDs.

## Full analysis inputs

The [archived data](https://doi.org/10.5281/zenodo.22759511),
`HML2_derived_data_2026-09-14.zip`, contain the 14 September catalog,
extracted HML-2 sequences, supplementary tables and retained analysis inputs.
That immutable archive predates the current review corrections. Its
`file_manifest.tsv` records the archived file identities.

For full-input reproduction, extract the archive into a separate directory,
then copy only files absent from this checkout. Do not overwrite the newer
repository tables. Run `python scripts/restore_corrected_catalog.py` to extract
the corrected catalog bundled in `data/`. The current figure and functional
refit builders read that corrected catalog. Complete discovery/follow-up
testing families, including all 1,289,856 MAGE discovery tests, are in
`Supplementary_Data/Table_S10*`.

Three historical source-provenance fields in the distributed catalog use
neutral `historical_source/` prefixes instead of local computer paths. These
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

`source_inventory.tsv` records the source-file identity used for this snapshot.
`CHANGES.md` lists the publication-specific corrections. The supplementary
manifest records both the original source checksum and the distributed checksum
where line endings or local-path metadata were normalized. Numeric values were
not changed by that normalization.

Original code is MIT-licensed and original derived data are CC BY 4.0.
Third-party sources retain their existing terms.
