# Structural polymorphism and population-variable coding capacity of HERV-K(HML-2) in human pangenomes

Source code and result tables for the HML-2 pangenome manuscript. Release 0.2.0
is frozen from the corrected manuscript inputs on 14 September 2026.

The analysis covers structural variation, tandem-copy validation, coding
potential, phylogenetic relationships, Type-I cassette variation, target-site
duplications, short-read comparisons, and exploratory functional associations.

## Reproduce the compact results

With Python 3.12, install `requirements.txt`, then run:

```sh
python scripts/check_release.py
python scripts/reproduce_compact_panels.py --output outputs/retained_panels
python project/manuscript/restore_figure3_20260914.py
python Supplementary_Data/Table_S13/analysis_code/type1_observed_loci_uniformity_20260914.py
```

The first command verifies the bundled tables, chromosome denominators,
structural-state tests, recovered NucFreq method, retained model formulas,
deterministic model evidence, and Python syntax. The second
regenerates Figure 7, Supplementary Figures S6 and S13, and supporting model
plots. The third command reproduces the current Figure 3, including all-locus
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
`HML2_derived_data_2026-09-14.zip`, contain the corrected catalog,
extracted HML-2 sequences, current supplementary tables and retained analysis
inputs. Extract it at the repository root to restore the relative paths.
Its `file_manifest.tsv` records the archived file identities.

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
