# HERV-K (HML-2) variation in human pangenomes

Source code and compact result tables for the HML-2 pangenome manuscript.

The analysis covers structural variation, tandem-copy validation, coding
potential, phylogenetic relationships, Type-I cassette variation, target-site
duplications, short-read comparisons, and exploratory functional associations.

## Reproduce the compact results

With Python 3.12, install `requirements.txt`, then run:

```sh
python scripts/check_release.py
python scripts/reproduce_compact_panels.py --output outputs/retained_panels
```

The first command verifies the bundled tables, chromosome denominators,
structural-state tests, recovered NucFreq method, retained model formulas,
deterministic model evidence, and Python syntax. The second
regenerates Figure 7, Supplementary Figures S7, S13, S14, S17 and S20, and the
duplicated-group schematic from the bundled compact data. It does not download
data or run a cluster job.

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
| `Supplementary_Data/` | The 18 manuscript source tables and checksum manifest |
| `data/` | Additional compact model and copy-validation evidence |

The current plot builders use 292 donor IDs. Reference assemblies are not
additional donors. A sequence-compatible coding call can include a terminal
frameshift where the stated analysis includes that category. It is not evidence
of protein expression or biological activity.

## Full analysis inputs

This repository is not a self-contained archive of the genomic inputs.
`DATA_REQUIREMENTS.md` lists the additional inputs needed for the full analysis.
The large ORF tables, assembly/graph FASTAs, sequencing alignments and raw
collaborator data are not bundled. Figure source code is included even where
those inputs are still needed.

Additional Python packages for full-input analyses are in
`requirements-analysis.txt`. The R sources declare their own package imports.
External executables such as MAFFT, minimap2, samtools, bedtools and odgi are
required by the relevant sequence workflows, not by the compact reproduction.

The targeted CNV workflow source retains its original manifest-driven design.
Paths shown as `/path/to/...` are installation examples, not an available
cluster. These workflows require their frozen manifests, sequencing inputs and
an explicitly configured compute allocation. They have not been rerun as part
of this source release.

The Type-I source-versus-effect and helper analyses are conditional sensitivity
models. They do not estimate a measured fitness effect. Their exact assumptions
and reproduction checks are in `MODEL_NOTES.md`.

## Provenance

`source_inventory.tsv` records the source-file identity used for this snapshot.
`CHANGES.md` lists the publication-specific corrections. The supplementary
manifest records both the original source checksum and the distributed checksum
where line endings or local-path metadata were normalized. Numeric values were
not changed by that normalization.

No third-party license has been added to collaborator data or figure artwork.
This repository does not include such artwork or assert redistribution rights
over the original external datasets.
