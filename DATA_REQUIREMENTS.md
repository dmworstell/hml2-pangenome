# Analysis inputs

The code and derived-data archive is at https://doi.org/10.5281/zenodo.22759511.
Extract `HML2_derived_data_2026-09-14.zip` at the repository root to restore
the input paths. Its manifest lists each file's size and SHA-256 checksum.

| Analysis | Inputs |
|---|---|
| Structural and coding catalog | `project/results/resolved_manuscript_catalog_20260914/combined_hml2_orf_analysis.RESOLVED.tsv`, source observations, sample metadata and locus assignments |
| Phylogenies and nucleotide sharing | Extracted HML-2 sequences, selected sequence clusters, alignments, Newick trees and source witnesses |
| Type-I cassette and conversion calculations | Direct locus calls and the inputs and analysis code in `Supplementary_Data/Table_S13/` |
| Functional refits | Retained person-level marker and phenotype inputs, sampling-frame covariates and expression matrices from the sources in `DATA_SOURCES.md` |
| Genomic extraction and copy validation | Source whole-genome assemblies or graphs, long-read alignments, locus coordinates, copy-state controls and workflow manifests |
| Fiber-seq analyses | Collaborator-provided Fiber-seq inputs |

`external_inputs.tsv` records input identities for the principal figure builders.
Full sequence workflows additionally require MAFFT, minimap2, samtools, bedtools
or odgi as specified by each workflow. R scripts declare their package imports.
Python analysis dependencies are in `requirements-analysis.txt`.
