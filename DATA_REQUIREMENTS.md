# Analysis inputs

Use the derived-data asset matching the release version. Archived versions
are available at https://doi.org/10.5281/zenodo.22759510. Version 0.4.0 is
https://doi.org/10.5281/zenodo.22879194.
It includes the corrected Rec catalog and complete testing families.
Extract `HML2_derived_data_v0.4.0.zip` at the matching code repository root.
Do not overwrite a later checkout with an older archive. Run
`python scripts/restore_corrected_catalog.py` to extract both the current
catalog and the preceding catalog used only for replay.

| Analysis | Inputs |
|---|---|
| Structural and coding catalog | `project/results/rec_exon_boundary_correction_20260915/combined_hml2_orf_analysis.RESOLVED.REC_CORRECTED.tsv`, source observations, sample metadata and locus assignments |
| Rec boundary correction | The preceding short-ORF-corrected catalog, KCON reference, retained spliced alignment columns, paired call files and provenance receipts bundled under `project/results/rec_exon_boundary_correction_20260915/` |
| Phylogenies and nucleotide sharing | Extracted HML-2 sequences, selected sequence clusters, alignments, Newick trees and source witnesses |
| Type-I cassette and conversion calculations | Direct locus calls and the inputs and analysis code in `Supplementary_Data/Table_S13/` |
| Functional refits | Retained person-level marker and phenotype inputs, sampling-frame covariates and expression matrices from the sources in `DATA_SOURCES.md` |
| Genomic extraction and copy validation | Source whole-genome assemblies or graphs, long-read alignments, locus coordinates, copy-state controls and workflow manifests |
| Fiber-seq analyses | Collaborator-provided Fiber-seq inputs |

`external_inputs.tsv` records input identities for the principal figure builders.
Full sequence workflows additionally require MAFFT, minimap2, samtools, bedtools
or odgi as specified by each workflow. R scripts declare their package imports.
Python analysis dependencies are in `requirements-analysis.txt`.

The September 21 revision includes compact sequence inputs under
`analysis/september2026_revision/inputs/` and alignments in Tables S15–S17.
See that analysis README for the tested commands and original extraction inputs.
