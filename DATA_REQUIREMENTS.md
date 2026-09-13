# Additional data required for full reproduction

The compact checks and panels described in the README run from this repository.
The complete assembly-to-manuscript analysis needs the additional inputs below.
The code and compact table release is not a claim that these inputs have been
deposited or that a full genomic rerun has been completed.

| Analysis | Required inputs not bundled |
|---|---|
| Catalog and complete coding/structural figures | The 176 MB CNV-weighted, biologically annotated v3 ORF table, its summary, the v3r1 precursor, sample frame, locus/subfamily authority and exact artifact/copy-state decisions |
| Genomic extraction and ORF scoring | Source assembly/graph FASTAs, locus coordinates and sample universe, Type-I/Type-II reference sequences, alignment tools and frozen extraction manifests |
| Targeted ONT copy validation | The retained ONT alignments or source reads, assembly regions, expected-copy controls, CNV input manifests and calibration receipts |
| 8q11.23 network re-derivation | The original provirus and extracted solo-LTR FASTAs, plus sample superpopulation metadata. The compact table can be replotted without them |
| Phylogenies and mutation maps | Processed locus sequences, KCON-anchored alignments, Newick trees, locus/subfamily assignments and mutation tables |
| Type-I cassette analysis | The full ORF master, delta-292 summary/record files, processed and graph-locus FASTAs, and Type-II KCON reference |
| Functional refits | Person-level direct marker/phenotype matrices, sampling-frame/pedigree covariates, MAGE and GEUVADIS expression data, anti-CD20 outcomes and original source-result files |
| Fiber-seq and fusion-ORF panels | Collaborator-authorized Fiber-seq inputs and the original fusion-ORF inputs. They are not redistributed here |
| Original artwork/panel assembly | The original author-owned or licensed graphical assets, which are not included in this code repository |

`external_inputs.tsv` records exact local input identities for the principal
retained figure builders where available. Place these under the same
repository-relative layout before using the full builders. Input files not
listed there may still be requested by downstream scripts. Those scripts fail
on missing inputs rather than fabricate a result.

Public source studies retain their own data access and redistribution terms.
An archival accession/DOI for the complete reproducibility bundle is not yet
assigned. A license for the newly released source code also remains for the
authors to select.
