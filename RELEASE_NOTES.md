# HML-2 pangenome v0.4.3 preparation

This repository revision contains the corrected public code, data and Figures 2–7 for the 22 September 2026 manuscript revision.

- Replace the consensus-derived whole-ORF recovery calculation with comparison of ORF-associated small variants to native high-coverage 1000 Genomes genotypes. The panel has 20,500 carrier observations: 13,717 supported alternate, 118 supported reference, 10 other alternate and 6,655 unresolved. Missing or inadequate genotypes remain unresolved; these tests do not establish complete ORF sequence or phase.
- Include all 583 retained 8q11.23 solo-LTR haplotypes, their source sequences and the reference-alignment evidence. Fourteen sequences are observed; the dominant sequence has 547 observations.
- Correct the Type I deletion-class model's conservative Delta292 count from 10 to 16; the inclusive count remains 17. Retain prior model code as historical analysis.
- Include the reviewed public figure assets, revised Figure S13 label recipe and relative-viability label for Figure 7C.
- Preserve immutable biological inputs. New text provenance uses neutral historical_source prefixes and retains original/distributed checksums; frozen native VCF bytes retain their retrieval headers.

The prior v0.4.2 release remains available at https://doi.org/10.5281/zenodo.22903723. Its consensus-derived Table 2 calculation is superseded by this revision. All other retained public inputs remain in the matching full derived-data archive.

The main manuscript, private collaborator Fiber-seq inputs, cover letters, reviewer/contact documents and licensed BioRender artwork are excluded. No ONT or Dating production pipeline was run for this release.

Archive status: prepared locally; no new Zenodo draft or DOI exists. The public v0.4.2 archive remains unchanged.
