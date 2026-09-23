# Structural polymorphism and population-variable coding capacity of HERV-K(HML-2) in human pangenomes

Source code, curated result tables and public figure exports for the HML-2
pangenome manuscript. This code/data revision, prepared for v0.4.3, contains the corrected
native variant comparison, complete 8q11.23 solo-LTR sequence set and revised
Type-I model inputs. The successor archive is not published and has no reserved DOI. The
[published v0.4.2 archive](https://doi.org/10.5281/zenodo.22903723) remains an
immutable historical snapshot and does not contain these corrections.

The analyses cover structural variation, tandem-copy validation, coding
potential, regional phylogenies, Type-I cassette variation, target-site
duplications, short-read comparisons and exploratory functional associations.
The retained HML-2 catalog includes 59,656 analysis records from 292 donors.

## Current results and inputs

- Main Table 2 and the ORF-associated variant part of Table S3 use
  [`Supplementary_Data/Table_S3/ORF_variant_comparison/`](Supplementary_Data/Table_S3/ORF_variant_comparison/README.md).
  Native short-read genotypes replace the superseded consensus-derived
  coding/structural recovery calculation. The current comparison tests selected
  variants; it does not reconstruct complete ORFs.
- The corrected 8q11.23 analysis includes all 583 solo-LTR donor-haplotypes:
  461 assembly-derived and 122 graph-derived sequences. There are 14 sequence
  haplotypes, including 547 observations of the dominant sequence. The data and
  verification script are in
  [`Supplementary_Data/Table_S16/`](Supplementary_Data/Table_S16/README.txt),
  which supplies current Table S4, Figure 6B and Figure S15.
- The current Type-I deletion-class model uses 16 resolved ancestral
  Δ292-bearing insertions and one possible additional insertion. Its
  conservative count vector is `[16,1,1,1,1,1,1]`; the inclusive vector is
  `[17,1,1,1,1,1,1]`. Inputs and reproduction code are in
  [`Supplementary_Data/TypeI_deletion_models/`](Supplementary_Data/TypeI_deletion_models/README.md).
  These results replace the earlier lower-count analysis. The former
  helper-population and mechanism-ranking results are not current manuscript
  results.
- Main Figure 6 describes 8q11.23. Main Figure 7 describes cellular phenotypes;
  panel C reports adjusted differences in relative viability. Public rendered
  figures are stored in `Figures/`, including `Figure_7.png` and `Figure_7.pdf`.

See [FIGURE_TABLE_MAP.md](FIGURE_TABLE_MAP.md) for the seven main figures,
Figures S1–S17 and their historical source filenames. The public collection
excludes the main manuscript, BioRender artwork and private collaborator
Fiber-seq source data.

## Native ORF-associated variant comparison

The shared catalog panel contains 282 donors and 83 loci. Native regional
genotypes were recovered for the 79 loci that have mapped hg38 element
intervals. The selected target set contains 564 normalized alleles at 36 loci
and 20,500 donor–locus–variant observations. It includes single-nucleotide
premature-stop gains/losses in canonical reading frames and indels shorter than
50 bp whose length changes are not divisible by three. Programmed Gag/Pro/Pol
frame transitions are not sequence defects. Type-I Env denotes a theoretical
N-terminally truncated product whose translation has not been demonstrated.

Passing calls require complete GT, DP ≥10, GQ ≥20 and passing site/sample
filters. The results contain 13,717 supported alternate calls, 118 supported
reference calls, 10 other alternate calls and 6,655 unresolved observations.
Concordance is 99.0755% among 13,845 callable pairs; 32.4634% of the selected
panel remains unresolved. Absent records and inadequate calls are unresolved,
not reference calls. These unphased carrier states establish neither linkage
across an ORF nor callability at unreported bases.

The intact-copy analysis keeps tested and untested observations separate.
Among 2,330 exact-`Intact` copy–gene observations carrying targeted variants,
2,166 have at least one unresolved call. Another 10,934 intact copy–gene
observations have no qualifying target and are untested by this panel.

From `Supplementary_Data/Table_S3/ORF_variant_comparison/`, use Python 3.10+
with pysam and Biopython:

```sh
python compare_variant_genotypes.py --targets Data/orf_variant_targets.tsv --raw-dir Data/raw_gt --reference-manifest Data/reference_manifest.tsv --out-prefix reproduced_variant_comparison
python -m unittest -v test_variant_genotypes.py
python quantify_intact_orf_evidence_gaps.py
```

The subtree includes the frozen targets, native sequences, regional VCFs,
reference mappings and call-level results. Its README describes target
rediscovery and the upstream provenance scripts. The old consensus-derived
whole-ORF recovery percentages and zero multi-copy/21% solo-LTR statements are
superseded. The separate Illumina ensemble structural-variant diagnostics in
`Supplementary_Data/Table_S6_short_read_vcf_locus_diagnostics.tsv` still supply
Figure S1; they are not inputs to the current main Table 2.

## Other current reproduction entry points

From `Supplementary_Data/Table_S16/`, run the standard-library verification:

```sh
python3 verify_8q11_23_diversity.py
```

It recalculates diversity from the supplied 583-sequence alignment and checks
the summary, comparison row and haplotype counts. The former assembly-only
script is historical provenance and does not reproduce the corrected result.

From `Supplementary_Data/TypeI_deletion_models/`, use Python with NumPy:

```sh
python reproduce_source_model.py
```

The model supplies Figure 5D and Figure S14. Its weights describe relative
cumulative contributions of deletion classes, not individual ancestral source
loci or their active periods. The model does not estimate ancient
viral-population frequencies.

For Table S13, the retained recurrent-conversion analysis is documented in
`Supplementary_Data/Table_S13/`. For regional phylogenies, tandem-array
nucleotide variation, cassette alignments and terminal-junction measurements,
consult [`analysis/september2026_revision/README.md`](analysis/september2026_revision/README.md)
and the current supplementary source readmes. That workflow retains earlier
numbering and inputs for some panels; it is not an end-to-end reproducer of
this snapshot. In particular, its older solo-LTR inputs and source-model
outputs must not replace the corrected subtrees above.

## Code and archive provenance

| Location | Contents |
|---|---|
| `Supplementary_Data/` | Current manuscript tables, retained sequences, alignments, trees, focused reproduction code and file manifests |
| `Figures/` | Public figure exports; see the figure map for content and source identities |
| `project/manuscript/` | Retained figure builders, artifact filtering, functional refits and historical mechanism models |
| `project/working/biological_orf_annotation_20260802/` | Host-flank-supported locus assignment and biological ORF annotation |
| `project/working/cnv_copy_state_reinterpretation_v1/` and `project/cluster_workflows/cnv_*/` | Copy-state calibration and targeted ONT copy-validation sources |
| `manuscript_figures/python/analysis/` and `manuscript_figures/R/` | Retained analysis and plotting sources; filenames can predate the current figure order |
| `HML2_ProjectResources/cluster_pipeline_source/shared/` | Upstream ORF and flank-anchoring source |
| `analysis/september2026_revision/` | Earlier revision workflows and retained inputs; observe the current-result limits above |
| `data/` | Compact catalog, model and copy-validation evidence retained from the published archive |

The v0.4.2 `HML2_derived_data_v0.4.2.zip` contains that release's catalog,
sequences, supplementary tables and retained inputs. Use matching code and
data archives to reproduce historical releases. Do not overwrite the current
GitHub snapshot with an older archive. The compact-panel runner retains
historical mechanism/helper plots under explicit `Historical_` output names.
The release checks and individual plot runners do not constitute a complete
reproduction of every current manuscript result or figure.

`scripts/restore_corrected_catalog.py` restores the two catalog versions
bundled in `data/`. `REC_CORRECTED` includes the complete Rec stop codon;
`SHORT_ORF_CORRECTED` is the retained input to that documented correction.
`scripts/check_rec_boundary.py` replays that catalog-boundary correction. These
catalog operations do not recreate the current native genotype comparison.

`DATA_REQUIREMENTS.md`, `DATA_SOURCES.md` and `MODEL_NOTES.md` describe retained
workflow inputs, source datasets and model provenance. Historical commands and
model interpretations in those documents must be read against the current
source readmes above. `source_inventory.tsv`, `CHANGES.md` and supplementary
manifests record the distributed source identities and corrections. Neutral
`historical_source/` labels identify provenance, not additional bundled files.
Older archived versions remain available through
[Zenodo](https://doi.org/10.5281/zenodo.22759510).

Original code is MIT-licensed and original derived data are CC BY 4.0.
Third-party sources retain their existing terms.
