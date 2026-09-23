# v0.4.3 manuscript corrections, 22 September 2026

- Replace the consensus-derived coding recovery calculation with the native ORF-associated variant comparison and its public inputs. Missing or inadequate genotypes remain unresolved.
- Include all 583 retained 8q11.23 solo-LTR haplotypes and their reference-aligned sequences. The dominant haplotype contains 547 observations.
- Update the deletion-class model to a conservative count of 16 and inclusive count of 17 Delta292-bearing insertions. Retain the earlier helper/mechanism code as historical analyses.
- Include reviewed Figures 2–7. Correct Figure 7C to relative viability and the Figure S13 label recipe.
- Normalize local path prefixes in new text provenance; preserve source and distributed hashes. Native VCF bytes retain their original public-source retrieval headers.
- Preserve the published v0.4.2 archive. The corrected v0.4.3 archives are prepared separately; no new archive DOI exists.

# HML-2 pangenome v0.4.2

This release accompanies the final Table 2 revision of 22 September 2026.

- Add separate intact and combined ORF recovery counts for the matched 282-donor, 83-locus comparison of short-read VCF-derived annotations and long-read calls.
- Separate Type-II Env from the theoretical N-terminally truncated Type-I Env annotation. Translation of the latter product has not been demonstrated.
- Include a standard-library reproduction script, complete donor–locus observations, per-locus counts and input provenance. The original seven aggregate totals reproduce exactly. The frozen input catalogs are unchanged.
- Include the post-v0.4.1 documentation corrections and Figure 6C palette correction already present on GitHub.

All earlier scientific results and public data files are retained. The normal programmed gag/pro/pol frame transitions are permitted by the intact screen. The combined screen additionally retains the caller's stop-free altered-frame candidates; sequence classification alone does not establish protein activity.

No ORF-break dating analysis, private collaborator source data or licensed BioRender artwork is added.

Archive DOI: https://doi.org/10.5281/zenodo.22903723

---

## 0.4.1 — 22 September 2026

This release accompanies the manuscript revision of 22 September 2026.

- Restore the matched short-read reconstruction comparison as main Table 2. The manuscript retains seven main figures. The analysis includes 282 donors and 83 loci after artifact filtering.
- Preserve the catalog inclusion flag in the comparison and exclude rejected assembly records. Add a regression test for an excluded copy that would otherwise create a false coding carrier.
- Include compact input catalogs and exact reproduction commands for Table S6.
- Move the copy-count legend below the array panel (Figure 2B).
- Give each category in the cassette figure a distinct color.
- Complete missing population metadata for four donors in the 8q11.23 network. All 461 observations and 13 sequence haplotypes are unchanged. Restore the accepted two-panel design in the portable renderer.
- Update final figure references after the restored main figure and removal of four redundant supplementary displays. Their source data remain available.

All earlier scientific results are retained. No ORF-break dating is included. Private collaborator data and licensed BioRender artwork are excluded from the public archive.

Archive DOI: https://doi.org/10.5281/zenodo.22888353

## 21 September 2026 — v0.4.0

Code and supporting data for the manuscript revision of 21 September 2026.

- Nucleotide and gap counts for 617 copies from 289 tandem arrays, with conditional duplication-time calculations.
- Type I cassette comparisons using matched Type II subfamilies, human and ape alignments, and the complete supplementary sequence alignment.
- Solo-LTR diversity comparisons and expected nucleotide differences.
- Phylogenetic bootstrap support and nucleotide-difference quantification for the revised Figure 3.
- Terminal-junction evidence with candidate TSD differences evaluated across 4–6-base lengths. The unsupported 3q12.3 paired-TSD population result is withdrawn.
- The public ONT read-source index and updated provenance.
- Portable reproduction commands and retained sequence inputs for the revised analyses.

The corrected catalog and earlier review corrections are retained. This release adds no ORF-break dating analysis. Original code is MIT licensed. Original derived data are CC BY 4.0. Third-party sources retain their existing terms.

# Publication snapshot changes

## Release 0.3.0, 16 September 2026

- Archived the corrected code, complete testing families and Rec catalog at
  https://doi.org/10.5281/zenodo.22783894.
- Retained the previously archived public HML-2 sequence and phenotype inputs.
  Replaced superseded supplementary tables and catalog copies with the current
  release. Version 0.2.0 remains separately accessible as the earlier snapshot.
- Refreshed three input identity records to match the distributed source files,
  including normalized provenance and the corrected 7p22.1 missing call.

## 15 September 2026

- Corrected the Type-II Rec second-exon endpoint from 8466 to 8467
  (zero-based, half-open), retaining the complete stop codon. Recalled all
  36,073 Rec-bearing records from the existing alignments. The unchanged-boundary
  replay had no unexplained discrepancies. Paired call files account for the
  three alignments regenerated after the older catalog snapshot.
- Removed false terminal-frameshift labels while retaining sequence-supported
  frameshifts and premature stops. Updated Rec counts and Figure 4. Functional
  exposures and P/q values did not change after the Rec correction.
- Included minimal alignment slices, reference sequence, paired calls and a
  complete public replay test. `REC_CORRECTED` is now the analysis catalog.
  `SHORT_ORF_CORRECTED` remains only as the replay input.
- Applied the existing short-product rule to formerly exempt `Intact_FS_End`
  calls and synchronized coding summaries. Type-I Env calls are unchanged.
- Refit functional burdens from the current resolved catalog. Corrected q-value
  reporting to 174 finite P values among 237 attempted models, preserving all
  missing tests as missing. Added complete MAGE and GEUVADIS testing families.
- Kept the uncalled 7p22.1 haplotype missing. The array numerator is 247 among
  583 called haplotypes, from 584 eligible haplotypes. The two known zero-copy
  haplotypes are solo-LTR calls.
- Restricted the 10q24.2 regional nearest-neighbor comparison to 64 loci shared
  between the LTR and Pol trees. Existing phylogenetic trees were not refit.
- Published the corrected catalog as a compressed derived-data file in GitHub.
  These changes are not present in the immutable 14 September Zenodo archive.
- Updated the Fiber-seq source to the collaborator-confirmed minimum of ten
  reads per peak and removed the redundant peak flag filter. Collaborator data
  and artwork remain excluded from the public repository.
- Added the coordinate-matched Figure S17 locus crosswalk. Assay labels and
  reference intervals are preserved beside the catalog identifiers.

- Matched the copy-review figure label to the current manuscript, using
  "Other supported copy" in place of "Later duplication". No calls, counts,
  model results or archived data changed. The archived v0.2.0 release remains
  available unchanged.

## Release 0.2.0, 14 September 2026

- Updated the manuscript title to specify HERV-K(HML-2).
- Incorporated the resolved catalog, 2,235 public acrocentric assignments,
  recovered source calls, and exclusion of the two HML-11 comparator loci.
  Figure 1 now covers 99 physical loci and four separately counted copy groups.
- Replaced the distributed supplementary data with the current numbered tables,
  alignment/tree inputs, direct Type-I calls and source-sequence evidence.
  The 36 previously unresolved internal-bearing calls are now assigned Type I,
  giving 9,733 Type-I calls and no Type-II calls at the 20 surveyed Type-I loci.
- Added the current Figure 3 builder, including chromosome 4 and the distinct
  exact nucleotide-sharing versus regional host-flank comparisons.
- Added all 54 recurrent-conversion/drift scenarios and portable reproduction
  code, with the frozen input table and numerical checks.
- Updated the retained R analyses for tandem arrays, variant-position counts
  and Fiber-seq summaries. Fiber-seq data and artwork are not redistributed.
- Packaged full derived sequence and table inputs separately from Git.
- Applied MIT to original code and CC BY 4.0 to original derived data, as
  authorized by the author. Source-specific third-party terms remain intact.

## Earlier snapshot

- Corrected Figure 1 and Tables S1/S5 using one shared observation summary.
  Empty physical-locus cells are Unknown. Pipeline noncarrier sentinels are
  labelled Noncarrier call, and unlocalized copy groups are not ranked as
  single loci. The descriptive panel uses 584 autosomal, 432 X and 144 Y
  copies, with unknown-sex donors excluded only from sex-chromosome counts.
  Added compact sex/partition inputs and focused conservation tests.
  The bundled source-path copy of S5 is checked against the distributed table.
- Recovered the original NucFreq reducer and documented its >=3 secondary
  reads and >=0.15 top-two base-fraction rule. Added executable checks of all
  416 regional depth-reliability flags and 351 numeric site fractions. The
  frozen regional counts are unchanged. Historical numeric-window semantics
  are retained explicitly rather than silently shifted.
- Preserved retained analysis source from the local research workspaces in a
  new repository, without copying their Git history or operational logs.
- Updated the main figure builder to retain only the structure and measured
  solo-LTR variation at 8q11.23. The network now uses the checked 461-sequence
  source table with 13 haplotypes and a 432-sequence major haplotype. The verified
  CACAC target-site duplications are labelled at both proviral flanks.
- Replaced the Figure 3 chronology schematic with an un-timed duplicated-group
  schematic. The phylogeny and chromosome-location code remains in place.
- Kept only the completed helper and neutral-founder sensitivity panels in S13,
  with panels relettered A and B.
- Kept only the nearest-neighbor, linked-substitution and orthology-aware
  deletion-count panels in S17, with panels relettered A, B and C.
- Reduced S20 to the observed Type-I/Type-II counts. Zero Type-II calls are not
  plotted as a nonzero point on a logarithmic axis.
- Changed Figure 5 and S14A labels to Bayes factor, matching marginal-evidence
  integration. Changed S14D to model weight under equal priors.
- Changed the S4A axis to frequency meeting ORF criteria. Its original six
  accepted status categories and exclusion of Undetermined calls are unchanged.
- Clarified Figure 2A/B as locus copy counts and frequencies. Preserved every
  count and observed color, removing only the unused five-copy legend entry.
  Documented the separate 14q11.2 fragment and S4B's reference-inclusive
  array-bearing-haplotype denominator, without changing S4B values.
- Added a compact, directly runnable helper-sensitivity calculation instead of
  a dependency on the former broader simulation workflow.
- Retained the completed targeted ONT copy validation and published short-read
  population comparisons.
- Converted source roots and installation-specific metadata to repository-
  relative or explicit example paths. Updated the figure-style description to
  the current manuscript target without altering the retained color scheme.
- Updated R configuration loading to work from the new repository root.
