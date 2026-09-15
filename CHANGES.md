# Publication snapshot changes

## 15 September 2026

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
