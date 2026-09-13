# Compact supporting data

`cnv_depth_summary.tsv` is the retained 416-region copy-validation summary.
351 rows have numeric `nucfreq_het_frac` values. That field equals the number of
sites counted as heterozygous divided by the number of evaluated sites, rounded
to four decimals. It is not a median second-base fraction. All 351 numeric
fractions are below 0.05, with a maximum of 0.0135. The original July reducer
was recovered. A qualifying site has at least three reads for its second-most
abundant base and a second-base fraction of at least 0.15 using the two most
abundant base counts. Low-baseline regions retain their NucFreq counts, while
their depth-support ratios are omitted. The source provenance and historical
coordinate handling are documented in
`project/cluster_workflows/cnv_v2/NucFreq_method_provenance.md`.

`Figure_1_structural_source_data.tsv` and
`Figure_1_sex_chromosome_eligibility.tsv` contain the corrected structural
observations and chromosome eligibility. Missing records are Unknown, not
noncarrier calls. The 101 physical-locus labels use 584 autosomal, 432 X and
144 Y copies. Four unlocalized copy groups have observation counts but no
single-locus frequency denominator. Compact sex and X/Y partition inputs are
under `project/inputs/`.

`8q11_23_provirus_annotation.json` records the sole retained provirus row used to
check the two CACAC target-site duplications in Figure 7. Its full catalog input
is identified in `external_inputs.tsv`.

`figure_2_s4_array_scope.json` records the checked copy-count summaries and
specific catalog entries needed to test the Figure 2/S4B definition difference.
It includes the exact input identities used for this derived regression fixture.

`additional_manifest.tsv` identifies the supplementary compact model and CNV
inputs. Extracted JSON blocks contain only the fields consumed by retained
figure builders. Original source checksums are recorded separately from the
distributed-file checksums.
