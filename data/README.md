# Compact supporting data

`cnv_depth_summary.tsv` is the retained 416-region copy-validation summary.
351 rows have numeric `nucfreq_het_frac` values. That field equals the number of
sites counted as heterozygous divided by the number of evaluated sites, rounded
to four decimals. It is not a median second-base fraction. All 351 numeric
fractions are below 0.05, with a maximum of 0.0135. The original rule for calling
an individual site heterozygous is not established by this summary table alone.

`8q11_23_provirus_annotation.json` records the sole retained provirus row used to
check the two CACAC target-site duplications in Figure 7. Its full catalog input
is identified in `external_inputs.tsv`.

`additional_manifest.tsv` identifies the supplementary compact model and CNV
inputs. Extracted JSON blocks contain only the fields consumed by retained
figure builders. Original source checksums are recorded separately from the
distributed-file checksums.
