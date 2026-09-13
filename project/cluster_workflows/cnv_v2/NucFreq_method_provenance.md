# NucFreq summary method

The retained regional table is byte-identical to the original July 5, 2026 output, SHA-256 `566cec9fc714c5b3f37f6cd892d800c87bcd46cd9b5b66fbe5a17497d0d12923`.

The source is the original `04_cnv_summary.py`, July 5 copy SHA-256 `a8a7d560fc19c68d232c778983345a9c1a6a3e4050a03032c35049f5be317abe`. Its NucFreq parsing and site-counting functions are identical to the June 12 copy, SHA-256 `78b1ba732ab21f3bd31ec8df8b97e19ef66facb72957147b18772f6fb30a78c1`.

At each covered position the four A/C/G/T counts are ranked. If the two largest counts are d1 and d2, a mixed-base site requires d2 >= 3 and d2 / (d1 + d2) >= 0.15. The regional fraction is the number of mixed-base sites divided by the number of positions with d1 + d2 > 0. There is no separate minimum-total-depth threshold. A zero denominator produces a missing fraction, not zero.

Per-sample/per-locus fractions sum the site numerators and denominators across regions. The cross-locus figure plots those sample-level fractions and their per-locus medians. A 0.05 regional-burden reference is separate from the 0.15 per-site base-count threshold.

Depth support is body median / (flank median / 2), with 0.6 as the minimum support value. The sample reference is the median of its positive per-region flank medians. A region is low baseline if its flank median is nonpositive or below 0.4 times that reference. The original code clears depth ratios for those regions, but retains their NucFreq summaries.

The retained table contains 416 regions, 381 marked NucFreq available, and 351 numeric regional fractions. All 351 are below 0.05, with maximum 0.0135. Thirty-six of those numeric rows are flagged low baseline. The public checker recomputes all recorded fractions, reference medians, low-baseline flags, and available all-read haploid-copy depth ratios from the table.

## Historical interval semantics

The original reducer derives body bounds from the minimum and maximum positions in the depth track by removing 10,000 from each flank, using the full observed span when the resulting lower bound is not below the upper bound. It applies those numeric bounds inclusively to the NucFreq BED start column. The reproduction module deliberately preserves this behavior and does not silently move boundaries.

Depth positions are one-based, while the original BED documentation describes zero-based starts. Therefore a strict coordinate conversion may shift the historical window by one base. For a file with one row per genomic position, changing [b,e] to [b-1,e-1] removes e and adds b-1. An upper bound for a currently numeric region is (H+1)/(N-1), where H and N are its retained site counts. Across the 351 currently numeric rows this bound is at most 0.018182. This is a conditional boundary-sensitivity calculation, not a raw-data recount. It does not address currently missing fractions or duplicate-coordinate BED records.

## Invocation evidence

The recovered code and contemporaneous README specify defaults `NUC_HET_THRESH=0.15` and `NUC_HET_MIN_ALT=3`. Both could be overridden by environment variables. No saved override has been found. The original July 5 rendered cross-locus figure explicitly labels the NucFreq axis as "alt >= 15%" and reports 348/348 depth-supported regions. Its PNG SHA-256 is `3a5216d42cd5654e27bd373950a9e2ad6ea883ba5581b9bdaddcc55caee7ec52`. This corroborates the per-site fraction threshold in the original output, alongside the preserved site-counting code. A raw-BED recount is distinct from the table arithmetic checks and is not claimed here.
