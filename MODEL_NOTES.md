# Retained analysis definitions

## Locus counts and tandem-array sizes in Figures 2 and S4B

Figure 2A shows locus copy counts among carrier haplotypes. Figure 2B divides
the count in each displayed multicopy class by all 584 population haplotypes.
The observed multicopy classes are 2, 3, 4 and 6. No five-copy class is observed.
At 14q11.2, NA20282 h2 has three locus copies: one separate fragment and two
tandem proviruses. Its Figure 2A/B locus count is three, not a three-provirus
tandem. The separate fragment does not enter Figure 2C–F, whose array-unit
filter includes only the two `_part` proviruses. Their 290-array denominator
and within-array results are unchanged.

Figure S4B instead shows the distribution of tandem-array sizes among
array-bearing haplotypes at each locus, using the maximum `_part` index. Its
7p22.1 denominator is 249 array-bearing haplotypes, including the GCA and CHM13
reference arrays. These two two-copy arrays are not population observations.
The corresponding population-only count in Figure 2 is 247 arrays. At 14q11.2,
S4B shows the two-copy tandem in its one array-bearing haplotype, excluding the
separate fragment from the array-size count. S4B source values are preserved.

## Positional ORF criterion in Supplementary Figure S4A

The positional tandem-array panel accepts the historical status values
`intact`, `no_stop`, `no_stop_fs_end`, `frameshift_at_end`, `intact_fs_end` and
`intact_fs_end_premature_stop`, after lowercasing. `Undetermined` calls are
excluded from each mean. This is the plotted frequency meeting ORF criteria,
not an Intact-only frequency or a measured functional fraction. It is distinct
from the three-category sequence-compatible burden used in the Np9 refit.

## Helper availability

For a focal fraction `f` and mean number of viral genomes per cell `m`, helper
availability is `1 - exp[-m(1-f)]`. The minimum compensating replication gain is
its reciprocal. The two focal fractions are `10/16` and `17/23`, and the evaluated
values of `m` are `0.5, 1, 2, 3, 5, 10`.

For `c` starting copies in an effective viral population of size `N`, the neutral
frequency-martingale upper bound for ever reaching `f` is `min(1,c/(N*f))`.
The code evaluates `c=1,5` and `N=10,50,100,1000,10000,100000`. These are
analytical sensitivity bounds, not fitted biological parameters.

## Source opportunity versus a focal effect

The two count vectors are `[10,1,1,1,1,1,1]` and `[17,1,1,1,1,1,1]`. The seven
source probabilities have a symmetric Dirichlet prior with per-class alpha
`0.1,0.25,0.5,1,2,5,20`. The plotted heterogeneity scale `1/sqrt(alpha)` is the CV
of the underlying Gamma source weights, not the exact marginal CV of normalized
Dirichlet probabilities.

The focal-effect model adds a log-uniform multiplier from 1 to 50. Both models
integrate over a log-uniform focal ascertainment multiplier, either 0.67–1.5 or
0.25–4.0. Marginal evidence is integrated deterministically using Gauss–Legendre
quadrature, with 192 probability nodes and 80 nodes per log-uniform prior.
Evidence ratios are Bayes factors. The source tables retain their historical
column name `likelihood_ratio_effect_vs_source_only` for compatibility, but the
figures label the quantity correctly.

The separate gain-budget model averages likelihoods over its explicit component
priors. Its `equal_model_prior_share` is a normalized model-evidence weight
under equal prior model probabilities. It is not a count of retained simulation
draws. No component weight establishes that the corresponding mechanism occurs.
