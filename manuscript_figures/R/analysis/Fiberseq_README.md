# Fiber-seq analysis and sample counts

`Fiberseq_haplotype_reanalysis.R` supplies the analysis for final Figure S21.
Historical `Figure_S17_*` output names are retained as provenance. The figure
map at the repository root records the final manuscript labels.

The 21 collaborator peak files contain 39 distinct donors and 45 donor/run
groups. Each group has two haplotypes, represented by 90 distinct `sample_id`
values. Filtering `is_primary_sample=true` retains 39 donor/run groups and
78 run-haplotype sample IDs. All 39 biological donors remain. The removed
records are alternate technical runs, not additional donors.

| Stage | Peak rows | Biological donors | Run-haplotype sample IDs | Haplotype–locus combinations |
|---|---:|---:|---:|---:|
| All input records | 3,960 | 39 | 90 | 1,638 |
| Primary records | 3,432 | 39 | 78 | 1,638 |
| Coverage ≥10, all runs | 3,243 | 39 | 90 | 1,358 |
| Coverage ≥10, primary runs | 2,782 | 39 | 78 | 1,344 |

The 1,344 retained haplotype–locus observations span 21 loci. Each locus has
26–77 haplotypes from 22–39 donors. Every observation has exactly one primary
run. Peak FIRE counts and total coverage are summed within each primary
haplotype and locus before taking their ratio.

A source recount on 21 September 2026 reproduced the total coverage, FIRE
coverage, peak count and pooled fraction for all 1,344 observations exactly.
All 21 input files matched the retained input-manifest checksums. The source
files contain 382 records with both coverage fields missing. No finite FIRE
coverage is negative or greater than total coverage. The nonnegative check in
the script is an input-validity check, not a biological exclusion criterion.
The primary filter and its resulting values required no correction.

Reproduction requires the collaborator peak files, which are not redistributed:

```sh
Rscript manuscript_figures/R/analysis/Fiberseq_haplotype_reanalysis.R \
  --peaks-dir /path/to/collaborator/peak-files \
  --out-dir outputs/fiberseq_reanalysis
```

The default minimum peak coverage is 10. The coordinate-matched locus crosswalk
is supplied beside the script. Donor-level evidence and collaborator peak rows
remain outside this public repository.
