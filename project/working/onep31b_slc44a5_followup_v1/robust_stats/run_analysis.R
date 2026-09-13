#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
root <- normalizePath(".", mustWork = TRUE)
base <- file.path(root, "project/working/onep31b_slc44a5_followup_v1/robust_stats")
out <- file.path(base, "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  markers = file.path(root, "project/working/locus_marker_expansion_v1/association_ready_matrix.tsv"),
  states = file.path(root, "project/working/locus_marker_expansion_v1/person_locus_states.tsv"),
  copies = file.path(root, "project/working/array_function_burden_v1/person_copy_function.tsv"),
  expression = file.path(root, "HML2_ProjectResources/data/mage/mage_v1/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz"),
  covariates = file.path(root, "HML2_ProjectResources/data/mage/mage_v1/eQTL_covariates.tab.gz"),
  correction = file.path(root, "HML2_ProjectResources/data/mage/mage_v1/mage_v1_sample_identity_correction.tsv"),
  metadata = file.path(root, "HML2_ProjectResources/data/mage/mage_v1/sample.metadata.MAGE.v1.0.txt"),
  pedigree = file.path(root, "HML2_ProjectResources/data/ref/1kGP.3202_samples.pedigree_info.txt")
)
stopifnot(all(vapply(paths, file.exists, logical(1))))

# Tabix extracts the single pre-specified SLC44A5 TSS row while retaining the
# complete MAGE sample header.
cmd <- sprintf("tabix -h %s chr1:75611115-75611116", shQuote(paths$expression))
expr <- read.delim(pipe(cmd), check.names = FALSE)
stopifnot(nrow(expr) == 1, expr[[4]] == "ENSG00000137968.16")
cov0 <- read.delim(gzfile(paths$covariates), check.names = FALSE)
markers <- read.delim(paths$markers, check.names = FALSE)
states <- read.delim(paths$states, check.names = FALSE)
copies <- read.delim(paths$copies, check.names = FALSE)
corr <- read.delim(paths$correction, check.names = FALSE)
meta <- read.delim(paths$metadata, check.names = FALSE)
ped <- read.table(paths$pedigree, header = TRUE)

raw_axis <- names(expr)[5:ncol(expr)]
stopifnot(identical(raw_axis, names(cov0)[2:ncol(cov0)]))
bio_axis <- raw_axis
swap <- setNames(corr$corrected_biological_identity, corr$raw_mage_row_label)
was_corrected <- bio_axis %in% names(swap)
bio_axis[was_corrected] <- unname(swap[bio_axis[was_corrected]])
stopifnot(length(unique(bio_axis)) == 731)

sex_axis <- factor(as.character(unlist(cov0[cov0$id == "sex", -1], use.names = FALSE)))
cn <- cov0[cov0$id != "sex", , drop = FALSE]
rownames(cn) <- cn$id
cn$id <- NULL
covar <- t(cn)
storage.mode(covar) <- "numeric"
rownames(covar) <- bio_axis

mi <- match(bio_axis, markers$sample)
keep <- which(!is.na(mi))
stopifnot(length(keep) == 39, sum(was_corrected[keep]) == 0)
sample_id <- bio_axis[keep]

st <- states[states$locus == "1p31.1b", ]
st <- st[match(sample_id, st$sample), ]
cp <- copies[copies$locus == "HML-2_1p31.1b", ]
cp <- cp[match(sample_id, cp$sample), ]
stopifnot(all(st$sample == sample_id), all(cp$sample == sample_id))

exposure <- as.numeric(markers[mi[keep], "1p31.1b::internal_presence::all_people"])
dose <- as.numeric(st$internal_element_haplotype_dose)
stopifnot(all(exposure == as.numeric(dose > 0)), sum(exposure) == 13,
          unname(table(factor(dose, levels=0:2))) == c(26, 12, 1))
stopifnot(all(st$retained_element_haplotype_dose == 2),
          all(st$fragment_haplotype_dose == dose),
          all(st$solo_ltr_haplotype_dose + dose == 2),
          all(as.numeric(cp$has_array) == 0),
          all(as.numeric(cp$diploid_physical_copy_number) == dose))

y <- as.numeric(expr[1, 5:ncol(expr)])[keep]
direct_sp <- markers$superpopulation[mi[keep]]
direct_pop <- markers$population[mi[keep]]
meta_i <- match(raw_axis[keep], meta$sample_kgpID)
ped_i <- match(sample_id, ped$sampleID)

d <- data.frame(
  sample = sample_id, raw_mage_label = raw_axis[keep], y = y,
  exposure = exposure, internal_haplotype_dose = dose,
  solo_ltr_haplotype_dose = as.numeric(st$solo_ltr_haplotype_dose),
  retained_haplotype_dose = as.numeric(st$retained_element_haplotype_dose),
  has_array = as.numeric(cp$has_array), superpopulation = factor(direct_sp),
  population = direct_pop, sex = sex_axis[keep],
  batch = factor(meta$batch[meta_i]), RIN = as.numeric(meta$RIN[meta_i]),
  log_num_reads = log(as.numeric(meta$numReads[meta_i])),
  covar[keep, c(paste0("PC", 1:5), paste0("PEER", 1:5)), drop = FALSE],
  father = ped$fatherID[ped_i], mother = ped$motherID[ped_i],
  stringsAsFactors = FALSE, check.names = FALSE
)

rel <- data.frame(sample1 = character(), sample2 = character(), relationship = character())
if (nrow(d) > 1) for (i in 1:(nrow(d)-1)) for (j in (i+1):nrow(d)) {
  p1 <- setdiff(c(d$father[i], d$mother[i]), c("0", 0, NA))
  p2 <- setdiff(c(d$father[j], d$mother[j]), c("0", 0, NA))
  typ <- character()
  if (d$sample[i] %in% p2 || d$sample[j] %in% p1) typ <- c(typ, "parent_child")
  if (length(intersect(p1, p2))) typ <- c(typ, "shared_recorded_parent")
  if (length(typ)) rel <- rbind(rel, data.frame(sample1=d$sample[i], sample2=d$sample[j], relationship=paste(typ, collapse=";")))
}

hc3 <- function(fit, term = "exposure") {
  X <- model.matrix(fit); e <- residuals(fit); h <- hatvalues(fit)
  inv <- solve(crossprod(X))
  V <- inv %*% crossprod(X, X * (e / (1-h))^2) %*% inv
  se <- sqrt(V[term, term]); b <- coef(fit)[term]
  c(beta=unname(b), se=unname(se), statistic=unname(b/se),
    p=unname(2*pt(-abs(b/se), df.residual(fit))))
}

fit_row <- function(id, formula, dat, notes) {
  f <- lm(formula, data=dat)
  z <- coef(summary(f))["exposure", ]; r <- hc3(f)
  data.frame(model_id=id, n=nobs(f), residual_df=df.residual(f),
    beta=z[["Estimate"]], conventional_se=z[["Std. Error"]], conventional_p=z[["Pr(>|t|)"]],
    hc3_se=r[["se"]], hc3_p=r[["p"]], notes=notes, stringsAsFactors=FALSE)
}

primary_formula <- y ~ exposure + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PEER1 + PEER2 + PEER3 + PEER4 + PEER5
models <- rbind(
  fit_row("primary_reproduction", primary_formula, d, "Exact original model"),
  fit_row("pc_only", y ~ exposure + sex + PC1 + PC2 + PC3 + PC4 + PC5, d, "No PEER factors"),
  fit_row("superpopulation_peer", y ~ exposure + sex + superpopulation + PEER1 + PEER2 + PEER3 + PEER4 + PEER5, d, "Explicit superpopulation replaces genotype PCs"),
  fit_row("superpopulation_only", y ~ exposure + sex + superpopulation, d, "No PEER or genotype PCs"),
  fit_row("technical_augmented", y ~ exposure + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PEER1 + PEER2 + PEER3 + PEER4 + PEER5 + batch + RIN + log_num_reads, d, "Adds batch/RIN/read depth; high-dimensional sensitivity"),
  fit_row("dose_primary", primary_formula, transform(d, exposure=internal_haplotype_dose), "0/1/2 internal-fragment haplotype dose"),
  fit_row("exclude_homozygote", primary_formula, d[d$internal_haplotype_dose < 2, ], "Drops sole internal-fragment homozygote")
)

pf <- lm(primary_formula, data=d)
loo <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  f <- lm(primary_formula, data=d[-i, ])
  z <- coef(summary(f))["exposure", ]
  data.frame(omitted_sample=d$sample[i], omitted_exposure=d$exposure[i], beta=z[["Estimate"]],
             se=z[["Std. Error"]], p=z[["Pr(>|t|)"]])
}))
infl <- data.frame(sample=d$sample, exposure=d$exposure, cooks_distance=cooks.distance(pf),
                   exposure_dfbetas=dfbetas(pf)[,"exposure"])

anc <- do.call(rbind, lapply(levels(d$superpopulation), function(s) {
  q <- d[d$superpopulation == s, ]; n1 <- sum(q$exposure == 1); n0 <- sum(q$exposure == 0)
  data.frame(superpopulation=s, n=nrow(q), carriers=n1, solo_ltr_references=n0,
             carrier_fraction=mean(q$exposure), mean_expression_carrier=if(n1) mean(q$y[q$exposure==1]) else NA,
             mean_expression_reference=if(n0) mean(q$y[q$exposure==0]) else NA,
             unadjusted_difference=if(n1 && n0) mean(q$y[q$exposure==1])-mean(q$y[q$exposure==0]) else NA)
}))

set.seed(73115)
tab <- table(d$superpopulation, d$exposure)
ct <- suppressWarnings(chisq.test(tab, simulate.p.value=TRUE, B=20000))
anc_test <- data.frame(test="superpopulation_by_exposure_monte_carlo_chisquare", statistic=unname(ct$statistic),
                       p_value=ct$p.value, replicates=20000)

permute_fl <- function(full_formula, null_formula, dat, B=20000) {
  full <- lm(full_formula, dat); null <- lm(null_formula, dat)
  X <- model.matrix(full); inv <- solve(crossprod(X)); df <- df.residual(full)
  obs <- abs(coef(summary(full))["exposure", "t value"]); groups <- split(seq_len(nrow(dat)), dat$superpopulation)
  exceed <- 0L
  for (b in seq_len(B)) {
    rp <- residuals(null)
    for (ii in groups) rp[ii] <- sample(rp[ii])
    yy <- fitted(null) + rp; cc <- inv %*% crossprod(X, yy); rr <- yy - X %*% cc
    tt <- cc["exposure",] / sqrt(sum(rr^2)/df * inv["exposure","exposure"])
    exceed <- exceed + as.integer(abs(tt) >= obs)
  }
  data.frame(method="Freedman-Lane residual permutation within superpopulation",
    replicates=B, observed_abs_t=obs, exceedances=exceed,
    corrected_p=(exceed+1)/(B+1), seed=73115)
}
perm <- permute_fl(primary_formula,
  y ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PEER1 + PEER2 + PEER3 + PEER4 + PEER5, d)

alignment <- d[, c("sample","raw_mage_label","superpopulation","population","sex","exposure",
  "internal_haplotype_dose","solo_ltr_haplotype_dose","retained_haplotype_dose","has_array","y")]
names(alignment)[names(alignment)=="y"] <- "SLC44A5_inverse_normal_expression"
write.table(alignment, file.path(out,"sample_alignment.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(rel, file.path(out,"recorded_related_pairs.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(models, file.path(out,"model_sensitivities.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(loo, file.path(out,"leave_one_out.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(infl, file.path(out,"influence_diagnostics.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(anc, file.path(out,"ancestry_strata.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="")
write.table(anc_test, file.path(out,"ancestry_exposure_test.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(perm, file.path(out,"blocked_permutation.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

manifest <- data.frame(resource=names(paths), path=unlist(paths), md5=unname(tools::md5sum(unlist(paths))))
write.table(manifest, file.path(out,"source_manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

p0 <- models[models$model_id=="primary_reproduction",]
json <- c("{", sprintf('  "n": %d,', nrow(d)), sprintf('  "carriers": %d,', sum(d$exposure)),
  sprintf('  "solo_ltr_reference_people": %d,', sum(d$exposure==0)),
  sprintf('  "internal_dose_counts": {"0": %d, "1": %d, "2": %d},', sum(d$internal_haplotype_dose==0), sum(d$internal_haplotype_dose==1), sum(d$internal_haplotype_dose==2)),
  sprintf('  "array_carriers": %d,', sum(d$has_array)), sprintf('  "corrected_identity_swaps_in_overlap": %d,', sum(was_corrected[keep])),
  sprintf('  "recorded_related_pairs": %d,', nrow(rel)), sprintf('  "primary_beta": %.12g,', p0$beta),
  sprintf('  "primary_p": %.12g,', p0$conventional_p), sprintf('  "primary_hc3_p": %.12g,', p0$hc3_p),
  sprintf('  "loo_beta_min": %.12g,', min(loo$beta)), sprintf('  "loo_beta_max": %.12g,', max(loo$beta)),
  sprintf('  "blocked_permutation_p": %.12g,', perm$corrected_p),
  '  "contrast": "Type-I internal fragment carriage versus two solo-LTR alleles; not insertion presence/absence",',
  '  "conclusion": "positive association survives prespecified robustness checks but remains exploratory and non-causal"', "}")
writeLines(json, file.path(out,"summary.json"))
cat(sprintf("n=%d; beta=%.6f; HC3 p=%.3g; blocked permutation p=%.3g\n", nrow(d), p0$beta, p0$hc3_p, perm$corrected_p))
