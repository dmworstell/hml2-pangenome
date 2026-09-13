options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(".")
result_root <- file.path(root, "project/working/direct_ebv_fitness_screen_v1/results")
read_tsv <- function(path) read.delim(path, sep = "\t", quote = "", check.names = FALSE,
                                       na.strings = c("NA"), comment.char = "")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
                                            quote = FALSE, na = "NA")
compatible <- c("Intact", "Frameshift_at_end", "Intact_FS_End")

cluster_fit <- function(data, exposure, outcome, adjustment_ids = character(),
                        ancestry = c("superpopulation", "population")) {
  ancestry <- match.arg(ancestry)
  d <- data.frame(y = data[[outcome]], x = data[[exposure]], sex = factor(data$sex),
                  ancestry = factor(if (ancestry == "superpopulation") data$superpopulation else data$population),
                  cluster = data$pedigree_component)
  if (length(adjustment_ids)) {
    for (j in seq_along(adjustment_ids)) d[[paste0("adj", j)]] <- data[[adjustment_ids[j]]]
  }
  d <- d[complete.cases(d), ]
  if (nrow(d) < 20L || length(unique(d$x)) < 2L || sd(d$x) == 0) {
    return(list(status = "not_estimable", n = nrow(d)))
  }
  adjustment_terms <- if (length(adjustment_ids)) paste0("adj", seq_along(adjustment_ids)) else character()
  fit <- lm(as.formula(paste("y ~", paste(c("x", adjustment_terms, "ancestry", "sex"), collapse = " + "))), data = d)
  X <- model.matrix(fit)
  pivot <- fit$qr$pivot[seq_len(fit$rank)]
  Xr <- X[, pivot, drop = FALSE]
  if (!("x" %in% colnames(Xr))) return(list(status = "exposure_aliased", n = nrow(d)))
  bread <- tryCatch(solve(crossprod(Xr)), error = function(e) NULL)
  if (is.null(bread)) return(list(status = "singular", n = nrow(d)))
  ix <- match("x", colnames(Xr))
  residual <- residuals(fit)
  cluster_rows <- split(seq_len(nrow(d)), d$cluster)
  scores <- lapply(cluster_rows, function(rows) colSums(Xr[rows, , drop = FALSE] * residual[rows]))
  meat <- Reduce(`+`, lapply(scores, tcrossprod))
  clusters <- length(scores)
  correction <- (clusters / (clusters - 1)) * ((nrow(d) - 1) / (nrow(d) - fit$rank))
  vc <- bread %*% meat %*% bread * correction
  se <- sqrt(vc[ix, ix])
  beta <- coef(fit)["x"]
  p <- 2 * pt(abs(beta / se), df = clusters - 1, lower.tail = FALSE)
  list(status = "fit", n = nrow(d), beta = beta, se = se, p = p, clusters = clusters,
       ci_low = beta + qt(0.025, clusters - 1) * se,
       ci_high = beta + qt(0.975, clusters - 1) * se)
}

fit_row <- function(fit, label, ancestry) {
  getv <- function(field) if (is.null(fit[[field]])) NA else fit[[field]]
  data.frame(label = label, ancestry_model = ancestry, status = fit$status, n = fit$n,
             beta = getv("beta"), se = getv("se"), ci95_low = getv("ci_low"),
             ci95_high = getv("ci_high"), p_pedigree_cluster = getv("p"),
             clusters = getv("clusters"), stringsAsFactors = FALSE)
}

# Build locus-resolved functional-unit and physical-unit truth directly from the
# authoritative catalog. Absent catalog rows contribute zero units.
catalog <- read_tsv("HML2_ProjectResources/data/catalog/combined_hml2_orf_analysis.tsv")
matrix <- read_tsv(file.path(result_root, "person_level_direct_matrix.tsv"))
direct_ids <- matrix$sample
catalog <- catalog[catalog$ID %in% direct_ids, ]
prov <- catalog[catalog$Structure %in% c("Provirus", "Provirus_from_Multi") &
                  catalog$provirus_type %in% c("type1", "type2"), ]
locus_types <- aggregate(provirus_type ~ Locus, prov, function(z) paste(sort(unique(z)), collapse = ";"))
stopifnot(all(locus_types$provirus_type %in% c("type1", "type2")))
type1_loci <- sort(locus_types$Locus[locus_types$provirus_type == "type1"])
stopifnot(length(type1_loci) == 17L)

make_unit_matrix <- function(rows, loci) {
  out <- matrix(0, nrow = length(direct_ids), ncol = length(loci),
                dimnames = list(direct_ids, loci))
  if (nrow(rows)) {
    observed <- as.data.frame.matrix(xtabs(~ ID + Locus, rows))
    common_ids <- intersect(rownames(observed), rownames(out))
    common_loci <- intersect(colnames(observed), colnames(out))
    out[common_ids, common_loci] <- as.matrix(observed[common_ids, common_loci, drop = FALSE])
  }
  out
}

type1_np9_rows <- prov[prov$provirus_type == "type1" & prov$np9 %in% compatible, ]
type1_np9_by_locus <- make_unit_matrix(type1_np9_rows, type1_loci)
stopifnot(all(rowSums(type1_np9_by_locus) == matrix$type1_compatible_np9_units))

# Houldcroft: decompose the catalog-wide Type-I Np9 signal by locus and audit
# leave-one-person-out stability.
houldcroft_dir <- file.path(result_root, "Houldcroft2014_EBV_qPCR")
h <- read_tsv(file.path(houldcroft_dir, "person_level_outcome_matrix.tsv"))
np9_columns <- paste0("np9__", make.names(type1_loci))
np9_frame <- data.frame(sample = direct_ids, type1_np9_by_locus, check.names = FALSE)
names(np9_frame)[-1] <- np9_columns
h <- merge(h, np9_frame, by = "sample", all.x = TRUE)
stopifnot(nrow(h) == 49L)

decomp <- list()
for (i in seq_along(type1_loci)) {
  locus <- type1_loci[i]
  locus_col <- np9_columns[i]
  h$leave_locus_out_np9 <- h$type1_compatible_np9_units - h[[locus_col]]
  fits <- list(
    cluster_fit(h, "leave_locus_out_np9", "log2_relative_ebv_copy_number_qpcr", "type1_units", "superpopulation"),
    cluster_fit(h, "leave_locus_out_np9", "log2_relative_ebv_copy_number_qpcr", "type1_units", "population"),
    cluster_fit(h, locus_col, "log2_relative_ebv_copy_number_qpcr", "type1_units", "superpopulation"),
    cluster_fit(h, locus_col, "log2_relative_ebv_copy_number_qpcr", "type1_units", "population")
  )
  modes <- c("leave_locus_out", "leave_locus_out", "single_locus", "single_locus")
  ancestries <- c("superpopulation", "population", "superpopulation", "population")
  for (j in seq_along(fits)) {
    row <- fit_row(fits[[j]], locus, ancestries[j])
    row$mode <- modes[j]
    row$locus_unit_sum <- sum(h[[locus_col]], na.rm = TRUE)
    row$locus_unit_carriers <- sum(h[[locus_col]] > 0, na.rm = TRUE)
    row$locus_unit_unique_values <- length(unique(h[[locus_col]]))
    row$locus_unit_counts <- paste(names(table(h[[locus_col]])), as.integer(table(h[[locus_col]])),
                                   sep = ":", collapse = ";")
    row$correlation_with_total_np9 <- if (sd(h[[locus_col]]) > 0) {
      cor(h[[locus_col]], h$type1_compatible_np9_units)
    } else NA_real_
    decomp[[length(decomp) + 1L]] <- row
  }
}
decomp <- do.call(rbind, decomp)
write_tsv(decomp[, c("label", "mode", "ancestry_model", "locus_unit_sum", "locus_unit_carriers",
                     "locus_unit_unique_values", "locus_unit_counts", "correlation_with_total_np9",
                     "status", "n", "beta", "se", "ci95_low", "ci95_high",
                     "p_pedigree_cluster", "clusters")],
          file.path(houldcroft_dir, "type1_np9_locus_decomposition.tsv"))

h_complete <- h[complete.cases(h[, c("log2_relative_ebv_copy_number_qpcr", "type1_compatible_np9_units",
                                     "type1_units", "sex", "superpopulation", "population",
                                     "pedigree_component")]), ]
np9_loo <- lapply(h_complete$sample, function(sid) {
  d <- h[h$sample != sid, ]
  primary <- cluster_fit(d, "type1_compatible_np9_units", "log2_relative_ebv_copy_number_qpcr",
                         "type1_units", "superpopulation")
  population <- cluster_fit(d, "type1_compatible_np9_units", "log2_relative_ebv_copy_number_qpcr",
                            "type1_units", "population")
  data.frame(dropped_sample = sid, dropped_np9_units = h$type1_compatible_np9_units[h$sample == sid],
             dropped_outcome = h$log2_relative_ebv_copy_number_qpcr[h$sample == sid],
             beta_superpopulation = primary$beta, p_superpopulation = primary$p,
             beta_population = population$beta, p_population = population$p)
})
write_tsv(do.call(rbind, np9_loo), file.path(houldcroft_dir, "type1_np9_leave_one_person_out.tsv"))

# Im: audit the 6q14.1 structural association by population and individual, and
# determine whether the general Type-II physical burden persists without 6q14.1.
im_dir <- file.path(result_root, "Im2012_intrinsic_growth")
im <- read_tsv(file.path(im_dir, "person_level_outcome_matrix.tsv"))
six_id <- "HML-2_6q14.1::structural::solo_ltr_vs_any_provirus"
stopifnot(nrow(im) == 34L, sum(im[[six_id]] > 0, na.rm = TRUE) == 7L)

group_summary <- aggregate(im$intrinsic_growth_rate_per_10000,
                           by = list(exposure = im[[six_id]], superpopulation = im$superpopulation,
                                     population = im$population),
                           FUN = function(z) c(n = length(z), mean = mean(z), sd = if (length(z) > 1) sd(z) else NA))
group_summary <- data.frame(exposure = group_summary$exposure,
                            superpopulation = group_summary$superpopulation,
                            population = group_summary$population,
                            n = group_summary$x[, "n"], mean_growth_per_10000 = group_summary$x[, "mean"],
                            sd_growth_per_10000 = group_summary$x[, "sd"])
write_tsv(group_summary, file.path(im_dir, "growth_6q14_group_summary.tsv"))

im_complete <- im[complete.cases(im[, c("intrinsic_growth_rate_per_10000", six_id, "sex",
                                        "superpopulation", "population", "pedigree_component")]), ]
six_loo <- lapply(im_complete$sample, function(sid) {
  d <- im[im$sample != sid, ]
  primary <- cluster_fit(d, six_id, "intrinsic_growth_rate_per_10000", character(), "superpopulation")
  population <- cluster_fit(d, six_id, "intrinsic_growth_rate_per_10000", character(), "population")
  data.frame(dropped_sample = sid, dropped_exposure = im[[six_id]][im$sample == sid],
             dropped_outcome = im$intrinsic_growth_rate_per_10000[im$sample == sid],
             beta_superpopulation = primary$beta, p_superpopulation = primary$p,
             beta_population = population$beta, p_population = population$p)
})
write_tsv(do.call(rbind, six_loo), file.path(im_dir, "growth_6q14_leave_one_person_out.tsv"))

type2_loci <- sort(locus_types$Locus[locus_types$provirus_type == "type2"])
physical_by_locus <- make_unit_matrix(prov[prov$provirus_type == "type2", ], type2_loci)
stopifnot("HML-2_6q14.1" %in% colnames(physical_by_locus))
six_physical <- data.frame(sample = direct_ids,
                           sixq14_physical_units = physical_by_locus[, "HML-2_6q14.1"])
im <- merge(im, six_physical, by = "sample", all.x = TRUE)
im$type2_units_without_6q14 <- im$type2_units - im$sixq14_physical_units
burden_audit <- rbind(
  fit_row(cluster_fit(im, "type2_units", "intrinsic_growth_rate_per_10000", character(), "superpopulation"),
          "all_type2_physical_units", "superpopulation"),
  fit_row(cluster_fit(im, "type2_units", "intrinsic_growth_rate_per_10000", character(), "population"),
          "all_type2_physical_units", "population"),
  fit_row(cluster_fit(im, "type2_units_without_6q14", "intrinsic_growth_rate_per_10000", character(), "superpopulation"),
          "type2_physical_units_without_6q14", "superpopulation"),
  fit_row(cluster_fit(im, "type2_units_without_6q14", "intrinsic_growth_rate_per_10000", character(), "population"),
          "type2_physical_units_without_6q14", "population")
)
write_tsv(burden_audit, file.path(im_dir, "type2_burden_6q14_exclusion_audit.tsv"))

cat("Houldcroft Type-I Np9 leave-one-person-out ranges:\n")
print(sapply(do.call(rbind, np9_loo)[, c("beta_superpopulation", "p_superpopulation",
                                        "beta_population", "p_population")], range))
cat("Im 6q14.1 leave-one-person-out ranges:\n")
print(sapply(do.call(rbind, six_loo)[, c("beta_superpopulation", "p_superpopulation",
                                        "beta_population", "p_population")], range))
cat("Im Type-II burden exclusion audit:\n")
print(burden_audit[, c("label", "ancestry_model", "beta", "p_pedigree_cluster")])
