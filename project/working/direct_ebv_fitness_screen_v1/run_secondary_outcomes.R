options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(".")
result_root <- file.path(root, "project/working/direct_ebv_fitness_screen_v1/results")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(path) read.delim(path, sep = "\t", quote = "", check.names = FALSE,
                                       na.strings = c("NA"), comment.char = "")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
                                            quote = FALSE, na = "NA")

matrix <- read_tsv(file.path(result_root, "person_level_direct_matrix.tsv"))
manifest <- read_tsv(file.path(result_root, "feature_manifest.tsv"))
stopifnot(nrow(matrix) == 292L, nrow(manifest) == 79L)

outcomes <- list(
  Houldcroft2014_EBV_qPCR = list(
    path = "project/working/data_search_direct_cellular_phenotypes_v1/derived/Houldcroft2014_EBV_qPCR_unique_direct292_join.tsv",
    raw = "relative_ebv_copy_number_qpcr_mean",
    transform = function(x) log2(x),
    transformed = "log2_relative_ebv_copy_number_qpcr",
    expected_n = 49L,
    units = "log2 relative qPCR EBV copy number"
  ),
  Im2012_intrinsic_growth = list(
    path = "project/working/data_search_direct_cellular_phenotypes_v1/derived/Im2012_LCL_intrinsic_growth_direct292_join.tsv",
    raw = "intrinsic_growth_rate",
    transform = function(x) x / 10000,
    transformed = "intrinsic_growth_rate_per_10000",
    expected_n = 34L,
    units = "10,000 published intrinsic-growth-rate units"
  )
)

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
  if (nrow(d) < 20L || length(unique(d$x)) < 2L || sd(d$x) == 0) return(list(status = "not_estimable", n = nrow(d)))
  adjustment_terms <- if (length(adjustment_ids)) paste0("adj", seq_along(adjustment_ids)) else character()
  fit <- lm(as.formula(paste("y ~", paste(c("x", adjustment_terms, "ancestry", "sex"), collapse = " + "))), data = d)
  X <- model.matrix(fit)
  pivot <- fit$qr$pivot[seq_len(fit$rank)]
  Xr <- X[, pivot, drop = FALSE]
  if (!("x" %in% colnames(Xr))) return(list(status = "exposure_aliased", n = nrow(d)))
  bread <- tryCatch(solve(crossprod(Xr)), error = function(e) NULL)
  if (is.null(bread)) return(list(status = "singular", n = nrow(d)))
  ix <- match("x", colnames(Xr))
  beta <- coef(fit)["x"]
  residual <- residuals(fit)
  cluster_rows <- split(seq_len(nrow(d)), d$cluster)
  scores <- lapply(cluster_rows, function(rows) colSums(Xr[rows, , drop = FALSE] * residual[rows]))
  meat <- Reduce(`+`, lapply(scores, tcrossprod))
  clusters <- length(scores)
  correction <- (clusters / (clusters - 1)) * ((nrow(d) - 1) / (nrow(d) - fit$rank))
  vc <- bread %*% meat %*% bread * correction
  se <- sqrt(vc[ix, ix])
  stat <- beta / se
  p <- 2 * pt(abs(stat), df = clusters - 1, lower.tail = FALSE)
  list(status = "fit", n = nrow(d), beta = beta, se = se, p = p, clusters = clusters,
       ci_low = beta + qt(0.025, clusters - 1) * se,
       ci_high = beta + qt(0.975, clusters - 1) * se,
       standardized_beta = beta * sd(d$x) / sd(d$y))
}

all_models <- list()
all_descriptives <- list()
for (outcome_name in names(outcomes)) {
  spec_outcome <- outcomes[[outcome_name]]
  source <- read_tsv(spec_outcome$path)
  stopifnot(nrow(source) == spec_outcome$expected_n, !anyDuplicated(source$sample),
            all(is.finite(source[[spec_outcome$raw]])), all(source[[spec_outcome$raw]] > 0))
  joined <- merge(matrix, source, by = "sample", all = FALSE, suffixes = c("", "_outcome"))
  stopifnot(nrow(joined) == spec_outcome$expected_n)
  joined[[spec_outcome$transformed]] <- spec_outcome$transform(joined[[spec_outcome$raw]])

  outcome_dir <- file.path(result_root, outcome_name)
  dir.create(outcome_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(joined, file.path(outcome_dir, "person_level_outcome_matrix.tsv"))

  result_rows <- list()
  desc_rows <- list()
  keys <- character(nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    exposure <- manifest$exposure_id[i]
    adjustment_ids <- character()
    if (!is.na(manifest$adjustment[i]) && nzchar(manifest$adjustment[i])) adjustment_ids <- manifest$adjustment[i]
    x <- joined[[exposure]]
    known <- !is.na(x)
    zero <- sum(x[known] == 0)
    nonzero <- sum(x[known] != 0)
    counts <- table(x[known])
    rare <- manifest$model_scale[i] %in% c("binary", "dosage") && min(zero, nonzero) < 5L
    gate <- if (sum(known) < 30L) "underpowered_n_lt30" else if (length(counts) < 2L) "invariant" else if (rare) "underpowered_arm_lt5" else "fit"
    desc_rows[[i]] <- data.frame(
      outcome = outcome_name, exposure_id = exposure, family = manifest$family[i],
      n = sum(known), n_missing = sum(!known), n_zero = zero, n_nonzero = nonzero,
      n_unique = length(counts), arm_counts = paste(names(counts), as.integer(counts), sep = ":", collapse = ";"),
      screening_gate = gate
    )
    primary <- if (gate == "fit") cluster_fit(joined, exposure, spec_outcome$transformed, adjustment_ids, "superpopulation") else list(status = gate, n = sum(known))
    population <- if (gate == "fit") cluster_fit(joined, exposure, spec_outcome$transformed, adjustment_ids, "population") else list(status = gate, n = sum(known))
    getv <- function(object, field) if (is.null(object[[field]])) NA else object[[field]]
    result_rows[[i]] <- data.frame(
      outcome = outcome_name, outcome_units = spec_outcome$units,
      exposure_id = exposure, family = manifest$family[i], locus = manifest$locus[i],
      endpoint = manifest$endpoint[i], adjustment = paste(adjustment_ids, collapse = ";"),
      high_expression = manifest$high_expression[i], n = getv(primary, "n"),
      model_status = primary$status, beta = getv(primary, "beta"), se_pedigree_cluster = getv(primary, "se"),
      ci95_low = getv(primary, "ci_low"), ci95_high = getv(primary, "ci_high"),
      p_pedigree_cluster = getv(primary, "p"), standardized_beta = getv(primary, "standardized_beta"),
      population_model_status = population$status, population_beta = getv(population, "beta"),
      population_se_cluster = getv(population, "se"), population_p_cluster = getv(population, "p"),
      stringsAsFactors = FALSE
    )
    adjustment_vector <- if (length(adjustment_ids)) joined[[adjustment_ids[1]]] else rep(NA_real_, nrow(joined))
    keys[i] <- paste(manifest$model_scale[i],
                     paste(ifelse(is.na(x), "NA", format(x, scientific = FALSE, trim = TRUE)), collapse = ","),
                     paste(ifelse(is.na(adjustment_vector), "NA", format(adjustment_vector, scientific = FALSE, trim = TRUE)), collapse = ","), sep = "|")
  }
  results <- do.call(rbind, result_rows)
  desc <- do.call(rbind, desc_rows)
  order_index <- order(match(keys, unique(keys)), manifest$prespecified_priority, manifest$exposure_id)
  rep_flag <- logical(nrow(manifest))
  rep_flag[order_index] <- !duplicated(keys[order_index])
  results$observed_vector_group <- match(keys, unique(keys))
  results$multiplicity_representative <- rep_flag
  results$q_bh_outcome_global <- NA_real_
  results$p_holm_outcome_global <- NA_real_
  fit_ix <- results$multiplicity_representative & results$model_status == "fit" & is.finite(results$p_pedigree_cluster)
  results$q_bh_outcome_global[fit_ix] <- p.adjust(results$p_pedigree_cluster[fit_ix], method = "BH")
  results$p_holm_outcome_global[fit_ix] <- p.adjust(results$p_pedigree_cluster[fit_ix], method = "holm")
  for (g in unique(results$observed_vector_group)) {
    rows <- which(results$observed_vector_group == g)
    rep_row <- rows[results$multiplicity_representative[rows]][1]
    results$q_bh_outcome_global[rows] <- results$q_bh_outcome_global[rep_row]
    results$p_holm_outcome_global[rows] <- results$p_holm_outcome_global[rep_row]
  }
  write_tsv(results[order(results$q_bh_outcome_global, results$p_pedigree_cluster, na.last = TRUE), ],
            file.path(outcome_dir, "model_results.tsv"))
  write_tsv(desc, file.path(outcome_dir, "exposure_descriptives.tsv"))
  all_models[[outcome_name]] <- results
  all_descriptives[[outcome_name]] <- desc
}

combined <- do.call(rbind, all_models)
combined$q_bh_secondary_suite <- NA_real_
combined$p_holm_secondary_suite <- NA_real_
ix <- combined$multiplicity_representative & combined$model_status == "fit" & is.finite(combined$p_pedigree_cluster)
combined$q_bh_secondary_suite[ix] <- p.adjust(combined$p_pedigree_cluster[ix], method = "BH")
combined$p_holm_secondary_suite[ix] <- p.adjust(combined$p_pedigree_cluster[ix], method = "holm")
write_tsv(combined[order(combined$q_bh_secondary_suite, combined$p_pedigree_cluster, na.last = TRUE), ],
          file.path(result_root, "secondary_outcome_model_results.tsv"))

# A second correction spans the complete three-outcome screen, including the
# Mandage discovery phenotype, so outcome-specific and suite-wide evidence are
# never conflated.
mandage <- read_tsv(file.path(result_root, "model_results.tsv"))
mandage_keep <- mandage[mandage$multiplicity_representative & mandage$model_status == "fit" &
                          is.finite(mandage$p_pedigree_cluster), ]
all_three <- rbind(
  data.frame(outcome = "Mandage2017_EBV_in_silico", exposure_id = mandage_keep$exposure_id,
             n = mandage_keep$n, beta = mandage_keep$beta_log2_ebv,
             p_pedigree_cluster = mandage_keep$p_pedigree_cluster,
             population_p_cluster = mandage_keep$population_p_cluster),
  combined[ix, c("outcome", "exposure_id", "n", "beta", "p_pedigree_cluster", "population_p_cluster")]
)
all_three$q_bh_three_outcome_suite <- p.adjust(all_three$p_pedigree_cluster, method = "BH")
all_three$p_holm_three_outcome_suite <- p.adjust(all_three$p_pedigree_cluster, method = "holm")
write_tsv(all_three[order(all_three$q_bh_three_outcome_suite, all_three$p_pedigree_cluster), ],
          file.path(result_root, "all_three_outcome_multiplicity.tsv"))

summary_rows <- lapply(names(all_models), function(name) {
  x <- all_models[[name]]
  data.frame(outcome = name, people = outcomes[[name]]$expected_n,
             exposures = nrow(x), unique_vectors = length(unique(x$observed_vector_group)),
             fitted_representatives = sum(x$model_status == "fit" & x$multiplicity_representative),
             nominal_p_lt_0.05 = sum(x$p_pedigree_cluster < .05 & x$multiplicity_representative, na.rm = TRUE),
             outcome_bh_lt_0.05 = sum(x$q_bh_outcome_global < .05, na.rm = TRUE))
})
summary <- do.call(rbind, summary_rows)
write_tsv(summary, file.path(result_root, "secondary_outcome_summary.tsv"))
print(summary)
cat("Top secondary-outcome results:\n")
print(head(combined[ix, c("outcome", "exposure_id", "n", "beta", "p_pedigree_cluster", "q_bh_secondary_suite")][order(combined$p_pedigree_cluster[ix]), ], 15))
