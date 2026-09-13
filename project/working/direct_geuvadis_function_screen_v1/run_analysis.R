#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(".", mustWork = TRUE)
work <- file.path(root, "project/working/direct_geuvadis_function_screen_v1")
data_dir <- file.path(work, "data")
outdir <- file.path(work, "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, sep = "\t", quote = "", check.names = FALSE,
                                      na.strings = c("NA", ""), comment.char = "")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
                                           quote = FALSE, na = "NA")

expr <- read_tsv(file.path(data_dir, "direct33_expression.tsv"))
exposure <- read_tsv(file.path(data_dir, "direct33_exposures.tsv"))
manifest <- read_tsv(file.path(data_dir, "feature_manifest.tsv"))
coordinates <- read_tsv(file.path(data_dir, "locus_coordinates_grch38.tsv"))
type1_locus <- read_tsv(file.path(data_dir, "type1_locus_units.tsv"))
mage_overlap <- read_tsv(file.path(data_dir, "mage_overlap_audit.tsv"))

stopifnot(nrow(expr) == 23722L, nrow(exposure) == 33L, nrow(manifest) == 79L,
          identical(exposure$sample, type1_locus$sample), identical(exposure$sample, mage_overlap$sample),
          !anyDuplicated(exposure$sample), sum(mage_overlap$in_prior_MAGE_direct39) == 5L)
sample_columns <- exposure$sample
stopifnot(all(sample_columns %in% names(expr)))
Y <- t(as.matrix(expr[, sample_columns, drop = FALSE]))
storage.mode(Y) <- "numeric"
stopifnot(nrow(Y) == 33L, ncol(Y) == 23722L, all(is.finite(Y)))

gene_id <- expr$gene_id_base
gene_symbol <- expr$gene_symbol
gene_chr <- expr$chromosome_grch38
gene_tss <- suppressWarnings(as.numeric(expr$tss_grch38))
coord_start <- setNames(coordinates$start, coordinates$locus)
coord_end <- setNames(coordinates$end, coordinates$locus)
coord_chr <- setNames(coordinates$chromosome, coordinates$locus)

# Collapse exact observed exposure plus adjustment vectors in this cohort. This
# avoids counting multiple biological labels for the same numerical test twice.
vector_key <- character(nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  id <- manifest$exposure_id[i]
  adj <- manifest$adjustment[i]
  x_key <- paste(ifelse(is.na(exposure[[id]]), "NA", exposure[[id]]), collapse = ",")
  adj_key <- if (is.na(adj) || !nzchar(adj)) "" else {
    paste(ifelse(is.na(exposure[[adj]]), "NA", exposure[[adj]]), collapse = ",")
  }
  vector_key[i] <- paste(manifest$model_scale[i], x_key, adj_key, sep = "|")
}
manifest$observed_direct33_vector_group <- match(vector_key, unique(vector_key))
representative_order <- order(manifest$observed_direct33_vector_group,
                              manifest$prespecified_priority, manifest$exposure_id)
manifest$direct33_multiplicity_representative <- FALSE
manifest$direct33_multiplicity_representative[representative_order] <-
  !duplicated(manifest$observed_direct33_vector_group[representative_order])

exposure_gate <- function(i) {
  id <- manifest$exposure_id[i]
  x <- exposure[[id]]
  adj <- manifest$adjustment[i]
  needed <- is.finite(x)
  if (!is.na(adj) && nzchar(adj)) needed <- needed & is.finite(exposure[[adj]])
  x_known <- x[needed]
  if (length(x_known) < 20L) return("underpowered_n_lt20")
  if (length(unique(x_known)) < 2L) return("invariant")
  if (manifest$model_scale[i] %in% c("binary", "dosage")) {
    arm <- c(sum(x_known == 0), sum(x_known != 0))
    if (min(arm) < 5L) return("underpowered_arm_lt5")
  } else if (length(unique(x_known)) < 3L) {
    return("insufficient_continuous_variation")
  }
  "fit"
}
manifest$screening_gate <- vapply(seq_len(nrow(manifest)), exposure_gate, character(1))

genes_for_group <- function(group_rows) {
  loci <- unique(manifest$locus[group_rows])
  if ("catalog-wide" %in% loci) {
    return(list(rows = seq_len(nrow(expr)), lane = "general_transcriptome",
                loci = "catalog-wide", missing = ""))
  }
  known_loci <- intersect(loci, coordinates$locus)
  missing <- setdiff(loci, known_loci)
  rows <- integer()
  for (locus in known_loci) {
    rows <- union(rows, which(gene_chr == coord_chr[locus] & is.finite(gene_tss) &
                              gene_tss >= coord_start[locus] - 1000000 &
                              gene_tss <= coord_end[locus] + 1000000))
  }
  list(rows = sort(rows), lane = "locus_cis_1Mb", loci = paste(sort(loci), collapse = ";"),
       missing = paste(sort(missing), collapse = ";"))
}

fit_matrix_hc3 <- function(y, x, adjustment = NULL, ancestry = c("population", "superpopulation")) {
  ancestry <- match.arg(ancestry)
  keep <- is.finite(x)
  if (!is.null(adjustment)) keep <- keep & is.finite(adjustment)
  y <- y[keep, , drop = FALSE]
  x <- x[keep]
  d <- data.frame(x = x, sex = factor(exposure$sex[keep]),
                  ancestry = factor(exposure[[ancestry]][keep]))
  if (!is.null(adjustment)) d$adjustment <- adjustment[keep]
  # Reduce redundant nuisance columns first (for example, an invariant dosage
  # adjustment), then require the exposure to add rank. This preserves an
  # estimable ORF contrast but still blocks a functional count identical to its
  # intended copy-number adjustment.
  nuisance_formula <- if (is.null(adjustment)) ~ ancestry + sex else ~ adjustment + ancestry + sex
  nuisance <- model.matrix(nuisance_formula, data = d)
  nuisance_qr <- qr(nuisance)
  nuisance <- nuisance[, nuisance_qr$pivot[seq_len(nuisance_qr$rank)], drop = FALSE]
  X <- cbind(nuisance, x = x)
  qr_x <- qr(X)
  if (qr_x$rank < ncol(X)) {
    return(list(status = "exposure_aliased", n = nrow(X)))
  }
  inverse <- solve(crossprod(X))
  coefficients <- inverse %*% crossprod(X, y)
  residuals <- y - X %*% coefficients
  leverage <- rowSums((X %*% inverse) * X)
  coefficient_weights <- (inverse %*% t(X))[match("x", colnames(X)), ]
  adjusted_residuals <- sweep(residuals, 1, 1 - leverage, "/")
  score <- sweep(adjusted_residuals, 1, coefficient_weights, "*")
  se <- sqrt(colSums(score^2))
  beta <- as.numeric(coefficients[match("x", colnames(X)), ])
  statistic <- beta / se
  df <- nrow(X) - qr_x$rank
  p <- 2 * pt(abs(statistic), df = df, lower.tail = FALSE)
  y_sd <- apply(y, 2, sd)
  list(status = "fit", n = nrow(X), beta = beta, se = se, statistic = statistic,
       p = p, df = df, ci_low = beta + qt(0.025, df) * se,
       ci_high = beta + qt(0.975, df) * se,
       standardized_beta = beta * sd(x) / y_sd)
}

result_list <- list()
sensitivity_only_list <- list()
status_list <- list()
groups <- sort(unique(manifest$observed_direct33_vector_group))
for (group in groups) {
  group_rows <- which(manifest$observed_direct33_vector_group == group)
  representative <- group_rows[manifest$direct33_multiplicity_representative[group_rows]][1]
  id <- manifest$exposure_id[representative]
  routing <- genes_for_group(group_rows)
  gate <- manifest$screening_gate[representative]
  if (gate == "fit" && !length(routing$rows)) gate <- "no_mapped_cis_gene"
  status_list[[length(status_list) + 1L]] <- data.frame(
    observed_direct33_vector_group = group, representative_exposure_id = id,
    alias_count = length(group_rows), aliases = paste(sort(manifest$exposure_id[group_rows]), collapse = ";"),
    family = manifest$family[representative], lane = routing$lane, tested_loci = routing$loci,
    missing_coordinate_loci = routing$missing, selected_gene_count = length(routing$rows),
    screening_gate = gate, stringsAsFactors = FALSE
  )
  if (gate != "fit") next
  adjustment_id <- manifest$adjustment[representative]
  adjustment <- if (is.na(adjustment_id) || !nzchar(adjustment_id)) NULL else exposure[[adjustment_id]]
  x <- exposure[[id]]
  primary <- fit_matrix_hc3(Y[, routing$rows, drop = FALSE], x, adjustment, "population")
  sensitivity <- fit_matrix_hc3(Y[, routing$rows, drop = FALSE], x, adjustment, "superpopulation")
  if (primary$status != "fit") {
    status_list[[length(status_list)]][["screening_gate"]] <- if (sensitivity$status == "fit") {
      paste0(primary$status, "_superpopulation_sensitivity_fit")
    } else primary$status
    if (sensitivity$status == "fit") {
      sensitivity_only_list[[length(sensitivity_only_list) + 1L]] <- data.frame(
        observed_direct33_vector_group = group, exposure_id = id, alias_count = length(group_rows),
        family = manifest$family[representative], lane = routing$lane,
        locus = manifest$locus[representative], endpoint = manifest$endpoint[representative],
        adjustment = ifelse(is.null(adjustment), "", adjustment_id),
        high_expression = manifest$high_expression[representative],
        gene_id = gene_id[routing$rows], gene_symbol = gene_symbol[routing$rows],
        chromosome_grch38 = gene_chr[routing$rows], tss_grch38 = gene_tss[routing$rows],
        n = sensitivity$n, beta_superpopulation_sex = sensitivity$beta,
        se_hc3 = sensitivity$se, ci95_low = sensitivity$ci_low,
        ci95_high = sensitivity$ci_high, standardized_beta = sensitivity$standardized_beta,
        statistic_hc3 = sensitivity$statistic,
        p_superpopulation_sex_hc3 = sensitivity$p, residual_df = sensitivity$df,
        stringsAsFactors = FALSE
      )
    }
    next
  }
  if (sensitivity$status != "fit") {
    sensitivity <- list(beta = rep(NA_real_, length(routing$rows)), p = rep(NA_real_, length(routing$rows)))
  }
  locus_for_distance <- unique(manifest$locus[group_rows])
  midpoint <- if (length(locus_for_distance) == 1L && locus_for_distance %in% coordinates$locus) {
    mean(c(coord_start[locus_for_distance], coord_end[locus_for_distance]))
  } else NA_real_
  result_list[[length(result_list) + 1L]] <- data.frame(
    observed_direct33_vector_group = group, exposure_id = id,
    alias_count = length(group_rows), family = manifest$family[representative],
    lane = routing$lane, locus = manifest$locus[representative], endpoint = manifest$endpoint[representative],
    adjustment = ifelse(is.null(adjustment), "", adjustment_id),
    high_expression = manifest$high_expression[representative],
    gene_id = gene_id[routing$rows], gene_symbol = gene_symbol[routing$rows],
    chromosome_grch38 = gene_chr[routing$rows], tss_grch38 = gene_tss[routing$rows],
    distance_to_locus_midpoint = gene_tss[routing$rows] - midpoint,
    n = primary$n, beta_peer50_rpkm = primary$beta, se_hc3 = primary$se,
    ci95_low = primary$ci_low, ci95_high = primary$ci_high,
    standardized_beta = primary$standardized_beta, statistic_hc3 = primary$statistic,
    p_population_sex_hc3 = primary$p, residual_df = primary$df,
    beta_superpopulation_sex = sensitivity$beta,
    p_superpopulation_sex_hc3 = sensitivity$p,
    stringsAsFactors = FALSE
  )
}

status <- do.call(rbind, status_list)
results <- do.call(rbind, result_list)
stopifnot(!is.null(results), nrow(results) > 0L)
finite <- is.finite(results$p_population_sex_hc3)
results$q_bh_global <- NA_real_
results$p_holm_global <- NA_real_
results$q_bh_global[finite] <- p.adjust(results$p_population_sex_hc3[finite], method = "BH")
results$p_holm_global[finite] <- p.adjust(results$p_population_sex_hc3[finite], method = "holm")
results$q_bh_lane <- ave(results$p_population_sex_hc3, results$lane,
                         FUN = function(p) p.adjust(p, method = "BH"))
results$q_bh_within_exposure <- ave(results$p_population_sex_hc3, results$exposure_id,
                                    FUN = function(p) p.adjust(p, method = "BH"))
results <- results[order(results$p_population_sex_hc3), ]

write_tsv(manifest, file.path(outdir, "exposure_alias_map.tsv"))
write_tsv(status, file.path(outdir, "screen_status.tsv"))
write_tsv(results, file.path(outdir, "broad_model_results.tsv"))
if (length(sensitivity_only_list)) {
  sensitivity_only <- do.call(rbind, sensitivity_only_list)
  sensitivity_only$q_bh_superpopulation_sensitivity <-
    p.adjust(sensitivity_only$p_superpopulation_sex_hc3, method = "BH")
  sensitivity_only <- sensitivity_only[order(sensitivity_only$p_superpopulation_sex_hc3), ]
} else {
  sensitivity_only <- data.frame(
    observed_direct33_vector_group = integer(), exposure_id = character(), alias_count = integer(),
    family = character(), lane = character(), locus = character(), endpoint = character(),
    adjustment = character(), high_expression = logical(), gene_id = character(), gene_symbol = character(),
    chromosome_grch38 = character(), tss_grch38 = numeric(), n = integer(),
    beta_superpopulation_sex = numeric(), se_hc3 = numeric(), ci95_low = numeric(),
    ci95_high = numeric(), standardized_beta = numeric(), statistic_hc3 = numeric(),
    p_superpopulation_sex_hc3 = numeric(), residual_df = numeric(),
    q_bh_superpopulation_sensitivity = numeric(), stringsAsFactors = FALSE
  )
}
write_tsv(sensitivity_only, file.path(outdir, "superpopulation_only_sensitivity.tsv"))

# Targeted SLC44A5 replication and Type-I-shared-versus-local heterogeneity.
slc_row <- which(expr$gene_id_base == "ENSG00000137968")
stopifnot(length(slc_row) == 1L, expr$gene_symbol[slc_row] == "SLC44A5")
slc_y <- Y[, slc_row]

fit_scalar_hc3 <- function(y, x, adjustment = NULL,
                           ancestry = c("population", "superpopulation")) {
  fit <- fit_matrix_hc3(matrix(y, ncol = 1L), x, adjustment, match.arg(ancestry))
  if (fit$status != "fit") {
    return(data.frame(status = fit$status, n = fit$n, beta = NA_real_, se_hc3 = NA_real_,
                      ci95_low = NA_real_, ci95_high = NA_real_, standardized_beta = NA_real_,
                      statistic_hc3 = NA_real_, p_two_sided_hc3 = NA_real_, residual_df = NA_real_))
  }
  data.frame(status = "fit", n = fit$n, beta = fit$beta, se_hc3 = fit$se,
             ci95_low = fit$ci_low, ci95_high = fit$ci_high,
             standardized_beta = fit$standardized_beta, statistic_hc3 = fit$statistic,
             p_two_sided_hc3 = fit$p, residual_df = fit$df)
}

target_rows <- list()
add_target <- function(model_id, focal_id, x, adjustment_id = "", adjustment = NULL,
                       family = "targeted", model_scale = "continuous") {
  known <- is.finite(x)
  if (!is.null(adjustment)) known <- known & is.finite(adjustment)
  x_known <- x[known]
  gate <- if (length(x_known) < 20L) "underpowered_n_lt20" else if (length(unique(x_known)) < 2L) {
    "invariant"
  } else if (model_scale %in% c("binary", "dosage") &&
             min(table(x_known)) < 5L) {
    "underpowered_arm_lt5"
  } else "fit"
  for (ancestry in c("population", "superpopulation")) {
    fit <- if (gate == "fit") fit_scalar_hc3(slc_y, x, adjustment, ancestry) else {
      data.frame(status = gate, n = sum(known), beta = NA_real_, se_hc3 = NA_real_,
                 ci95_low = NA_real_, ci95_high = NA_real_, standardized_beta = NA_real_,
                 statistic_hc3 = NA_real_, p_two_sided_hc3 = NA_real_, residual_df = NA_real_)
    }
    fit$model_id <- model_id
    fit$focal_exposure_id <- focal_id
    fit$adjustment <- adjustment_id
    fit$family <- family
    fit$ancestry_model <- ancestry
    target_rows[[length(target_rows) + 1L]] <<- fit[, c("model_id", "focal_exposure_id", "adjustment",
                                                        "family", "ancestry_model", setdiff(names(fit),
                                                        c("model_id", "focal_exposure_id", "adjustment",
                                                          "family", "ancestry_model")))]
  }
}

add_subset_target <- function(model_id, focal_id, x, keep, family) {
  # Force the focal vector before replacing the global exposure frame. R's
  # lazy argument evaluation would otherwise re-evaluate exposure[[focal_id]]
  # against the already-subset frame and then apply `keep` a second time.
  x <- as.numeric(x)
  old_exposure <- exposure
  exposure <<- exposure[keep, , drop = FALSE]
  for (ancestry in c("population", "superpopulation")) {
    fit <- fit_scalar_hc3(slc_y[keep], x[keep], NULL, ancestry)
    fit$model_id <- model_id
    fit$focal_exposure_id <- focal_id
    fit$adjustment <- ""
    fit$family <- family
    fit$ancestry_model <- ancestry
    target_rows[[length(target_rows) + 1L]] <<- fit[, c("model_id", "focal_exposure_id", "adjustment",
                                                        "family", "ancestry_model", setdiff(names(fit),
                                                        c("model_id", "focal_exposure_id", "adjustment",
                                                          "family", "ancestry_model")))]
  }
  exposure <<- old_exposure
}

local_id <- "HML-2_1p31.1b::internal_fragment_present"
type1_id <- "general::type1_physical_units"
add_target("prespecified_MAGE_replication", local_id, exposure[[local_id]],
           family = "prespecified_replication", model_scale = "binary")
add_subset_target("MAGE_nonoverlap_replication", local_id, exposure[[local_id]],
                  !mage_overlap$in_prior_MAGE_direct39, "independent_donor_sensitivity")
add_target("general_type1_physical_burden", type1_id, exposure[[type1_id]], family = "shared_type1")
add_target("local_state_adjusted_type1_burden", local_id, exposure[[local_id]], type1_id,
           exposure[[type1_id]], "local_vs_shared_type1")
add_target("type1_burden_adjusted_local_state", type1_id, exposure[[type1_id]], local_id,
           exposure[[local_id]], "local_vs_shared_type1")

onep_ids <- c("HML-2_1p31.1b::total_internal_fragment_cn", "HML-2_1p31.1b::any_array",
              "HML-2_1p31.1b::max_haplotype_cn", "HML-2_1p31.1b::compatible_env_units",
              "HML-2_1p31.1b::compatible_np9_units")
for (id in onep_ids) {
  row <- match(id, manifest$exposure_id)
  adj_id <- manifest$adjustment[row]
  adj <- if (is.na(adj_id) || !nzchar(adj_id)) NULL else exposure[[adj_id]]
  add_target(paste0("onep31b_feature::", id), id, exposure[[id]],
             ifelse(is.null(adj), "", adj_id), adj, "onep31b_secondary",
             manifest$model_scale[row])
}
for (id in c("general::type1_loci_present", "general::type1_compatible_gag_units",
             "general::type1_compatible_env_units", "general::type1_compatible_np9_units",
             "general::type2_physical_units")) {
  row <- match(id, manifest$exposure_id)
  adj_id <- manifest$adjustment[row]
  adj <- if (is.na(adj_id) || !nzchar(adj_id)) NULL else exposure[[adj_id]]
  add_target(paste0("general_feature::", id), id, exposure[[id]],
             ifelse(is.null(adj), "", adj_id), adj, "general_secondary",
             manifest$model_scale[row])
}

type1_names <- setdiff(names(type1_locus), "sample")
type1_matrix <- as.matrix(type1_locus[, type1_names, drop = FALSE])
storage.mode(type1_matrix) <- "numeric"
stopifnot(all(rowSums(type1_matrix) == exposure[[type1_id]]))
for (locus in type1_names) {
  add_target(paste0("type1_locus_units::", locus), locus, type1_locus[[locus]],
             family = "type1_locus_heterogeneity", model_scale = "dosage")
  leave_out <- rowSums(type1_matrix) - type1_locus[[locus]]
  add_target(paste0("type1_leave_one_locus_out::", locus), paste0("type1_units_without_", locus),
             leave_out, family = "type1_burden_leave_one_locus_out")
}
targeted <- do.call(rbind, target_rows)
targeted$p_one_sided_positive <- ifelse(targeted$status == "fit",
                                         pt(targeted$statistic_hc3, df = targeted$residual_df,
                                            lower.tail = FALSE), NA_real_)
targeted$q_bh_secondary_family <- NA_real_
secondary_ix <- targeted$status == "fit" & targeted$ancestry_model == "population" &
  targeted$family != "prespecified_replication"
targeted$q_bh_secondary_family[secondary_ix] <- p.adjust(targeted$p_two_sided_hc3[secondary_ix], "BH")
write_tsv(targeted, file.path(outdir, "slc44a5_targeted_results.tsv"))

# Population/exposure summaries and leave-one-person-out influence for the exact
# prespecified local-state replication.
group_summary <- aggregate(slc_y,
                           by = list(internal_fragment_present = exposure[[local_id]],
                                     population = exposure$population),
                           FUN = function(z) c(n = length(z), mean = mean(z), sd = if (length(z) > 1) sd(z) else NA))
group_summary <- data.frame(internal_fragment_present = group_summary$internal_fragment_present,
                            population = group_summary$population, n = group_summary$x[, "n"],
                            mean_peer50_rpkm = group_summary$x[, "mean"], sd_peer50_rpkm = group_summary$x[, "sd"])
write_tsv(group_summary, file.path(outdir, "slc44a5_group_summary.tsv"))

loo <- lapply(seq_len(nrow(exposure)), function(drop) {
  keep <- seq_len(nrow(exposure)) != drop
  old_exposure <- exposure
  exposure <<- exposure[keep, , drop = FALSE]
  pop <- fit_scalar_hc3(slc_y[keep], old_exposure[[local_id]][keep], NULL, "population")
  super <- fit_scalar_hc3(slc_y[keep], old_exposure[[local_id]][keep], NULL, "superpopulation")
  exposure <<- old_exposure
  data.frame(dropped_sample = old_exposure$sample[drop],
             dropped_state = old_exposure[[local_id]][drop], dropped_outcome = slc_y[drop],
             beta_population = ifelse(pop$status == "fit", pop$beta, NA),
             p_population_hc3 = ifelse(pop$status == "fit", pop$p_two_sided_hc3, NA),
             beta_superpopulation = ifelse(super$status == "fit", super$beta, NA),
             p_superpopulation_hc3 = ifelse(super$status == "fit", super$p_two_sided_hc3, NA))
})
loo <- do.call(rbind, loo)
write_tsv(loo, file.path(outdir, "slc44a5_leave_one_person_out.tsv"))

# Preserve carrier counts within each population while shuffling the local
# state. This addresses the most immediate population-composition alternative
# without treating the small cohort as if ancestry were continuous.
set.seed(20260716)
permutations <- 50000L
observed_permutation_fit <- fit_scalar_hc3(slc_y, exposure[[local_id]], NULL, "population")
permuted_statistics <- numeric(permutations)
population_rows <- split(seq_len(nrow(exposure)), exposure$population)
for (b in seq_len(permutations)) {
  permuted <- exposure[[local_id]]
  for (rows in population_rows) permuted[rows] <- sample(permuted[rows], replace = FALSE)
  permuted_fit <- fit_scalar_hc3(slc_y, permuted, NULL, "population")
  permuted_statistics[b] <- permuted_fit$statistic_hc3
}
permutation_result <- data.frame(
  permutations = permutations,
  observed_statistic_hc3 = observed_permutation_fit$statistic_hc3,
  empirical_p_two_sided = (1 + sum(abs(permuted_statistics) >= abs(observed_permutation_fit$statistic_hc3))) /
    (permutations + 1),
  empirical_p_one_sided_positive = (1 + sum(permuted_statistics >= observed_permutation_fit$statistic_hc3)) /
    (permutations + 1),
  permutation_blocks = "population; within-block carrier counts preserved",
  seed = 20260716L,
  stringsAsFactors = FALSE
)
write_tsv(permutation_result, file.path(outdir, "slc44a5_population_blocked_permutation.tsv"))

# Rank-inverse-normal sensitivity for the primary replication only.
slc_inverse_normal <- qnorm((rank(slc_y, ties.method = "average") - 0.5) / length(slc_y))
inverse_rows <- do.call(rbind, lapply(c("population", "superpopulation"), function(ancestry) {
  fit <- fit_scalar_hc3(slc_inverse_normal, exposure[[local_id]], NULL, ancestry)
  fit$ancestry_model <- ancestry
  fit
}))
write_tsv(inverse_rows, file.path(outdir, "slc44a5_inverse_normal_sensitivity.tsv"))

priority_pattern <- "^(general::|HML-2_1p31.1b::|HML-2_1q22::|HML-2_7p22.1::)"
focal_cis <- results[grepl("^HML-2_(1p31.1b|1q22|7p22.1)::", results$exposure_id), ]
general_or_high <- results[grepl("^general::", results$exposure_id) | results$high_expression, ]
priority <- unique(rbind(focal_cis, head(general_or_high, 500L),
                         results[results$gene_id == "ENSG00000137968", ]))
priority <- priority[order(priority$p_population_sex_hc3), ]
write_tsv(priority, file.path(outdir, "priority_results.tsv"))

summary <- data.frame(
  metric = c("official_expression_people", "exact_direct_people", "genes", "gencode_mapped_genes",
             "prior_MAGE_overlap_people", "GEUVADIS_only_people", "frozen_exposures",
             "unique_direct33_vectors", "fitted_vector_groups",
             "broad_unique_tests", "broad_nominal_p_lt_0.05", "broad_global_bh_lt_0.05",
             "broad_global_holm_lt_0.05"),
  value = c(462, nrow(exposure), nrow(expr), sum(is.finite(gene_tss)),
            sum(mage_overlap$in_prior_MAGE_direct39), sum(!mage_overlap$in_prior_MAGE_direct39),
            nrow(manifest), length(groups),
            sum(status$screening_gate == "fit"), nrow(results),
            sum(results$p_population_sex_hc3 < 0.05, na.rm = TRUE),
            sum(results$q_bh_global < 0.05, na.rm = TRUE),
            sum(results$p_holm_global < 0.05, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write_tsv(summary, file.path(outdir, "run_summary.tsv"))

cat("GEUVADIS direct expression screen complete\n")
print(summary)
cat("Prespecified SLC44A5 replication:\n")
print(targeted[targeted$model_id == "prespecified_MAGE_replication",
               c("ancestry_model", "n", "beta", "se_hc3", "p_two_sided_hc3", "p_one_sided_positive")])
cat("Top broad results:\n")
print(head(results[, c("exposure_id", "gene_symbol", "gene_id", "n", "standardized_beta",
                       "p_population_sex_hc3", "q_bh_global", "p_superpopulation_sex_hc3")], 15L))
