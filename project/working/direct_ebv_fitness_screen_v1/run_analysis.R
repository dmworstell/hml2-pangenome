options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(".")
outdir <- file.path(root, "project/working/direct_ebv_fitness_screen_v1/results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, sep = "\t", quote = "", check.names = FALSE,
                                       na.strings = c("NA"), comment.char = "")
write_tsv <- function(x, name) write.table(x, file.path(outdir, name), sep = "\t",
                                            row.names = FALSE, quote = FALSE, na = "NA")
compatible <- c("Intact", "Frameshift_at_end", "Intact_FS_End")
candidate_status <- c("eligible_primary_multi_ancestry", "eligible_primary_ancestry_restricted",
                      "eligible_secondary_structural", "eligible_exploratory",
                      "eligible_secondary_exploratory", "targeted_positive_control")

paths <- list(
  phenotype = "project/working/direct_ebv_fitness_screen_v1/data/mandage_ebv_load.tsv",
  truth = "project/inputs/phase0_functional_screen_v1/sample_contrast_truth.tsv",
  universe = "project/inputs/phase0_functional_screen_v1/contrast_universe.tsv",
  catalog = "HML2_ProjectResources/data/catalog/combined_hml2_orf_analysis.tsv",
  arrays = "project/working/sevenp22_tandem_array_agent/haplotype_copy_orf_reconciliation.tsv",
  arrays_7p22_cn_truth = "project/working/sevenp22_proxy_resolution_agent/haplotype_copy_number_truth.tsv",
  igsr = "project/inputs/phase0_functional_screen_v1/upstream_igsr_samples.tsv",
  high_expression = "project/working/functional_state_population_agent/results/priority_high_expression_screen.tsv"
)
stopifnot(all(file.exists(unlist(paths))))

phenotype <- read_tsv(paths$phenotype)
stopifnot(nrow(phenotype) == 1753L, !anyDuplicated(phenotype$samples),
          all(is.finite(phenotype$`EBV load`)), all(phenotype$`EBV load` > 0))
names(phenotype)[names(phenotype) == "samples"] <- "sample"
names(phenotype)[names(phenotype) == "EBV load"] <- "ebv_load"

arrays <- read_tsv(paths$arrays)
direct_ids <- sort(unique(arrays$sample))
stopifnot(length(direct_ids) == 292L)

# 7p22.1 copy number is taken from the authoritative corrected truth, which is
# reproduced independently by array_function_burden_v1 and array_structural_reaudit_v1.
# The reconciliation table above counts an assembly-gap duplicate lying outside the
# canonical locus (allele_class single_provirus_asmdup_flag) as a real copy; the truth
# table calls it CN0.
#
# This overrides rather than repoints. The truth table is 7p22.1-only, carries no
# locus column and no ORF unit columns, so it cannot serve aggregate_array(), which
# also aggregates HML-2_1p31.1b off arrays$locus/copy_number and the compatible_*_units
# columns. Repointing paths$arrays at it would silently drop 1p31.1b from the screen.
arrays_cn_truth <- read_tsv(paths$arrays_7p22_cn_truth)
stopifnot(nrow(arrays_cn_truth) == 584L,
          !anyDuplicated(paste(arrays_cn_truth$sample, arrays_cn_truth$haplotype)))

is_7p22 <- arrays$locus == "HML-2_7p22.1"
stopifnot(sum(is_7p22) == 584L)
key_recon <- paste(arrays$sample[is_7p22], arrays$haplotype[is_7p22])
key_truth <- paste(arrays_cn_truth$sample, arrays_cn_truth$haplotype)
stopifnot(setequal(key_recon, key_truth))

cn_old <- arrays$copy_number[is_7p22]
cn_new <- arrays_cn_truth$array_copy_number[match(key_recon, key_truth)]
diverged <- which(cn_old != cn_new)

# Fail loudly on any divergence that is not an assembly-duplicate flag. A silent CN
# change here propagates into every 7p22.1 exposure and into person_level_direct_matrix.tsv,
# which direct_geuvadis_function_screen_v1 consumes.
stopifnot(all(grepl("asmdup", arrays$allele_class[is_7p22][diverged])))
# A corrected haplotype must contribute no compatible ORF units, or total_cn and the
# unit sums would disagree about whether the copy exists.
stopifnot(all(arrays$fully_compatible_units[is_7p22][diverged] == 0L),
          all(arrays$compatible_gag_units[is_7p22][diverged] == 0L))
message(sprintf("7p22.1 CN reconciliation: %d haplotype(s) corrected (%s)",
                length(diverged),
                paste(key_recon[diverged], cn_old[diverged], "->", cn_new[diverged],
                      collapse = "; ")))
arrays$copy_number[is_7p22] <- cn_new
# Corrected 7p22.1 totals: 852 copies over 584 haplotypes, CN distribution
# 3/334/230/12/4/1 for CN 0/1/2/3/4/6. The uncorrected table gives 853 and 2/335/...
stopifnot(sum(arrays$copy_number[is_7p22]) == 852L,
          identical(as.integer(table(factor(arrays$copy_number[is_7p22],
                                            levels = c(0, 1, 2, 3, 4, 6)))),
                    c(3L, 334L, 230L, 12L, 4L, 1L)))

meta <- read_tsv(paths$igsr)
meta <- meta[!duplicated(meta$`Sample name`), c("Sample name", "Sex", "Population code", "Superpopulation code")]
names(meta) <- c("sample", "sex", "population_meta", "superpopulation")

truth <- read_tsv(paths$truth)
universe <- read_tsv(paths$universe)
ped <- unique(truth[, c("sample", "pedigree_component")])
ped <- ped[!duplicated(ped$sample), ]

base <- merge(data.frame(sample = direct_ids), phenotype, by = "sample", all.x = TRUE)
base <- merge(base, meta, by = "sample", all.x = TRUE)
base <- merge(base, ped, by = "sample", all.x = TRUE)
base$pedigree_component[is.na(base$pedigree_component)] <- paste0("SINGLETON::", base$sample[is.na(base$pedigree_component)])
base$population <- ifelse(is.na(base$population_meta), base$pop, base$population_meta)
base$log2_ebv_load <- log2(base$ebv_load)
direct_overlap <- base[!is.na(base$ebv_load), ]
stopifnot(nrow(direct_overlap) == 116L, !any(is.na(direct_overlap$sex)),
          !any(is.na(direct_overlap$superpopulation)))

same_component <- table(direct_overlap$pedigree_component)
related_component_count <- sum(same_component > 1L)

catalog <- read_tsv(paths$catalog)
catalog <- catalog[grepl("^(HG|NA)[0-9]+$", catalog$ID) & catalog$ID %in% direct_ids, ]
stopifnot(length(unique(catalog$ID)) == 292L)

# Every locus has a fixed Type-I or fixed Type-II identity in the direct catalog;
# no row at a locus contributes zero physical units, following the project rule
# that a missing insertion call denotes absence rather than an unknown genotype.
prov <- catalog[catalog$Structure %in% c("Provirus", "Provirus_from_Multi") &
                  catalog$provirus_type %in% c("type1", "type2"), ]
locus_types <- aggregate(provirus_type ~ Locus, prov, function(z) paste(sort(unique(z)), collapse = ";"))
stopifnot(all(locus_types$provirus_type %in% c("type1", "type2")))
type_map <- setNames(locus_types$provirus_type, locus_types$Locus)

unit_table <- as.data.frame.matrix(xtabs(~ ID + Locus, prov))
for (sid in setdiff(direct_ids, rownames(unit_table))) unit_table[sid, ] <- 0
unit_table <- unit_table[direct_ids, , drop = FALSE]
type1_loci <- names(type_map)[type_map == "type1"]
type2_loci <- names(type_map)[type_map == "type2"]
stopifnot(length(type1_loci) == 17L, length(type2_loci) == 37L)

type_authority_audit <- data.frame(
  locus = names(type_map), observed_types = unname(type_map),
  type_class = ifelse(type_map == "type1", "fixed_type1", "fixed_type2"),
  stringsAsFactors = FALSE
)
write_tsv(type_authority_audit[order(type_authority_audit$locus), ], "type_authority_audit.tsv")

hap_labels <- unique(arrays[, c("sample", "haplotype")])
expected_cells <- length(unique(catalog$Locus)) * nrow(hap_labels)
observed_cells <- length(unique(paste(catalog$Locus, catalog$ID, catalog$Haplotype, sep = "|")))
write_tsv(data.frame(
  rule = c("expected_locus_haplotype_cells", "observed_catalog_cells", "no_row_cells",
           "catalog_wide_burden_no_row_value", "internal_orf_no_row_value"),
  value = c(expected_cells, observed_cells, expected_cells - observed_cells, 0, "not_applicable"),
  interpretation = c("98 loci x 584 direct haplotypes", "At least one catalog row exists",
                     "Treated as insertion absent for structural and burden truth",
                     "Zero proviral/functional units, never missing", "Excluded from internal-function denominators")
), "absence_semantics_audit.tsv")

normalize_structure <- function(x) {
  ifelse(x == "Provirus_from_Multi", "Multi-copy",
         ifelse(x == "Insertion_Absent", "Absent",
                ifelse(x %in% c("Provirus", "Solo-LTR", "Fragment", "Absent"), x, "Other")))
}
structure_precedence <- c(Absent = 0, Other = 1, Fragment = 2, `Solo-LTR` = 3, Provirus = 4, `Multi-copy` = 5)
catalog$normalized_structure <- normalize_structure(catalog$Structure)
cell_states <- tapply(catalog$normalized_structure,
                      paste(catalog$Locus, catalog$ID, catalog$Haplotype, sep = "|"),
                      function(z) names(which.max(structure_precedence[unique(z)])))
hap_by_person <- split(hap_labels$haplotype, hap_labels$sample)

corrected_structural_dosage <- function(locus, contrast_id) {
  authority <- universe[universe$locus == locus & universe$contrast_domain == "structural" &
                          universe$contrast_id == contrast_id, ]
  stopifnot(nrow(authority) == 1L)
  focal <- strsplit(authority$focal_states, ";", fixed = TRUE)[[1]]
  reference <- strsplit(authority$reference_states, ";", fixed = TRUE)[[1]]
  result <- setNames(rep(NA_real_, length(direct_ids)), direct_ids)
  for (sid in direct_ids) {
    states <- vapply(hap_by_person[[sid]], function(hap) {
      key <- paste(locus, sid, hap, sep = "|")
      if (key %in% names(cell_states)) unname(cell_states[key]) else "Absent"
    }, character(1))
    if (all(states %in% c(focal, reference))) result[sid] <- sum(states %in% focal)
  }
  result
}

person <- data.frame(sample = direct_ids)
person$type1_units <- rowSums(unit_table[, type1_loci, drop = FALSE])
person$type2_units <- rowSums(unit_table[, type2_loci, drop = FALSE])
person$type1_loci_present <- rowSums(unit_table[, type1_loci, drop = FALSE] > 0)
person$type2_loci_present <- rowSums(unit_table[, type2_loci, drop = FALSE] > 0)
person$type1_fraction <- person$type1_units / pmax(person$type1_units + person$type2_units, 1)
person$type1_type2_log_ratio <- log2((person$type1_units + 0.5) / (person$type2_units + 0.5))
person$type1_minus_type2_units <- person$type1_units - person$type2_units

for (tp in c("type1", "type2")) {
  sub <- prov[prov$provirus_type == tp, ]
  for (feature in c("gag", "pro", "pol", "env", "np9", "rec")) {
    v <- tapply(sub[[feature]] %in% compatible, sub$ID, sum)
    person[[paste0(tp, "_compatible_", feature, "_units")]] <- as.numeric(v[person$sample])
    person[[paste0(tp, "_compatible_", feature, "_units")]][is.na(person[[paste0(tp, "_compatible_", feature, "_units")]])] <- 0
  }
}

aggregate_array <- function(locus) {
  sub <- arrays[arrays$locus == locus, ]
  split_rows <- split(sub, sub$sample)
  rows <- lapply(direct_ids, function(sid) {
    d <- split_rows[[sid]]
    data.frame(
      sample = sid,
      total_cn = sum(d$copy_number),
      max_haplotype_cn = max(d$copy_number),
      array_haplotypes = sum(d$copy_number >= 2),
      high_haplotypes_ge3 = sum(d$copy_number >= 3),
      fully_compatible_units = sum(d$fully_compatible_units),
      compatible_gag_units = sum(d$compatible_gag_units),
      compatible_pro_units = sum(d$compatible_pro_units),
      compatible_pol_units = sum(d$compatible_pol_units),
      compatible_env_units = sum(d$compatible_env_units),
      compatible_np9_units = sum(d$compatible_np9_units),
      compatible_rec_units = sum(d$compatible_rec_units)
    )
  })
  do.call(rbind, rows)
}

seven <- aggregate_array("HML-2_7p22.1")
onep <- aggregate_array("HML-2_1p31.1b")
names(seven)[-1] <- paste0("sevenp22_", names(seven)[-1])
names(onep)[-1] <- paste0("onep31b_", names(onep)[-1])
person <- merge(person, seven, by = "sample", all.x = TRUE)
person <- merge(person, onep, by = "sample", all.x = TRUE)
person$sevenp22_any_array <- as.integer(person$sevenp22_array_haplotypes > 0)
person$sevenp22_any_high_haplotype_ge3 <- as.integer(person$sevenp22_high_haplotypes_ge3 > 0)
person$sevenp22_gag_compatible_fraction <- person$sevenp22_compatible_gag_units / pmax(person$sevenp22_total_cn, 1)
person$onep31b_internal_fragment_present <- as.integer(person$onep31b_total_cn > 0)
person$onep31b_any_array <- as.integer(person$onep31b_array_haplotypes > 0)

feature_specs <- data.frame(
  exposure_id = character(), family = character(), locus = character(), endpoint = character(),
  model_scale = character(), adjustment = character(), high_expression = logical(),
  prespecified_priority = integer(), source = character(), stringsAsFactors = FALSE
)
feature_values <- list()
feature_provirus <- list()

add_feature <- function(id, values, family, locus = "catalog-wide", endpoint = id,
                        scale = "continuous", adjustment = "", high = FALSE,
                        priority = 100L, source = "derived") {
  stopifnot(length(values) == length(direct_ids))
  feature_values[[id]] <<- setNames(as.numeric(values), direct_ids)
  feature_specs[nrow(feature_specs) + 1L, ] <<- list(id, family, locus, endpoint, scale,
                                                     adjustment, high, priority, source)
}

truth_keep <- truth[truth$status %in% candidate_status &
                      truth$statistical_hypothesis_role == "representative", ]
priority_table <- read_tsv(paths$high_expression)
priority_ids <- unique(priority_table$candidate_id[priority_table$high_expression_flag == "TRUE"])

truth_groups <- split(truth_keep, truth_keep$truth_equivalence_id)
for (group in truth_groups) {
  stopifnot(length(unique(group$locus)) == 1L, length(unique(group$contrast_id)) == 1L,
            !anyDuplicated(group$sample))
  locus <- unique(group$locus)
  domain <- unique(group$contrast_domain)
  endpoint <- unique(group$contrast_id)
  id <- paste(locus, domain, endpoint, sep = "::")
  x <- setNames(rep(NA_real_, length(direct_ids)), direct_ids)
  x[group$sample] <- group$focal_dosage
  if (domain == "structural") x <- corrected_structural_dosage(locus, endpoint)
  candidate_id <- paste0(locus, "::", domain, "::", endpoint)
  lp <- unname(type_map[locus])
  family <- paste0("locus_", ifelse(is.na(lp), "untyped", lp), "_", domain)
  provirus_adjustment <- if (domain == "orf") paste0(id, "::provirus_dosage_adjustment") else ""
  add_feature(id, x, family, locus, endpoint, "dosage", provirus_adjustment,
              candidate_id %in% priority_ids,
              ifelse(candidate_id == "HML-2_1q22::orf::orf_gag", 1L, 20L),
              "frozen Phase-0 representative truth")
  if (domain == "orf") {
    pv <- setNames(rep(NA_real_, length(direct_ids)), direct_ids)
    pv[group$sample] <- group$provirus_dosage
    feature_provirus[[id]] <- pv
  }
}

# Explicit catalog-wide Type-I/Type-II hypotheses.
add_feature("general::type1_physical_units", person$type1_units, "general_type_burden", priority = 2L)
add_feature("general::type2_physical_units", person$type2_units, "general_type_burden", priority = 3L)
add_feature("general::type1_loci_present", person$type1_loci_present, "general_type_burden", priority = 4L)
add_feature("general::type2_loci_present", person$type2_loci_present, "general_type_burden", priority = 5L)
add_feature("general::type1_fraction", person$type1_fraction, "general_type_burden", priority = 6L)
add_feature("general::type1_type2_log_ratio", person$type1_type2_log_ratio, "general_type_burden", priority = 7L)
add_feature("general::type1_minus_type2_units", person$type1_minus_type2_units, "general_type_burden", priority = 8L)
for (feature in c("gag", "pro", "pol", "env", "np9")) {
  column <- paste0("type1_compatible_", feature, "_units")
  add_feature(paste0("general::", column), person[[column]], "general_typeI_function",
              endpoint = column, adjustment = "general::type1_physical_units", priority = 10L)
}
for (feature in c("gag", "pro", "pol", "env", "rec")) {
  column <- paste0("type2_compatible_", feature, "_units")
  add_feature(paste0("general::", column), person[[column]], "general_typeII_function",
              endpoint = column, adjustment = "general::type2_physical_units", priority = 11L)
}

# Copy-resolved requested loci. Functional-unit effects are conditioned on physical CN.
add_feature("HML-2_7p22.1::total_physical_cn", person$sevenp22_total_cn, "sevenp22_array", "HML-2_7p22.1", priority = 1L)
add_feature("HML-2_7p22.1::array_haplotype_count", person$sevenp22_array_haplotypes, "sevenp22_array", "HML-2_7p22.1", scale = "dosage", priority = 2L)
add_feature("HML-2_7p22.1::any_array", person$sevenp22_any_array, "sevenp22_array", "HML-2_7p22.1", scale = "binary", priority = 3L)
add_feature("HML-2_7p22.1::max_haplotype_cn", person$sevenp22_max_haplotype_cn, "sevenp22_array", "HML-2_7p22.1", priority = 4L)
add_feature("HML-2_7p22.1::any_haplotype_cn_ge3", person$sevenp22_any_high_haplotype_ge3, "sevenp22_array", "HML-2_7p22.1", scale = "binary", priority = 5L)
add_feature("HML-2_7p22.1::compatible_gag_units", person$sevenp22_compatible_gag_units, "sevenp22_function", "HML-2_7p22.1", adjustment = "HML-2_7p22.1::total_physical_cn", priority = 6L)
add_feature("HML-2_7p22.1::fully_compatible_units", person$sevenp22_fully_compatible_units, "sevenp22_function", "HML-2_7p22.1", adjustment = "HML-2_7p22.1::total_physical_cn", priority = 7L)
add_feature("HML-2_7p22.1::gag_compatible_fraction", person$sevenp22_gag_compatible_fraction, "sevenp22_function", "HML-2_7p22.1", priority = 8L)

add_feature("HML-2_1p31.1b::internal_fragment_present", person$onep31b_internal_fragment_present, "onep31b_array", "HML-2_1p31.1b", scale = "binary", priority = 1L)
add_feature("HML-2_1p31.1b::total_internal_fragment_cn", person$onep31b_total_cn, "onep31b_array", "HML-2_1p31.1b", priority = 2L)
add_feature("HML-2_1p31.1b::any_array", person$onep31b_any_array, "onep31b_array", "HML-2_1p31.1b", scale = "binary", priority = 3L)
add_feature("HML-2_1p31.1b::max_haplotype_cn", person$onep31b_max_haplotype_cn, "onep31b_array", "HML-2_1p31.1b", priority = 4L)
add_feature("HML-2_1p31.1b::compatible_env_units", person$onep31b_compatible_env_units, "onep31b_function", "HML-2_1p31.1b", adjustment = "HML-2_1p31.1b::total_internal_fragment_cn", priority = 5L)
add_feature("HML-2_1p31.1b::compatible_np9_units", person$onep31b_compatible_np9_units, "onep31b_function", "HML-2_1p31.1b", adjustment = "HML-2_1p31.1b::total_internal_fragment_cn", priority = 6L)

# Materialize the auditable person-by-feature matrix.
feature_matrix <- data.frame(sample = direct_ids)
for (id in feature_specs$exposure_id) feature_matrix[[id]] <- feature_values[[id]][direct_ids]
for (id in names(feature_provirus)) {
  adjustment_id <- paste0(id, "::provirus_dosage_adjustment")
  feature_matrix[[adjustment_id]] <- feature_provirus[[id]][direct_ids]
}
analysis_matrix <- merge(base, person, by = "sample", all.x = TRUE)
analysis_matrix <- merge(analysis_matrix, feature_matrix, by = "sample", all.x = TRUE)
write_tsv(analysis_matrix, "person_level_direct_matrix.tsv")

hc3_fit <- function(data, exposure, adjustment_ids = character(), ancestry = c("superpopulation", "population")) {
  ancestry <- match.arg(ancestry)
  d <- data.frame(y = data$log2_ebv_load, x = data[[exposure]], sex = factor(data$sex),
                  ancestry = factor(if (ancestry == "superpopulation") data$superpopulation else data$population),
                  cluster = data$pedigree_component)
  if (length(adjustment_ids)) {
    for (j in seq_along(adjustment_ids)) d[[paste0("adj", j)]] <- data[[adjustment_ids[j]]]
  }
  d <- d[complete.cases(d), ]
  if (nrow(d) < 20L || length(unique(d$x)) < 2L || sd(d$x) == 0) return(list(status = "not_estimable", n = nrow(d)))
  adjustment_terms <- if (length(adjustment_ids)) paste0("adj", seq_along(adjustment_ids)) else character()
  rhs <- c("x", adjustment_terms, "ancestry", "sex")
  fit <- lm(as.formula(paste("y ~", paste(rhs, collapse = " + "))), data = d)
  X <- model.matrix(fit)
  pivot <- fit$qr$pivot[seq_len(fit$rank)]
  Xr <- X[, pivot, drop = FALSE]
  if (!("x" %in% colnames(Xr))) return(list(status = "exposure_aliased", n = nrow(d)))
  bread <- tryCatch(solve(crossprod(Xr)), error = function(e) NULL)
  if (is.null(bread)) return(list(status = "singular", n = nrow(d)))
  h <- rowSums((Xr %*% bread) * Xr)
  w <- residuals(fit)^2 / pmax((1 - h)^2, 1e-12)
  meat <- crossprod(Xr, Xr * w)
  vc <- bread %*% meat %*% bread
  ix <- match("x", colnames(Xr))
  beta <- coef(fit)["x"]
  se <- sqrt(vc[ix, ix])
  stat <- beta / se
  df <- fit$df.residual
  cluster_rows <- split(seq_len(nrow(d)), d$cluster)
  cluster_scores <- lapply(cluster_rows, function(rows) colSums(Xr[rows, , drop = FALSE] * residuals(fit)[rows]))
  cluster_meat <- Reduce(`+`, lapply(cluster_scores, tcrossprod))
  cluster_count <- length(cluster_scores)
  cluster_correction <- (cluster_count / (cluster_count - 1)) * ((nrow(d) - 1) / (nrow(d) - fit$rank))
  cluster_vc <- bread %*% cluster_meat %*% bread * cluster_correction
  cluster_se <- sqrt(cluster_vc[ix, ix])
  cluster_stat <- beta / cluster_se
  cluster_p <- 2 * pt(abs(cluster_stat), df = cluster_count - 1, lower.tail = FALSE)
  conventional <- summary(fit)$coefficients["x", ]
  list(status = "fit", n = nrow(d), beta = beta, se_hc3 = se, t_hc3 = stat,
       p_hc3 = 2 * pt(abs(stat), df = df, lower.tail = FALSE),
       se_cluster = cluster_se, p_cluster = cluster_p, clusters = cluster_count,
       cluster_ci_low = beta + qt(0.025, cluster_count - 1) * cluster_se,
       cluster_ci_high = beta + qt(0.975, cluster_count - 1) * cluster_se,
       ci_low = beta + qt(0.025, df) * se, ci_high = beta + qt(0.975, df) * se,
       se_conventional = conventional["Std. Error"], p_conventional = conventional["Pr(>|t|)"],
       df = df, standardized_beta = beta * sd(d$x) / sd(d$y))
}

overlap <- analysis_matrix[!is.na(analysis_matrix$ebv_load), ]
descriptives <- list()
model_rows <- list()
for (i in seq_len(nrow(feature_specs))) {
  spec <- feature_specs[i, ]
  id <- spec$exposure_id
  x <- overlap[[id]]
  known <- !is.na(x)
  nonzero <- sum(x[known] != 0)
  zero <- sum(x[known] == 0)
  counts <- table(x[known])
  rare_binary <- spec$model_scale %in% c("binary", "dosage") && min(zero, nonzero) < 5L
  status_gate <- if (sum(known) < 30L) "underpowered_n_lt30" else if (length(counts) < 2L) "invariant" else if (rare_binary) "underpowered_arm_lt5" else "fit"
  descriptives[[length(descriptives) + 1L]] <- data.frame(
    exposure_id = id, family = spec$family, locus = spec$locus, endpoint = spec$endpoint,
    model_scale = spec$model_scale, adjustment = spec$adjustment, high_expression = spec$high_expression,
    n = sum(known), n_missing = sum(!known), n_zero = zero, n_nonzero = nonzero,
    n_unique = length(counts), min = ifelse(sum(known), min(x[known]), NA),
    median = ifelse(sum(known), median(x[known]), NA), max = ifelse(sum(known), max(x[known]), NA),
    arm_counts = paste(names(counts), as.integer(counts), sep = ":", collapse = ";"),
    screening_gate = status_gate
  )
  adjustment_ids <- character()
  if (nzchar(spec$adjustment)) adjustment_ids <- spec$adjustment
  if (id %in% names(feature_provirus)) {
    adj_name <- paste0(id, "::provirus_dosage_adjustment")
    overlap[[adj_name]] <- feature_provirus[[id]][overlap$sample]
    if (!(adj_name %in% adjustment_ids)) adjustment_ids <- c(adjustment_ids, adj_name)
  }
  primary <- if (status_gate == "fit") hc3_fit(overlap, id, adjustment_ids, "superpopulation") else list(status = status_gate, n = sum(known))
  sensitivity <- if (status_gate == "fit") hc3_fit(overlap, id, adjustment_ids, "population") else list(status = status_gate, n = sum(known))
  getv <- function(object, field) if (is.null(object[[field]])) NA else object[[field]]
  model_rows[[length(model_rows) + 1L]] <- data.frame(
    exposure_id = id, family = spec$family, locus = spec$locus, endpoint = spec$endpoint,
    model_scale = spec$model_scale, adjustment = paste(adjustment_ids, collapse = ";"),
    high_expression = spec$high_expression, n = getv(primary, "n"), model_status = primary$status,
    beta_log2_ebv = getv(primary, "beta"), se_hc3 = getv(primary, "se_hc3"),
    se_pedigree_cluster = getv(primary, "se_cluster"), p_pedigree_cluster = getv(primary, "p_cluster"),
    pedigree_clusters = getv(primary, "clusters"),
    ci95_low = getv(primary, "cluster_ci_low"), ci95_high = getv(primary, "cluster_ci_high"),
    p_hc3 = getv(primary, "p_hc3"), p_conventional = getv(primary, "p_conventional"),
    standardized_beta = getv(primary, "standardized_beta"),
    population_model_status = sensitivity$status, population_beta = getv(sensitivity, "beta"),
    population_se_cluster = getv(sensitivity, "se_cluster"), population_p_cluster = getv(sensitivity, "p_cluster"),
    stringsAsFactors = FALSE
  )
}
descriptives <- do.call(rbind, descriptives)
models <- do.call(rbind, model_rows)

# Collapse exact direct-overlap exposure/adjustment truth vectors before global correction.
vector_key <- function(id) {
  x <- overlap[[id]]
  spec <- feature_specs[feature_specs$exposure_id == id, ]
  adj <- spec$adjustment
  if (nzchar(adj)) adjustment_vector <- overlap[[adj]] else adjustment_vector <- rep(NA_real_, nrow(overlap))
  paste(spec$model_scale,
        paste(ifelse(is.na(x), "NA", format(x, scientific = FALSE, trim = TRUE)), collapse = ","),
        paste(ifelse(is.na(adjustment_vector), "NA", format(adjustment_vector, scientific = FALSE, trim = TRUE)), collapse = ","), sep = "|")
}
keys <- setNames(vapply(feature_specs$exposure_id, vector_key, character(1)), feature_specs$exposure_id)
feature_specs$observed_vector_group <- match(keys, unique(keys))
feature_specs <- feature_specs[order(feature_specs$observed_vector_group, feature_specs$prespecified_priority, feature_specs$exposure_id), ]
feature_specs$multiplicity_representative <- !duplicated(feature_specs$observed_vector_group)
representative_ids <- feature_specs$exposure_id[feature_specs$multiplicity_representative]

models <- merge(models, feature_specs[, c("exposure_id", "observed_vector_group", "multiplicity_representative")], by = "exposure_id", all.x = TRUE)
models$q_bh_global <- NA_real_
models$p_holm_global <- NA_real_
rep_fit <- models$multiplicity_representative & models$model_status == "fit" & is.finite(models$p_pedigree_cluster)
models$q_bh_global[rep_fit] <- p.adjust(models$p_pedigree_cluster[rep_fit], method = "BH")
models$p_holm_global[rep_fit] <- p.adjust(models$p_pedigree_cluster[rep_fit], method = "holm")
for (g in unique(models$observed_vector_group)) {
  rows <- which(models$observed_vector_group == g)
  rep_row <- rows[models$multiplicity_representative[rows]][1]
  if (!is.na(rep_row)) {
    models$q_bh_global[rows] <- models$q_bh_global[rep_row]
    models$p_holm_global[rows] <- models$p_holm_global[rep_row]
  }
}
models$q_bh_family <- NA_real_
for (fam in unique(models$family)) {
  ix <- which(models$family == fam & models$multiplicity_representative & models$model_status == "fit" & is.finite(models$p_pedigree_cluster))
  if (length(ix)) models$q_bh_family[ix] <- p.adjust(models$p_pedigree_cluster[ix], method = "BH")
}

alias_map <- merge(feature_specs[, c("exposure_id", "observed_vector_group", "multiplicity_representative")],
                   feature_specs[feature_specs$multiplicity_representative, c("observed_vector_group", "exposure_id")],
                   by = "observed_vector_group", suffixes = c("", "_representative"), all.x = TRUE)

# Does a catalog-wide Type-I/II unit burden reduce to one locus? Report each
# variable locus's contribution and the leave-one-locus-out burden result.
decomposition <- list()
for (tp in c("type1", "type2")) {
  loci <- if (tp == "type1") type1_loci else type2_loci
  total_id <- paste0("general::", tp, "_physical_units")
  total <- setNames(person[[paste0(tp, "_units")]], person$sample)
  for (locus in loci) {
    contribution <- setNames(unit_table[, locus], rownames(unit_table))
    overlap[["leave_one_burden"]] <- total[overlap$sample] - contribution[overlap$sample]
    fit <- hc3_fit(overlap, "leave_one_burden", character(), "superpopulation")
    decomposition[[length(decomposition) + 1L]] <- data.frame(
      burden = total_id, locus = locus, locus_type = tp,
      mean_locus_units = mean(contribution[overlap$sample]),
      sd_locus_units = sd(contribution[overlap$sample]),
      correlation_with_total = suppressWarnings(cor(contribution[overlap$sample], total[overlap$sample])),
      leave_one_model_status = fit$status,
      leave_one_beta = ifelse(is.null(fit$beta), NA, fit$beta),
      leave_one_p_cluster = ifelse(is.null(fit$p_cluster), NA, fit$p_cluster)
    )
  }
}
decomposition <- do.call(rbind, decomposition)

write_tsv(descriptives[order(descriptives$screening_gate, descriptives$family, descriptives$exposure_id), ], "exposure_descriptives.tsv")
write_tsv(models[order(models$q_bh_global, models$p_pedigree_cluster, na.last = TRUE), ], "model_results.tsv")
write_tsv(alias_map[order(alias_map$observed_vector_group, alias_map$exposure_id), ], "truth_vector_alias_map.tsv")
write_tsv(decomposition, "type_burden_leave_one_locus_out.tsv")
write_tsv(feature_specs, "feature_manifest.tsv")

limitations <- data.frame(
  item = c(
    "causal_attribution", "local_haplotype", "global_type_interpretation", "rare_high_7p22_copy",
    "1p31.1b_function_vs_copy", "pedigree", "outcome_transform", "absence_semantics",
    "Houldcroft_replication", "proliferation_replication"
  ),
  status = c(
    "untested", "untested", "descriptive_burden_only", "underpowered", "nonidentifiable",
    ifelse(related_component_count == 0, "no-repeated-components", "cluster-robust-adjustment-used"),
    "tested_log2", "implemented", "pending_acquisition", "pending_acquisition"
  ),
  reason_or_needed_source = c(
    "Sequence state can tag a linked host haplotype; association is not causal evidence.",
    "No local phased SNP/ancestry matrix is available in the 116-person phenotype subset.",
    "Type is fixed by locus, so cross-locus Type-I/Type-II burden is inseparable from locus composition.",
    "In this overlap, 7p22.1 person CN is 2-5 and only one person has CN5; no high-copy tail test is reliable.",
    "At 1p31.1b, Env/Np9-compatible-unit count equals internal-fragment CN in the current calls.",
    paste0("Repeated pedigree components in overlap: ", related_component_count, "; primary inference uses pedigree-cluster-robust covariance."),
    "All values are positive; primary outcome is log2(copies per cell). Raw outcome remains in the person matrix.",
    "No catalog row contributes zero units (provirus absent), not a missing genotype, for catalog-wide burdens.",
    "Acquire and QC the PLOS One 2014 qPCR supplement, resolve replicate rows, then rerun the frozen exposures.",
    "Extract the PLOS Genetics 2012 proliferation phenotype/sample identifiers, then rerun the frozen exposures."
  )
)
write_tsv(limitations, "limitations_and_pending.tsv")

summary <- data.frame(
  metric = c("source_phenotype_people", "direct_sequence_people", "direct_phenotype_overlap",
             "overlap_superpopulations", "overlap_populations", "related_pedigree_components",
             "exposures_total", "unique_direct_truth_vectors", "models_fit",
             "global_bh_lt_0.05", "global_holm_lt_0.05"),
  value = c(nrow(phenotype), length(direct_ids), nrow(overlap),
            length(unique(overlap$superpopulation)), length(unique(overlap$population)), related_component_count,
            nrow(feature_specs), length(unique(feature_specs$observed_vector_group)),
            sum(models$model_status == "fit" & models$multiplicity_representative),
            sum(models$q_bh_global < 0.05, na.rm = TRUE), sum(models$p_holm_global < 0.05, na.rm = TRUE))
)
write_tsv(summary, "run_summary.tsv")

cat("Mandage direct EBV screen complete\n")
print(summary)
cat("Top fitted representative results:\n")
print(head(models[models$multiplicity_representative & models$model_status == "fit",
                  c("exposure_id", "n", "beta_log2_ebv", "p_pedigree_cluster", "q_bh_global")][order(models$p_pedigree_cluster[models$multiplicity_representative & models$model_status == "fit"]), ], 12))
