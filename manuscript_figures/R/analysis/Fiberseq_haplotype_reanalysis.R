#!/usr/bin/env Rscript

# Reanalyse HML-2 Fiber-seq actuation at the biological unit of observation.
#
# The collaborator peak files contain one row per sample run x haplotype x
# consensus peak. Treating those rows as independent inflates the apparent
# sample size because a person can contribute two haplotypes, several peaks,
# and more than one sequencing run/platform. This script therefore:
#   1. pools peak-level FIRE and total coverage within each haplotype;
#   2. retains the collaborator-designated primary sample for inference;
#   3. summarises haplotypes within individuals and individuals within loci;
#   4. uses non-primary runs only for a technical/platform concordance audit.
#
# It intentionally makes no claim that chromatin actuation is transcription.

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, env, default = NULL) {
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) return(args[[hit + 1]])
  value <- Sys.getenv(env, unset = "")
  if (nzchar(value)) return(value)
  default
}

peaks_dir <- arg_value("--peaks-dir", "HML2_FIBERSEQ_PEAKS_DIR")
out_dir <- arg_value("--out-dir", "HML2_FIBERSEQ_REANALYSIS_OUT", "outputs/fiberseq_reanalysis")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
crosswalk_file <- arg_value(
  "--locus-crosswalk", "HML2_FIBERSEQ_LOCUS_CROSSWALK",
  file.path(dirname(script_file), "Fiberseq_locus_crosswalk.tsv")
)
locus_crosswalk <- read_tsv(crosswalk_file, show_col_types = FALSE, progress = FALSE)
crosswalk_keys <- paste(locus_crosswalk$assay_label, locus_crosswalk$assay_reference_interval)
if (anyDuplicated(crosswalk_keys)) stop("Duplicate assay intervals in the locus crosswalk.")
focus_locus <- arg_value("--focus-locus", "HML2_FIBERSEQ_FOCUS_LOCUS", "1q22")
minimum_coverage <- as.numeric(arg_value("--minimum-coverage", "HML2_FIBERSEQ_MINIMUM_COVERAGE", "10"))
if (!is.finite(minimum_coverage) || minimum_coverage < 1 || minimum_coverage != floor(minimum_coverage)) {
  stop("--minimum-coverage must be a positive integer.")
}
active_threshold <- as.numeric(arg_value(
  "--active-threshold", "HML2_FIBERSEQ_ACTIVE_THRESHOLD", "0.10"
))

if (is.null(peaks_dir) || !dir.exists(peaks_dir)) {
  stop("Provide the collaborator peak-file directory with --peaks-dir or HML2_FIBERSEQ_PEAKS_DIR.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

extract_assay_label <- function(path) {
  name <- basename(path)
  locus <- str_match(name, "_HML-2_(.+)_peaks_file\\.tsv$")[, 2]
  ifelse(is.na(locus), sub("_peaks_file\\.tsv$", "", name), locus)
}

extract_locus <- function(path) {
  assay_label <- extract_assay_label(path)
  interval <- str_match(basename(path), "^TestHML_(chr[^_]+)_([0-9]+)-([0-9]+)_")
  key <- paste(assay_label, paste0(interval[, 2], ":", interval[, 3], "-", interval[, 4]))
  index <- match(key, crosswalk_keys)
  if (anyNA(index)) stop("No coordinate-matched locus crosswalk entry for ", path)
  catalog_locus <- locus_crosswalk$catalog_locus[index]
  ifelse(is.na(catalog_locus) | !nzchar(catalog_locus), assay_label, catalog_locus)
}

extract_assay_interval <- function(path) {
  hit <- str_match(basename(path), "^TestHML_[^_]+_([0-9]+)-([0-9]+)_HML-2_")
  tibble(
    assay_start = suppressWarnings(as.integer(hit[, 2])),
    assay_end = suppressWarnings(as.integer(hit[, 3]))
  )
}

read_peak_file <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  interval <- extract_assay_interval(path)
  required <- c(
    "cons_start", "cons_end", "consensus_peak_id", "sample_id",
    "coverage", "fire_coverage", "Individual_ID", "Haplotype", "Platform"
  )
  absent <- setdiff(required, names(x))
  if (length(absent)) {
    stop("Missing required columns in ", path, ": ", paste(absent, collapse = ", "))
  }

  if (!"is_primary_sample" %in% names(x)) x$is_primary_sample <- TRUE
  if (!"meets_cutoff" %in% names(x)) x$meets_cutoff <- TRUE
  if (!"PS" %in% names(x)) x$PS <- NA_character_

  x %>%
    transmute(
      locus = extract_locus(path),
      assay_start = interval$assay_start[[1]],
      assay_end = interval$assay_end[[1]],
      Individual_ID = as.character(Individual_ID),
      Haplotype = as.character(Haplotype),
      Platform = as.character(Platform),
      PS = as.character(PS),
      sample_id = as.character(sample_id),
      is_primary_sample = as_flag(is_primary_sample),
      meets_cutoff = as_flag(meets_cutoff),
      consensus_peak_id = as.character(consensus_peak_id),
      cons_start = as.integer(cons_start),
      cons_end = as.integer(cons_end),
      coverage = as.numeric(coverage),
      fire_coverage = as.numeric(fire_coverage)
    )
}

files <- list.files(peaks_dir, pattern = "_peaks_file\\.tsv$", full.names = TRUE)
if (!length(files)) stop("No *_peaks_file.tsv inputs found in ", peaks_dir)

file_info <- file.info(files)
input_manifest <- tibble(
  input_file = basename(files),
  assay_label = map_chr(files, extract_assay_label),
  locus = map_chr(files, extract_locus),
  bytes = as.numeric(file_info$size),
  modified_utc = format(file_info$mtime, tz = "UTC", usetz = TRUE),
  sha256 = map_chr(files, ~ digest(.x, algo = "sha256", file = TRUE))
) %>%
  arrange(locus, input_file)
write_tsv(input_manifest, file.path(out_dir, "fiberseq_input_manifest.tsv"))
write_tsv(locus_crosswalk, file.path(out_dir, "Figure_S17_locus_crosswalk.tsv"))
write_tsv(
  tibble(
    parameter = c(
      "analysis_version", "focus_locus", "active_haplotype_threshold", "minimum_peak_coverage",
      "input_peak_files", "biological_sample_rule", "technical_sample_rule"
    ),
    value = c(
      "1.2.0", focus_locus, as.character(active_threshold), as.character(minimum_coverage),
      as.character(length(files)),
      "collaborator-designated primary runs collapsed peak-to-haplotype-to-individual",
      "non-primary runs excluded from biological counts and used only for concordance"
    )
  ),
  file.path(out_dir, "fiberseq_analysis_parameters.tsv")
)

message("Reading ", length(files), " Fiber-seq peak files...")
raw <- map_dfr(files, read_peak_file) %>%
  filter(
    is.finite(coverage), coverage >= minimum_coverage,
    is.finite(fire_coverage), fire_coverage >= 0,
    !is.na(Individual_ID), nzchar(Individual_ID),
    !is.na(Haplotype), nzchar(Haplotype)
  )

# Defensive collapse in case an input contains duplicate records for the same
# sample-run/haplotype/consensus-peak combination.
peak_run <- raw %>%
  group_by(
    locus, assay_start, assay_end, Individual_ID, Haplotype, Platform, PS, sample_id,
    is_primary_sample, consensus_peak_id, cons_start, cons_end
  ) %>%
  summarise(
    coverage = sum(coverage),
    fire_coverage = sum(fire_coverage),
    actuation = fire_coverage / coverage,
    .groups = "drop"
  )

run_summary <- peak_run %>%
  group_by(locus, Individual_ID, Haplotype, Platform, PS, sample_id, is_primary_sample) %>%
  summarise(
    n_consensus_peaks = n_distinct(consensus_peak_id),
    total_coverage = sum(coverage),
    total_fire_coverage = sum(fire_coverage),
    pooled_actuation = total_fire_coverage / total_coverage,
    mean_peak_actuation = mean(actuation),
    maximum_peak_actuation = max(actuation),
    .groups = "drop"
  )

primary_runs <- run_summary %>% filter(is_primary_sample)
if (!nrow(primary_runs)) stop("No rows are marked as primary samples.")

# There should be one primary run per locus/haplotype. Pooling here is a safe
# fallback if the source design contains more than one primary-labelled run.
haplotype_summary <- primary_runs %>%
  group_by(locus, Individual_ID, Haplotype) %>%
  summarise(
    n_primary_runs = n(),
    n_consensus_peaks = sum(n_consensus_peaks),
    total_coverage = sum(total_coverage),
    total_fire_coverage = sum(total_fire_coverage),
    pooled_actuation = total_fire_coverage / total_coverage,
    maximum_peak_actuation = max(maximum_peak_actuation),
    .groups = "drop"
  )

if (any(haplotype_summary$total_fire_coverage > haplotype_summary$total_coverage)) {
  stop("Primary-haplotype FIRE coverage exceeds total coverage.")
}

individual_summary <- haplotype_summary %>%
  group_by(locus, Individual_ID) %>%
  summarise(
    n_haplotypes = n_distinct(Haplotype),
    mean_haplotype_actuation = mean(pooled_actuation),
    maximum_haplotype_actuation = max(pooled_actuation),
    total_coverage = sum(total_coverage),
    total_fire_coverage = sum(total_fire_coverage),
    pooled_individual_actuation = total_fire_coverage / total_coverage,
    .groups = "drop"
  )

locus_summary <- haplotype_summary %>%
  group_by(locus) %>%
  summarise(
    n_individuals = n_distinct(Individual_ID),
    n_haplotypes = n(),
    total_coverage = sum(total_coverage),
    total_fire_coverage = sum(total_fire_coverage),
    pooled_locus_actuation = total_fire_coverage / total_coverage,
    median_haplotype_actuation = median(.data$pooled_actuation),
    q1_haplotype_actuation = quantile(.data$pooled_actuation, 0.25),
    q3_haplotype_actuation = quantile(.data$pooled_actuation, 0.75),
    maximum_haplotype_actuation = max(.data$pooled_actuation),
    active_haplotype_threshold = active_threshold,
    active_haplotype_fraction = mean(.data$pooled_actuation >= active_threshold),
    .groups = "drop"
  ) %>%
  left_join(
    individual_summary %>%
      group_by(locus) %>%
      summarise(
        median_individual_actuation = median(pooled_individual_actuation),
        q1_individual_actuation = quantile(pooled_individual_actuation, 0.25),
        q3_individual_actuation = quantile(pooled_individual_actuation, 0.75),
        .groups = "drop"
      ),
    by = "locus"
  ) %>%
  # Backward-compatible fields used by Public_expression_concordance.R.
  mutate(
    med = median_individual_actuation,
    mx = maximum_haplotype_actuation,
    pct_active = active_haplotype_fraction,
    n = n_individuals
  )

# Technical/platform audit. The primary designation is not used to cherry-pick
# a value here; it only establishes which run is compared with the alternate.
replicate_pairs <- run_summary %>%
  group_by(locus, Individual_ID, Haplotype) %>%
  filter(n() >= 2) %>%
  arrange(desc(is_primary_sample), Platform, sample_id, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  mutate(pair_member = c("primary", "alternate")) %>%
  ungroup() %>%
  select(
    locus, Individual_ID, Haplotype, pair_member, Platform, sample_id,
    pooled_actuation, total_coverage
  ) %>%
  pivot_wider(
    names_from = pair_member,
    values_from = c(Platform, sample_id, pooled_actuation, total_coverage),
    names_sep = "_"
  ) %>%
  filter(is.finite(pooled_actuation_primary), is.finite(pooled_actuation_alternate))

replicate_statistics <- if (nrow(replicate_pairs) >= 3) {
  within_individual <- replicate_pairs %>%
    group_by(Individual_ID) %>%
    filter(n() >= 3) %>%
    summarise(
      n_haplotype_locus_pairs = n(),
      spearman_rho = suppressWarnings(cor(
        pooled_actuation_primary, pooled_actuation_alternate, method = "spearman"
      )),
      median_absolute_difference = median(abs(
        pooled_actuation_primary - pooled_actuation_alternate
      )),
      .groups = "drop"
    )
  tibble(
    n_pairs = nrow(replicate_pairs),
    n_individuals = n_distinct(replicate_pairs$Individual_ID),
    spearman_rho = suppressWarnings(cor(
      replicate_pairs$pooled_actuation_primary,
      replicate_pairs$pooled_actuation_alternate,
      method = "spearman"
    )),
    pearson_r = suppressWarnings(cor(
      replicate_pairs$pooled_actuation_primary,
      replicate_pairs$pooled_actuation_alternate,
      method = "pearson"
    )),
    median_absolute_difference = median(abs(
      replicate_pairs$pooled_actuation_primary - replicate_pairs$pooled_actuation_alternate
    )),
    median_within_individual_spearman = median(within_individual$spearman_rho, na.rm = TRUE),
    minimum_within_individual_spearman = min(within_individual$spearman_rho, na.rm = TRUE),
    maximum_within_individual_spearman = max(within_individual$spearman_rho, na.rm = TRUE)
  )
} else {
  tibble(
    n_pairs = nrow(replicate_pairs),
    n_individuals = n_distinct(replicate_pairs$Individual_ID),
    spearman_rho = NA_real_,
    pearson_r = NA_real_,
    median_absolute_difference = NA_real_,
    median_within_individual_spearman = NA_real_,
    minimum_within_individual_spearman = NA_real_,
    maximum_within_individual_spearman = NA_real_
  )
}

write_tsv(haplotype_summary, file.path(out_dir, "fiberseq_haplotype_summary.tsv"))
write_tsv(individual_summary, file.path(out_dir, "fiberseq_individual_summary.tsv"))
write_tsv(locus_summary, file.path(out_dir, "fiberseq_locus_summary.tsv"))
write_tsv(haplotype_summary, file.path(out_dir, "Figure_S17_source_haplotype_summary.tsv"))
write_tsv(locus_summary, file.path(out_dir, "Figure_S17_source_locus_summary.tsv"))
write_tsv(replicate_pairs, file.path(out_dir, "fiberseq_technical_replicate_pairs.tsv"))
write_tsv(replicate_statistics, file.path(out_dir, "fiberseq_technical_replicate_statistics.tsv"))

# The manuscript supplement needs one distribution panel, using one point for
# each primary haplotype rather than treating peak tiles or technical runs as
# independent observations.
manuscript_locus_order <- locus_summary %>%
  arrange(median_haplotype_actuation, locus) %>%
  pull(locus)
manuscript_labels <- locus_summary %>%
  transmute(locus, label = paste0(locus, " (", n_haplotypes, ")"))
manuscript_panel <- haplotype_summary %>%
  left_join(manuscript_labels, by = "locus") %>%
  mutate(locus = factor(locus, levels = manuscript_locus_order)) %>%
  ggplot(aes(pooled_actuation, locus)) +
  geom_boxplot(
    width = 0.56, outlier.shape = NA, linewidth = 0.35,
    fill = "grey94", colour = "grey40"
  ) +
  geom_point(
    position = position_jitter(height = 0.13, width = 0, seed = 20260914),
    alpha = 0.55, size = 0.9, colour = "#426779"
  ) +
  scale_y_discrete(labels = setNames(manuscript_labels$label, manuscript_labels$locus)) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2), limits = c(0, 1),
    labels = function(x) paste0(round(x * 100), "%"),
    expand = expansion(mult = c(0.005, 0.01))
  ) +
  labs(
    title = "HML-2 chromatin accessibility",
    x = "FIRE coverage / total coverage per haplotype", y = NULL
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    text = element_text(colour = "black", family = "sans"),
    axis.text = element_text(colour = "black", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    plot.title = element_text(size = 12.5, face = "plain", hjust = 0.5),
    plot.margin = margin(6, 9, 6, 6)
  )
ggsave(
  file.path(out_dir, "hml2_fiberseq_locus_accessibility.png"),
  manuscript_panel, width = 7.1, height = 5.8, dpi = 450, bg = "white"
)
ggsave(
  file.path(out_dir, "hml2_fiberseq_locus_accessibility.pdf"),
  manuscript_panel, width = 7.1, height = 5.8, device = cairo_pdf, bg = "white"
)

# ---- Figure -----------------------------------------------------------------
theme_manuscript <- theme_minimal(base_size = 12.5) +
  theme(
    text = element_text(colour = "black", family = "sans"),
    axis.text = element_text(colour = "grey20"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.35),
    plot.title = element_text(size = 14),
    plot.subtitle = element_text(size = 10.5, colour = "grey30"),
    legend.position = "none",
    plot.margin = margin(8, 10, 8, 8)
  )

ordered_loci <- locus_summary %>%
  arrange(median_individual_actuation) %>%
  pull(locus)

p_a <- individual_summary %>%
  mutate(locus = factor(locus, levels = ordered_loci)) %>%
  ggplot(aes(pooled_individual_actuation * 100, locus)) +
  geom_boxplot(
    width = 0.58, outlier.shape = NA, linewidth = 0.35,
    fill = "grey92", colour = "grey35"
  ) +
  geom_point(
    aes(colour = locus == focus_locus),
    position = position_jitter(height = 0.12, width = 0),
    alpha = 0.45, size = 1.1
  ) +
  scale_colour_manual(values = c(`FALSE` = "#567189", `TRUE` = "#C54B3C")) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "A  Fiber-seq actuation after biological-unit collapse",
    subtitle = "One value per individual and locus; technical runs and peak tiles are not independent replicates",
    x = "Pooled FIRE coverage / total coverage",
    y = NULL
  ) +
  theme_manuscript

focus_peaks <- peak_run %>%
  filter(locus == focus_locus, is_primary_sample) %>%
  group_by(cons_start, cons_end, consensus_peak_id) %>%
  summarise(
    n_haplotypes = n_distinct(paste(Individual_ID, Haplotype)),
    coverage = sum(coverage),
    fire_coverage = sum(fire_coverage),
    pooled_actuation = fire_coverage / coverage,
    .groups = "drop"
  )

if (nrow(focus_peaks)) {
  provirus_start <- suppressWarnings(min(raw$assay_start[raw$locus == focus_locus], na.rm = TRUE))
  provirus_end <- suppressWarnings(max(raw$assay_end[raw$locus == focus_locus], na.rm = TRUE))
  focus_mid <- (focus_peaks$cons_start + focus_peaks$cons_end) / 2
  x_pad <- max(500, 0.08 * (provirus_end - provirus_start))
  focus_breaks <- pretty(c(provirus_start, provirus_end), n = 4)
  focus_breaks <- focus_breaks[focus_breaks >= provirus_start & focus_breaks <= provirus_end]

  p_b <- ggplot(focus_peaks, aes(focus_mid, pooled_actuation * 100)) +
    annotate(
      "rect", xmin = provirus_start, xmax = provirus_end,
      ymin = -Inf, ymax = Inf, fill = "#EAF0F4", alpha = 0.65
    ) +
    geom_segment(
      aes(x = cons_start, xend = cons_end, yend = pooled_actuation * 100),
      linewidth = 4.5, lineend = "butt", colour = "#C54B3C"
    ) +
    geom_point(size = 2.7, colour = "#8A2D24") +
    annotate(
      "text", x = provirus_end + x_pad * 0.85, y = Inf,
      label = "upstream GON4L/LTR12F\nregion not assayed",
      hjust = 0.5, vjust = 1.25, size = 3.1, colour = "grey30"
    ) +
    coord_cartesian(xlim = c(provirus_start - x_pad * 0.25, provirus_end + x_pad * 1.55)) +
    scale_x_continuous(
      breaks = focus_breaks,
      labels = scales::label_number(big.mark = ",")
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.02, 0.18))) +
    labs(
      title = paste0("B  ", focus_locus, " actuation is localized within the assayed interval"),
      subtitle = "Peak-level pooled coverage across primary haplotypes; the negative-strand upstream promoter falls beyond the window",
      x = "GRCh38 coordinate",
      y = "Pooled actuation"
    ) +
    theme_manuscript
} else {
  p_b <- ggplot() +
    annotate("text", x = 0, y = 0, label = paste("No peak data for", focus_locus)) +
    theme_void() +
    labs(title = paste0("B  ", focus_locus, " interval"))
}

if (nrow(replicate_pairs) >= 3) {
  qc <- replicate_statistics[1, ]
  p_c <- ggplot(
    replicate_pairs,
    aes(pooled_actuation_primary * 100, pooled_actuation_alternate * 100)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_point(alpha = 0.55, size = 1.6, colour = "#3C6E8F") +
    annotate(
      "label", x = Inf, y = -Inf,
      label = paste0(
        "n = ", qc$n_pairs, " haplotype-locus pairs\n",
        qc$n_individuals, " individuals; median within-person rho = ",
        formatC(qc$median_within_individual_spearman, digits = 2, format = "f"), "\n",
        "median |difference| = ",
        formatC(qc$median_absolute_difference * 100, digits = 1, format = "f"), " points"
      ),
      hjust = 1.05, vjust = -0.15, size = 3.2, linewidth = 0.2
    ) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = "C  Independent run/platform concordance",
      subtitle = "Alternate runs are reserved for technical validation and excluded from biological sample counts",
      x = "Primary-run pooled actuation",
      y = "Alternate-run pooled actuation"
    ) +
    theme_manuscript
} else {
  p_c <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Too few independent run pairs for concordance") +
    theme_void() +
    labs(title = "C  Independent run/platform concordance")
}

focus_haps <- haplotype_summary %>% filter(locus == focus_locus)
p_d <- ggplot(focus_haps, aes(pooled_actuation * 100)) +
  geom_histogram(binwidth = 2.5, boundary = 0, fill = "#C8D8E4", colour = "white") +
  geom_rug(alpha = 0.45, colour = "#8A2D24") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = paste0("D  ", focus_locus, " varies among haplotypes"),
    subtitle = paste0(
      "n = ", nrow(focus_haps), " primary haplotypes from ",
      n_distinct(focus_haps$Individual_ID), " individuals"
    ),
    x = "Pooled haplotype actuation",
    y = "Haplotypes"
  ) +
  theme_manuscript

right_column <- (p_b / p_c / p_d) +
  plot_layout(heights = c(0.82, 1, 0.78))

figure <- (p_a | right_column) +
  plot_layout(widths = c(0.95, 1.45)) +
  plot_annotation(
    title = "HML-2 Fiber-seq actuation with biological and technical replication separated",
    theme = theme(plot.title = element_text(size = 17, face = "plain", margin = margin(b = 8)))
  )

ggsave(
  file.path(out_dir, "hml2_fiberseq_haplotype_reanalysis.png"),
  figure, width = 15.5, height = 11.5, dpi = 320, bg = "white"
)
ggsave(
  file.path(out_dir, "hml2_fiberseq_haplotype_reanalysis.pdf"),
  figure, width = 15.5, height = 11.5, device = cairo_pdf, bg = "white"
)

message("Wrote Fiber-seq reanalysis to ", normalizePath(out_dir))
