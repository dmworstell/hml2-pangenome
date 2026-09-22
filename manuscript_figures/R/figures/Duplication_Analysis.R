source(Sys.getenv("HML2_CONFIG", file.path("R", "config.R")))
# Generates the tandem-array positional ORF and conditional array-size panels.

# HML-2 Tandem Duplication & ORF Integrity Analysis
#
# This script reads an aggregated TSV file of HML-2 ORF data to analyze loci
# containing tandem duplications (e.g., MULTI_part1, MULTI_part2).
# It generates two plots:
# 1. A stacked bar chart showing the frequency and size of tandem arrays per locus.
# 2. A faceted line plot showing how ORF integrity changes by position within the arrays.

# --- 1. Load Necessary Libraries ---
# Ensure you have these libraries installed: install.packages(c("tidyverse", "viridis", "patchwork"))
library(tidyverse)
library(viridis)
library(patchwork) # For combining plots

# --- 2. User-Defined Parameters ---
# IMPORTANT: Update this path to the location of your aggregated TSV file.
input_file_path <- HML2_ORF_TABLE

# Define the output directory where the plots will be saved.
output_dir <- HML2_FIG_DIR

# --- 3. Load and Prepare the Data ---
cat("Loading and preparing data from:", input_file_path, "\n")
full_data <- read_tsv(
  input_file_path,
  col_types = cols(.default = "c"), # This is the key change
  show_col_types = FALSE
)
full_data <- full_data %>% rename(locus = Locus)

required_columns <- c("analysis_include", "analysis_exclusion_reason", "ID", "Haplotype",
                      "ID_Full", "locus", "Structure", "gag", "pro", "pol", "env", "np9", "rec")
missing_columns <- setdiff(required_columns, names(full_data))
if (length(missing_columns)) {
  stop("The retained analysis catalog is required. Missing columns: ",
       paste(missing_columns, collapse = ", "))
}
public_data <- full_data %>% filter(str_detect(ID, "^(HG|NA)[0-9]+$"))
if (any(is.na(public_data$analysis_include)) ||
    any(!public_data$analysis_include %in% c("0", "1"))) {
  stop("Every public catalog row must explicitly declare analysis_include as 0 or 1.")
}
cleaned_data <- public_data %>%
  filter(analysis_include == "1", !is.na(locus), !is.na(ID_Full))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(
  public_data %>% filter(analysis_include == "0") %>%
    count(analysis_exclusion_reason, name = "excluded_rows"),
  file.path(output_dir, "hml2_tandem_excluded_public_rows.tsv")
)
write_tsv(tibble(
  input_catalog = normalizePath(input_file_path),
  input_sha256 = digest::digest(file = input_file_path, algo = "sha256"),
  input_rows = nrow(full_data),
  nonpublic_rows_excluded = nrow(full_data) - nrow(public_data),
  public_rows = nrow(public_data),
  included_rows = nrow(cleaned_data),
  included_haplotypes = nrow(distinct(cleaned_data, ID, Haplotype))
), file.path(output_dir, "hml2_tandem_input_provenance.tsv"))

cleaned_data <- cleaned_data %>%
  mutate(
    np9 = if_else(tolower(np9) == "n", NA_character_, np9),
    rec  = if_else(tolower(rec)  == "n", NA_character_, rec),
    provirus_type = case_when(
      !is.na(np9) ~ "type1",
      !is.na(rec)  ~ "type2",
      TRUE         ~ NA_character_
    )
  )

cat("Loaded and cleaned", nrow(cleaned_data), "rows.\n")

# --- 3b. Data Cleaning ---
cat("Performing data cleaning steps (e.g., renaming loci)...\n")
cleaned_data <- cleaned_data %>%
  mutate(
    # (FIXED) Correctly collapse 7p22.1a and 7p22.1b
    locus = case_when(
      locus %in% c("HML-2_7p22.1a", "HML-2_7p22.1b") ~ "HML-2_7p22.1",
      TRUE ~ locus
    )
  )

# --- 4. Analyze Allelic Structure Frequency & Size (for Figure 1) ---
cat("Analyzing allelic structure frequency...\n")

# First, calculate the total number of unique haplotypes for each locus. This is our denominator.
total_haplotypes_per_locus <- cleaned_data %>%
  distinct(locus, ID, Haplotype) %>%
  count(locus, name = "total_haplotypes")

# Next, identify and process tandem duplications. Tandem units carry a _part# suffix
# regardless of the preceding tag, so match _part# directly to catch _MULTI_part#,
# legacy _DOUBLE_part#, and nested _alt#_part# (a tandem inside a segdup copy). This
# (correctly) excludes _alt#-only segdup copies and _asmdup artifacts, which aren't tandems.
duplication_data <- cleaned_data %>%
  filter(str_detect(ID_Full, "_part\\d")) %>%
  mutate(part_num = as.integer(str_extract(ID_Full, "(?<=_part)\\d+")))
if (any(!duplication_data$Structure %in% c("Provirus", "Provirus_from_Multi"))) {
  stop("A retained tandem part is not a proviral unit. Inspect its structural record.")
}
array_membership <- duplication_data %>%
  group_by(locus, ID, Haplotype) %>%
  summarise(member_count = n(), distinct_positions = n_distinct(part_num),
            first_position = min(part_num), last_position = max(part_num), .groups = "drop")
if (any(array_membership$member_count != array_membership$distinct_positions) ||
    any(array_membership$first_position != 1L) ||
    any(array_membership$last_position != array_membership$member_count)) {
  stop("Tandem positions must be unique and contiguous within each locus-haplotype array.")
}
write_tsv(
  duplication_data %>% select(locus, ID, Haplotype, ID_Full, Structure, part_num,
                             provirus_type, gag, pro, pol, env, np9, rec),
  file.path(output_dir, "hml2_tandem_retained_members.tsv")
)

# Find the size of each tandem array.
array_size_summary <- duplication_data %>%
  group_by(locus, ID, Haplotype) %>%
  summarise(array_size = max(part_num, na.rm = TRUE), .groups = 'drop')

# Count how many arrays of each size exist for each locus.
array_size_counts <- array_size_summary %>%
  count(locus, array_size, name = "count") %>%
  mutate(category = factor(paste0(array_size, "x"))) %>%
  select(locus, category, count)

# **NOTE**: We are now only focusing on duplications for Figure 1.
# We will calculate their frequency relative to the total haplotype count.
duplication_freq_data <- array_size_counts %>%
  left_join(total_haplotypes_per_locus, by = "locus") %>%
  mutate(frequency = count / total_haplotypes)

# --- 5. Generate Figure 1: Duplication-Only Landscape Bar Chart ---
cat("Generating Figure 1: Duplication-Only Landscape Bar Chart...\n")

# Identify loci that have at least one tandem duplication.
loci_with_duplications <- unique(array_size_summary$locus)

# Filter the frequency data to only include these loci for Figure 1.
figure1_plot_data <- duplication_freq_data %>%
  filter(locus %in% loci_with_duplications)

# Order the filtered loci by their total frequency of duplications
locus_order_fig1 <- figure1_plot_data %>%
  group_by(locus) %>%
  summarise(total_freq = sum(frequency, na.rm = TRUE)) %>%
  arrange(desc(total_freq)) %>%
  pull(locus)

# Array-size counts use an ordered grayscale, distinct from gene identities.
all_categories <- c("2x", "3x", "4x", "5x", "6x")
present_categories <- intersect(all_categories, unique(figure1_plot_data$category))
size_palette <- c(`2x` = "#D9D9D9", `3x` = "#B0B0B0", `4x` = "#808080", `5x` = "#555555", `6x` = "#252525")
color_palette <- size_palette[present_categories]

# Split data for a faceted plot approach
data_fig1_top <- figure1_plot_data %>%
  filter(locus == "HML-2_7p22.1")

data_fig1_others <- figure1_plot_data %>%
  filter(locus != "HML-2_7p22.1") %>%
  mutate(locus = factor(locus, levels = locus_order_fig1[locus_order_fig1 != "HML-2_7p22.1"]))

# Plot A: Just the high-frequency locus
plot_A <- data_fig1_top %>%
  ggplot(aes(x = locus, y = frequency, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  # (MODIFIED) Add scale_x_discrete to remove prefix
  scale_x_discrete(labels = function(x) str_replace(x, "HML-2_", "")) +
  scale_fill_manual(name = "Tandem\narray\nsize", values = color_palette, limits = names(color_palette), drop = FALSE) +
  # (MODIFIED) Label is already correct, x is NULL
  labs(y = "Frequency among all haplotypes", x = NULL) +
  theme_pub(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(face = "plain"),
    legend.position = "right"
  )

# Plot B: The rest of the loci
plot_B <- data_fig1_others %>%
  ggplot(aes(x = locus, y = frequency, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  # (MODIFIED) Add scale_x_discrete to remove prefix
  scale_x_discrete(labels = function(x) str_replace(x, "HML-2_", "")) +
  scale_fill_manual(name = "Tandem array size", values = color_palette, limits = names(color_palette), drop = FALSE) +
  # (MODIFIED) Remove x-axis label
  labs(y = NULL, x = NULL) +
  theme_pub(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

# Combine plots with patchwork. No on-panel title/subtitle (NAR: description
# belongs in the figure legend). Panel letters are stamped by the composite.
plot_fig1 <- plot_A + plot_B +
  plot_layout(widths = c(1.5, 8), guides = "collect")


# --- 6. Analyze Positional ORF Integrity (Figure 2) ---
cat("Analyzing ORF integrity by position within tandem arrays...\n")

# Define intact criteria
intact_criteria <- c("intact", "no_stop", "no_stop_fs_end", "frameshift_at_end", "intact_fs_end", "intact_fs_end_premature_stop")

# Process the duplication data to get ORF status for each part
positional_orf_data <- duplication_data %>%
  mutate(
    rec = if("rec" %in% names(.)) rec else NA_character_,
    secondary_orf_status = coalesce(np9, rec),
    # Undetermined = uncallable (KCON-misalignment gate), NOT broken -> NA, dropped by the
    # existing `mean(is_intact, na.rm = TRUE)` below rather than scored as zero.
    gag = ifelse(tolower(gag) == "undetermined", NA, tolower(gag) %in% intact_criteria),
    pro = ifelse(tolower(pro) == "undetermined", NA, tolower(pro) %in% intact_criteria),
    pol = ifelse(tolower(pol) == "undetermined", NA, tolower(pol) %in% intact_criteria),
    env = ifelse(tolower(env) == "undetermined", NA, tolower(env) %in% intact_criteria),
    secondary_orf = ifelse(tolower(secondary_orf_status) == "undetermined", NA,
                           tolower(secondary_orf_status) %in% intact_criteria)
  ) %>%
  select(locus, part_num, provirus_type, gag, pro, pol, env, secondary_orf) %>%
  pivot_longer(
    cols = c(gag, pro, pol, env, secondary_orf),
    names_to = "orf",
    values_to = "is_intact"
  ) %>%
  mutate(
    orf = case_when(
      orf == "secondary_orf" & provirus_type == "type1" ~ "np9",
      orf == "secondary_orf" & provirus_type == "type2" ~ "K-rev",
      TRUE ~ orf
    )
  )

# Calculate the frequency of intactness for each ORF at each position, for each locus
positional_freq_summary <- positional_orf_data %>%
  group_by(locus, orf, part_num) %>%
  summarise(
    compatible_copies = sum(is_intact, na.rm = TRUE),
    assessed_copies = sum(!is.na(is_intact)),
    unassessed_copies = sum(is.na(is_intact)),
    intact_freq = ifelse(assessed_copies > 0, compatible_copies / assessed_copies, NA_real_),
    .groups = 'drop'
  )

# --- 7. Analyze Relative Duplication Size (Figure 3) ---
cat("Analyzing relative frequency of array sizes...\n")

# Re-use 'array_size_counts' from Section 4; the denominator is the total number
# of *duplicated* haplotypes per locus (so each locus's bars sum to 100%).
total_duplicated_haplotypes <- array_size_counts %>%
  group_by(locus) %>%
  summarise(total_duplicated = sum(count, na.rm = TRUE))

relative_freq_data <- array_size_counts %>%
  left_join(total_duplicated_haplotypes, by = "locus") %>%
  mutate(relative_frequency = count / total_duplicated)
write_tsv(array_size_summary, file.path(output_dir, "hml2_tandem_array_sizes.tsv"))
write_tsv(relative_freq_data, file.path(output_dir, "hml2_tandem_relative_size_summary.tsv"))
write_tsv(positional_freq_summary, file.path(output_dir, "hml2_tandem_positional_orf_summary.tsv"))
write_tsv(tibble(compatible_status = intact_criteria),
          file.path(output_dir, "hml2_tandem_positional_orf_criteria.tsv"))

# --- 8. Build the three figures ---
# Figure 1 (plot_fig1) was built in Section 5. Here we build Figure 2 (positional
# ORF integrity) and Figure 3 (relative array-size makeup) as DISTINCT objects --
# previously both were assigned to `plot_fig2`, so the first was silently clobbered
# and Figure 1 was saved twice. Each is now built once and saved once.

# Figure 2: Positional ORF integrity for the top loci with duplications.
cat("Generating Figure 2: Positional ORF Integrity Plot...\n")
top_loci_with_duplications <- head(locus_order_fig1[locus_order_fig1 %in% loci_with_duplications], 9)

plot_positional <- positional_freq_summary %>%
  filter(locus %in% top_loci_with_duplications) %>%
  mutate(locus = factor(locus, levels = top_loci_with_duplications)) %>%
  ggplot(aes(x = factor(part_num), y = intact_freq, group = orf, color = orf)) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 1.7, na.rm = TRUE) +
  facet_wrap(~ locus, labeller = labeller(locus = function(x) str_replace(x, "HML-2_", ""))) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_color_manual(name = "ORF", values = c(
    gag = "#527F4C", pro = "#444C56", pol = "#735394",
    env = "#519AC4", np9 = "#A75A37", `K-rev` = "#7D4C76")) +
  theme_pub(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(size = 11),
        axis.text = element_text(size = 10)) +
  labs(
    x = "Position in tandem array",
    y = "Copies meeting ORF criterion"
  )

# Figure 3: Relative array-size makeup of the duplicated alleles at each locus.
cat("Generating Figure 3: Relative Array Size Frequency Plot...\n")
plot_relsize <- relative_freq_data %>%
  filter(locus %in% loci_with_duplications) %>%
  mutate(locus = factor(locus, levels = locus_order_fig1)) %>%
  ggplot(aes(x = locus, y = relative_frequency, fill = category)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), oob = scales::squish,
                     expand = expansion(mult = c(0, 0.01))) +
  scale_x_discrete(limits = rev(locus_order_fig1), labels = function(x) {
    counts <- total_duplicated_haplotypes$total_duplicated[match(x, total_duplicated_haplotypes$locus)]
    paste0(str_replace(x, "HML-2_", ""), " (n = ", counts, ")")
  }) +
  scale_fill_manual(name = "Tandem array size", values = color_palette, limits = names(color_palette), drop = FALSE) +
  coord_flip() +
  theme_pub(base_size = 11) +
  labs(
    x = NULL,
    y = "Fraction of array-bearing haplotypes"
  ) +
  theme(legend.position = "bottom", axis.text = element_text(size = 10),
        plot.margin = margin(5.5, 14, 5.5, 5.5))

# --- 9. Save all figures (PNG + PDF) ---
save_fig(plot_fig1,      "hml2_duplication_landscape_ONLY",     width = 14, height = 8)
save_fig(plot_positional,"hml2_positional_orf_integrity",       width = 7.1, height = 5.4, dpi = 450)
save_fig(plot_relsize,   "hml2_duplication_relative_frequency", width = 7.1, height = 4.8, dpi = 450)

cat("--- Script Finished ---\n")
