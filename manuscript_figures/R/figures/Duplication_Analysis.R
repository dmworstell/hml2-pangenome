source(Sys.getenv("HML2_CONFIG", file.path("R", "config.R")))
# >>> Generates: Figures 3B, S4A, S4B <<<

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

cleaned_data <- full_data %>% filter(!is.na(locus) & !is.na(ID_Full))

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

# Create a consistent, colour-blind-safe palette for all possible duplication sizes.
# Array size is an ordered discrete class (2x..6x) -> Okabe-Ito categorical palette.
all_categories <- c("2x", "3x", "4x", "5x", "6x")
present_categories <- intersect(all_categories, unique(figure1_plot_data$category))
oi_seq <- unname(okabe_ito[c("sky_blue", "bluish_green", "orange", "vermillion", "reddish_purple")])
color_palette <- oi_seq[seq_along(present_categories)]
names(color_palette) <- present_categories

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
    intact_freq = mean(is_intact, na.rm = TRUE),
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
  geom_line(linewidth = 1) +
  geom_point(size = 2.6) +
  facet_wrap(~ locus, labeller = labeller(locus = function(x) str_replace(x, "HML-2_", ""))) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_color_manual(name = "ORF", values = unname(okabe_ito[c("blue", "orange", "bluish_green", "vermillion", "reddish_purple", "sky_blue")])) +
  theme_pub(base_size = 16) +
  labs(
    x = "Position in tandem array (part number)",
    y = "Frequency of intact ORF"
  )

# Figure 3: Relative array-size makeup of the duplicated alleles at each locus.
cat("Generating Figure 3: Relative Array Size Frequency Plot...\n")
plot_relsize <- relative_freq_data %>%
  filter(locus %in% loci_with_duplications) %>%
  mutate(locus = factor(locus, levels = locus_order_fig1)) %>%
  ggplot(aes(x = locus, y = relative_frequency, fill = category)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = function(x) str_replace(x, "HML-2_", "")) +
  scale_fill_manual(name = "Tandem array size", values = color_palette, limits = names(color_palette), drop = FALSE) +
  theme_pub(base_size = 16) +
  labs(
    x = "HML-2 locus",
    y = "Relative frequency"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- 9. Save all figures (PNG + PDF) ---
save_fig(plot_fig1,      "hml2_duplication_landscape_ONLY",     width = 14, height = 8)
save_fig(plot_positional,"hml2_positional_orf_integrity",       width = 12, height = 10)
save_fig(plot_relsize,   "hml2_duplication_relative_frequency", width = 14, height = 8)

cat("--- Script Finished ---\n")
