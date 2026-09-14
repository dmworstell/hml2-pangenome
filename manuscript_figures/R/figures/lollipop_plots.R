source(Sys.getenv("HML2_CONFIG", file.path("R", "config.R")))
# Generates the supplemental amino-acid substitution and frameshift-position panels.

# HML-2 Locus-Specific Mutation Lollipop Plot Generator
#
# This script reads an aggregated TSV file of HML-2 ORF data and generates
# a faceted lollipop plot to visualize the frequency and position of missense,
# nonsense, and frameshift mutations for a specified set of HML-2 loci.
#
# -- MODIFIED FOR POSTER (v2) --
# 1. Reduced target loci list.
# 2. Removed 'np9/K-rev' column from plots.
# 3. Increased all font and point sizes.
# 4. Removed subtitles and updated main title.
# 5. Adjusted output plot dimensions.
# 6. (NEW) Increased axis/legend font sizes.
# 7. (NEW) Forced x-axes to start at 0 and go to max ORF length.
# 8. (NEW) Increased ggrepel max.overlaps.
# -------------------------

# --- 1. Load Necessary Libraries ---
# Ensure you have these libraries installed: install.packages(c("tidyverse", "viridis", "ggrepel"))
library(tidyverse)
library(viridis)
library(ggrepel) # For non-overlapping text labels

# --- 2. User-Defined Parameters ---
# IMPORTANT: Update this path to the location of your aggregated TSV file.
input_file_path <- HML2_ORF_TABLE

# Define the output directory where the plots will be saved.
output_dir <- HML2_FIG_DIR

# (MODIFIED) Define the loci you want to visualize
target_loci <- c("HML-2_7p22.1", "HML-2_19p12b", "HML-2_10p12.1", "HML-2_5q33.3","HML-2_1q22",
                 "HML-2_22q11.21","HML-2_6q14.1","HML-2_19q12")

# --- 3. Helper Functions ---
# Parses the amino acid position from a mutation string (e.g., P108L -> 108)
parse_position <- function(mutation_str) {
  # This is a fully vectorized version of the function that works with dplyr.
  as.integer(stringr::str_extract(mutation_str, "\\d+"))
}



# --- 4. Load and Prepare the Data ---
cat("Loading and preparing data from:", input_file_path, "\n")
# (MODIFIED) Load all columns as character to avoid parsing warnings
full_data <- read_tsv(
  input_file_path,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
) %>% rename(locus = Locus)

required_columns <- c("analysis_include", "analysis_exclusion_reason", "ID", "ID_Full",
                      "locus", "Structure", "missense_gag", "missense_pro", "missense_pol",
                      "missense_env", "missense_np9", "np9", "rec")
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
input_row_count <- nrow(full_data)
full_data <- public_data %>% filter(analysis_include == "1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(tibble(
  input_catalog = normalizePath(input_file_path),
  input_sha256 = digest::digest(file = input_file_path, algo = "sha256"),
  input_rows = input_row_count,
  nonpublic_rows_excluded = input_row_count - nrow(public_data),
  public_rows_excluded = sum(public_data$analysis_include == "0"),
  retained_public_rows = nrow(full_data)
), file.path(output_dir, "hml2_variant_position_input_provenance.tsv"))

full_data <- full_data %>%
  mutate(
    np9 = if_else(tolower(np9) == "n", NA_character_, np9),
    rec  = if_else(tolower(rec)  == "n", NA_character_, rec),
    provirus_type = case_when(
      !is.na(np9) ~ "type1",
      !is.na(rec)  ~ "type2",
      TRUE         ~ NA_character_
    )
  )

# --- (NEW) Filter out rows with 100% ComparisonError_ProteinMissing ---
cat("Checking for 'ComparisonError_ProteinMissing' rows in missense columns...\n")
n_before <- nrow(full_data)

# (FIX) Use dplyr::filter() and replace_na()
full_data_filtered <- full_data %>%
  filter(!(
    replace_na(missense_gag, "") == "ComparisonError_ProteinMissing" &
      replace_na(missense_pro, "") == "ComparisonError_ProteinMissing" &
      replace_na(missense_pol, "") == "ComparisonError_ProteinMissing" &
      replace_na(missense_env, "") == "ComparisonError_ProteinMissing" &
      replace_na(missense_np9, "") == "ComparisonError_ProteinMissing"
  ))

n_after <- nrow(full_data_filtered)
cat("Removed", n_before - n_after, "rows where all missense columns were 'ComparisonError_ProteinMissing'.\n")
# --- End of new block ---

cleaned_data <- full_data_filtered %>%
  filter(!is.na(locus) & !is.na(ID_Full)) %>%
  mutate(locus = if_else(locus %in% c("HML-2_7p22.1a", "HML-2_7p22.1b"),
                         "HML-2_7p22.1", locus))

analysis_data <- cleaned_data %>%
  filter(locus %in% target_loci, Structure %in% c("Provirus", "Provirus_from_Multi"))
if (nrow(analysis_data) != nrow(distinct(analysis_data, locus, ID_Full))) {
  stop("A retained proviral copy occurs more than once in a plotted locus.")
}
cat("Filtered data to", nrow(analysis_data), "rows for target loci.\n")
write_tsv(analysis_data %>% select(locus, ID, Haplotype, ID_Full, Structure,
                                  starts_with("missense_")),
          file.path(output_dir, "hml2_variant_position_retained_members.tsv"))

# --- 5. Calculate Denominator for Frequencies ---
# Denominator is now the number of PROVIRUS alleles per locus.
total_provirus_alleles_per_locus <- analysis_data %>%
  count(locus, name = "total_alleles")
write_tsv(total_provirus_alleles_per_locus,
          file.path(output_dir, "hml2_variant_position_denominators.tsv"))
cat("Calculated provirus allele counts for frequency denominator.\n")

# --- 6. Parse, Classify, and Summarize Mutations ---
cat("--- Starting Mutation Analysis ---\n")
cat("Step 1: Pivoting data to long format...\n")
mutations_long <- analysis_data %>%
  select(locus, ID_Full, provirus_type, starts_with("missense_")) %>%
  pivot_longer(
    cols = starts_with("missense_"),
    names_to = "orf",
    values_to = "mutation_list",
    names_prefix = "missense_"
  ) %>%
  # (MODIFIED) Filter out NAs, blanks, AND individual protein errors
  filter(
    !is.na(mutation_list) &
      mutation_list != "" &
      mutation_list != "ComparisonError_ProteinMissing"
  ) %>%
  # (MODIFIED) Create combined np9/K-rev column first for filtering
  mutate(orf_group = if_else(orf %in% c("np9", "rec"), "np9/K-rev", orf))

mutations_long <- mutations_long %>%
  mutate(mutation_count = str_count(mutation_list, ",") + 1)

cat("Step 2: Parsing explicit substitutions and frameshift annotations...\n")

all_classified_mutations <- mutations_long %>%
  separate_rows(mutation_list, sep = ",\\s*") %>%
  rename(mutation = mutation_list) %>%
  filter(mutation != "") %>%
  # Skip masked-region markers ('Undetermined:START-END' / 'Undetermined:full'):
  # these are unresolved spans, not real mutations, so drop them before parsing/plotting.
  filter(!str_starts(str_to_lower(mutation), "undetermined")) %>%
  # These panels plot substitutions and explicitly annotated frameshift starts.
  # In-frame indel annotations remain in the exported source-member table.
  filter(str_detect(mutation, "^[A-Z*][0-9]+[A-Z*]$") |
           str_detect(mutation, "^Frameshift_at_[0-9]+(-Premature_Stop)?$")) %>%
  mutate(
    position = case_when(
      str_detect(mutation, "^Frameshift_at_") ~ as.integer(str_extract(mutation, "(?<=Frameshift_at_)\\d+")),
      TRUE ~ parse_position(mutation)
    ),
    type = case_when(
      str_detect(mutation, "^Frameshift_at_") ~ "Frameshift",
      str_detect(mutation, "\\*$")           ~ "Nonsense",
      TRUE                                    ~ "Missense"
    ),
    mutation = if_else(
      str_detect(mutation, "^Frameshift_at_"),
      paste0(position, "_fs"),
      mutation
    )
  ) %>%
  filter(!is.na(position)) %>%
  distinct(locus, ID_Full, orf, mutation, position, type, .keep_all = TRUE)

cat("Step 3: Filtering mutations downstream of explicit frameshifts...\n")
frameshift_positions <- all_classified_mutations %>%
  filter(type == "Frameshift") %>%
  group_by(ID_Full, orf) %>%
  summarise(frameshift_start_position = min(position), .groups = "drop")

mutations_filtered <- all_classified_mutations %>%
  left_join(frameshift_positions, by = c("ID_Full", "orf")) %>%
  filter(is.na(frameshift_start_position) | position <= frameshift_start_position)

n_removed <- nrow(all_classified_mutations) - nrow(mutations_filtered)
cat("Removed", n_removed, "mutations downstream of detected frameshifts.\n")

cat("Step 5: Summarizing final mutation counts and frequencies...\n")
mutation_summary <- mutations_filtered %>%
  count(locus, orf_group, orf, mutation, position, type, name = "count") %>%
  left_join(total_provirus_alleles_per_locus, by = "locus") %>%
  mutate(frequency = count / total_alleles)
if (any(mutation_summary$count > mutation_summary$total_alleles)) {
  stop("A variant numerator exceeds its retained proviral-copy denominator.")
}
write_tsv(mutations_filtered %>% select(locus, ID_Full, orf, mutation, position, type),
          file.path(output_dir, "hml2_variant_position_observations.tsv"))
write_tsv(mutation_summary, file.path(output_dir, "hml2_variant_position_summary.tsv"))

cat("Mutation parsing and frequency calculation complete.\n")


# --- 7. Special Analysis: Y195C in 7p22.1 pol ---
cat("\n--- Special Analysis ---\n")
y195c_data <- mutation_summary %>%
  filter(locus == "HML-2_7p22.1", orf == "pol", mutation == "Y195C")
if (nrow(y195c_data) > 0) {
  cat("Frequency of Y195C in HML-2_7p22.1 (pol):\n")
  cat(scales::percent(y195c_data$frequency, accuracy = 0.01)," (", y195c_data$count, " out of ", y195c_data$total_alleles, " provirus alleles)\n", sep = "")
} else {
  cat("Mutation Y195C was not found in HML-2_7p22.1 (pol).\n")
}
cat("-----------------------\n\n")

# --- 8. Generate and Save Plots ---
cat("Generating plots...\n")

# (MODIFIED) Define common properties
facet_order <- c("gag", "pro", "pol", "env") # Removed np9/K-rev
label_threshold <- 0.01

# Display a common origin and the maximum observed coordinate for each ORF.
# These are annotation coordinates, not approximate asserted protein lengths.
orf_max_lengths <- tibble(
  orf = facet_order
) %>% left_join(
  mutation_summary %>% group_by(orf) %>% summarise(max_pos = max(position), .groups = "drop"),
  by = "orf"
)
# Create a data frame with 0 and max_pos for every panel in the grid
orf_scale_data <- crossing(locus = target_loci, orf = facet_order) %>%
  left_join(orf_max_lengths, by = "orf") %>%
  pivot_longer(cols = c(max_pos), values_to = "position") %>%
  bind_rows(crossing(locus = target_loci, orf = facet_order, position = 0)) %>%
  select(locus, orf, position) %>%
  # (NEW) Remove HML-2_ prefix for plotting
  mutate(locus = str_replace(locus, "HML-2_", ""))


# --- Plot 1: Missense Mutations ---
missense_data <- mutation_summary %>%
  filter(type == "Missense", orf_group %in% facet_order) %>%
  # (NEW) Remove HML-2_ prefix for plotting
  mutate(locus = str_replace(locus, "HML-2_", ""))

plot_missense <- missense_data %>%
  mutate(locus = factor(locus, levels = str_replace(target_loci, "HML-2_", ""))) %>%
  ggplot(aes(x = position, y = frequency)) +
  geom_blank(data = orf_scale_data, aes(x = position, y = 0)) +
  geom_segment(aes(xend = position, yend = 0), color = "#0072B2", linewidth = 0.45) +
  geom_point(color = "#0072B2", size = 2.1) +
  # Exact substitution labels are retained in the source table. Labeling every
  # dense point in these panels would obscure their positions and frequencies.
  facet_grid(locus ~ factor(orf, levels = facet_order), scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = "Amino-acid substitutions", x = "KCON amino-acid position", y = "Fraction of proviral copies") +
  theme_pub(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.5),
        panel.spacing = grid::unit(2, "pt"),
        strip.background = element_rect(fill = "grey92", colour = "grey35", linewidth = 0.5),
        strip.text = element_text(face = "plain", size = 17))

# --- Plot 2: Disabling Mutations ---
disabling_data <- mutation_summary %>%
  filter(
    (type %in% c("Nonsense", "Frameshift") |
       (locus == "HML-2_7p22.1" & orf == "pol" & mutation == "Y195C")) &
      orf_group %in% facet_order
  ) %>%
  # (NEW) Remove HML-2_ prefix for plotting
  mutate(locus = str_replace(locus, "HML-2_", ""))

plot_disabling <- disabling_data %>%
  mutate(locus = factor(locus, levels = str_replace(target_loci, "HML-2_", ""))) %>%
  ggplot(aes(x = position, y = frequency)) +
  geom_blank(data = orf_scale_data, aes(x = position, y = 0)) +
  geom_segment(aes(xend = position, yend = 0, color = type), linewidth = 0.45) +
  geom_point(aes(color = type), size = 2.1) +
  geom_text_repel(
    data = . %>% group_by(locus, orf) %>% slice_max(frequency, n = 4, with_ties = FALSE) %>% ungroup(),
    aes(label = mutation),
    colour = "grey15", size = 5.2, min.segment.length = 0, max.overlaps = 6,
    segment.size = 0.25, segment.colour = "grey65", box.padding = 0.55, seed = 1
  ) +
  facet_grid(locus ~ factor(orf, levels = facet_order), scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.05, 0.12))) +
  scale_color_manual(name = NULL,
                     values = c("Nonsense" = "#D55E00", "Frameshift" = "#CC79A7", "Missense" = "#009E73"),
                     breaks = c("Frameshift", "Missense", "Nonsense"),
                     labels = c("Annotated frameshift", "Pol Y195C", "Stop-gain")) +
  labs(title = "Frameshift positions and Pol Y195C", x = "KCON amino-acid position", y = "Fraction of proviral copies") +
  # Full-width supplement panel (scaled ~0.67); bump base + the in-panel mutation labels
  # (*_fs / nonsense) which were flagged as the smallest text, ~5 pt / below.
  theme_pub(base_size = 17, legend = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.5),
        panel.spacing = grid::unit(2, "pt"),
        strip.background = element_rect(fill = "grey92", colour = "grey35", linewidth = 0.5),
        strip.text = element_text(face = "plain", size = 17))

# --- 9. Save Plots ---
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_fig(plot_missense,  "hml2_missense_lollipop_plot",            width = 12, height = 11, dir = output_dir)
save_fig(plot_disabling, "hml2_disabling_mutations_lollipop_plot", width = 12, height = 11, dir = output_dir)

cat("--- Script Finished ---\n")
