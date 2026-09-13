source(Sys.getenv("HML2_CONFIG", file.path("R", "config.R")))
# >>> Generates: Figures S3A, S3B <<<

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
                 "HML-2_22q11.21","HML-2_6q14.1a","HML-2_19q12")

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
  mutate(locus = if_else(locus == "HML-2_7p22.1a", "HML-2_7p22.1", locus))

analysis_data <- cleaned_data %>% filter(locus %in% target_loci)
cat("Filtered data to", nrow(analysis_data), "rows for target loci.\n")

# --- 5. Calculate Denominator for Frequencies ---
# Denominator is now the number of PROVIRUS alleles per locus.
total_provirus_alleles_per_locus <- analysis_data %>%
  # (MODIFIED) Handle NA in 'Structure' just in case
  filter(tolower(replace_na(Structure, "")) == "provirus") %>%
  count(locus, name = "total_alleles")
cat("Calculated provirus allele counts for frequency denominator.\n")

# --- 6. Parse, Classify, and Summarize Mutations ---
cat("--- Starting Mutation Analysis ---\n")
cat("Step 1: Pivoting data to long format...\n")
mutations_long <- analysis_data %>%
  filter(tolower(replace_na(Structure, "")) == "provirus") %>%
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

cat("Step 2: Classifying mutations and detecting missed frameshifts...\n")

all_classified_mutations <- mutations_long %>%
  separate_rows(mutation_list, sep = ",\\s*") %>%
  rename(mutation = mutation_list) %>%
  filter(mutation != "") %>%
  # Skip masked-region markers ('Undetermined:START-END' / 'Undetermined:full'):
  # these are unresolved spans, not real mutations, so drop them before parsing/plotting.
  filter(!str_starts(str_to_lower(mutation), "undetermined")) %>%
  mutate(
    position = case_when(
      str_detect(mutation, "^Frameshift_at_") ~ as.integer(str_extract(mutation, "(?<=Frameshift_at_)\\d+")),
      TRUE ~ parse_position(mutation)
    ),
    type = case_when(
      str_detect(mutation, "^Frameshift_at_") ~ "Frameshift",
      str_detect(mutation, fixed("*"))        ~ "Nonsense",
      TRUE                                    ~ "Missense"
    ),
    mutation = if_else(
      str_detect(mutation, "^Frameshift_at_"),
      paste0(position, "_fs"),
      mutation
    )
  ) %>%
  filter(!is.na(position))

# Fallback heuristic: for any ID_Full+orf combo with no explicit Frameshift_at_ tag,
# check for a dense cluster of mutations as a proxy for an untagged frameshift.
cat("Step 2b: Applying fallback frameshift heuristic to untagged ORFs...\n")

already_tagged <- all_classified_mutations %>%
  filter(type == "Frameshift") %>%
  distinct(ID_Full, orf)

untagged_orfs <- all_classified_mutations %>%
  anti_join(already_tagged, by = c("ID_Full", "orf")) %>%
  group_by(ID_Full, orf, locus) %>%
  summarise(mutation_list = paste(mutation, collapse = ","), .groups = "drop")

detect_frameshift_sequential_heuristic <- function(mutation_list_str) {
  if (is.na(mutation_list_str) || mutation_list_str == "") {
    return(list(is_frameshift = FALSE, frameshift_start_position = NA_integer_))
  }
  mutations <- str_split(mutation_list_str, ",\\s*")[[1]]
  mutations <- mutations[mutations != ""]
  positions <- as.integer(str_extract(mutations, "\\d+"))
  df <- tibble(mutation = mutations, position = positions) %>%
    filter(!is.na(position)) %>%
    arrange(position)
  if (nrow(df) < 12) return(list(is_frameshift = FALSE, frameshift_start_position = NA_integer_))
  for (i in 1:(nrow(df) - 11)) {
    window_end <- df$position[i] + 14
    if (sum(df$position >= df$position[i] & df$position <= window_end) >= 12) {
      return(list(is_frameshift = TRUE, frameshift_start_position = df$position[i]))
    }
  }
  return(list(is_frameshift = FALSE, frameshift_start_position = NA_integer_))
}

heuristic_results <- untagged_orfs %>%
  mutate(fs_info = map(mutation_list, detect_frameshift_sequential_heuristic)) %>%
  unnest_wider(fs_info) %>%
  filter(is_frameshift) %>%
  select(ID_Full, orf, frameshift_start_position)

cat("Heuristic detected", nrow(heuristic_results), "additional untagged frameshifts.\n")

# Inject synthetic frameshift markers into all_classified_mutations
synthetic_fs <- heuristic_results %>%
  mutate(
    mutation = paste0(frameshift_start_position, "_fs"),
    position = frameshift_start_position,
    type = "Frameshift"
  )

# Add locus and orf_group to synthetic rows by joining back
synthetic_fs <- synthetic_fs %>%
  left_join(distinct(all_classified_mutations, ID_Full, orf, locus, orf_group),
            by = c("ID_Full", "orf"))

# NEW: remove stop codons that coincide with heuristic frameshift starts
all_classified_mutations <- all_classified_mutations %>%
  anti_join(
    heuristic_results %>%
      rename(position = frameshift_start_position) %>%
      select(ID_Full, orf, position),
    by = c("ID_Full", "orf", "position")
  )

all_classified_mutations <- bind_rows(
  all_classified_mutations,
  synthetic_fs %>% select(-frameshift_start_position)
)

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

# (NEW) Create dummy data to force axis scales from 0 to max length
orf_max_lengths <- tibble(
  orf = c("gag", "pro", "pol", "env"),
  max_pos = c(560, 150, 950, 580) # Approximate full lengths
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
  # Individual missense labels aren't biologically meaningful here, so the panels
  # show the mutation *landscape* (position vs frequency) without per-point labels.
  facet_grid(locus ~ factor(orf, levels = facet_order), scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = "Missense mutations", x = "Amino-acid position", y = "Allele frequency") +
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
  scale_color_manual(name = "Mutation type",
                     values = c("Nonsense" = "#D55E00", "Frameshift" = "#CC79A7", "Missense" = "#009E73")) +
  labs(title = "Disabling mutations", x = "Amino-acid position", y = "Allele frequency") +
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
