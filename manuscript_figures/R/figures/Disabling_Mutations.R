source(Sys.getenv("HML2_CONFIG", file.path("manuscript_figures", "R", "config.R")))
# >>> Generates: Figure S10A <<<

# HML-2 Disabling Mutations Stacked Bar Chart Generator
#
# This script identifies the primary disabling mutation for each defective provirus.
# It uses explicit upstream labels (if available) and falls back to a sliding-window
# heuristic to detect missed frameshifts.

# --- 1. Load Necessary Libraries ---
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(data.table)

# --- 2. User-Defined Parameters ---
input_file_path <- HML2_ORF_TABLE
output_dir <- HML2_FIG_DIR

# --- 3. Helper Functions (The Backup Heuristic) ---
parse_mutation_position <- function(mutation_str) {
  as.integer(stringr::str_extract(mutation_str, "\\d+"))
}

# Backup heuristic: detects frameshifts based on high mutation density
detect_frameshift_sequential_heuristic <- function(mutation_list_str) {
  if (is.na(mutation_list_str) || mutation_list_str == "" || mutation_list_str == "ComparisonError_ProteinMissing") {
    return(list(is_frameshift = FALSE, frameshift_start_mutation = NA_character_))
  }

  mutations <- str_split(mutation_list_str, ",\\s*")[[1]]
  mutations <- mutations[mutations != ""]

  if (length(mutations) < 12) {
    return(list(is_frameshift = FALSE, frameshift_start_mutation = NA_character_))
  }

  mutation_data <- tibble(mutation = mutations, position = parse_mutation_position(mutations)) %>%
    filter(!is.na(position)) %>%
    arrange(position)

  if (nrow(mutation_data) < 12) {
    return(list(is_frameshift = FALSE, frameshift_start_mutation = NA_character_))
  }

  for (i in 1:(nrow(mutation_data) - 11)) {
    window_start_pos <- mutation_data$position[i]
    window_end_pos <- window_start_pos + 14
    muts_in_window <- mutation_data %>%
      filter(position >= window_start_pos & position <= window_end_pos)

    if (nrow(muts_in_window) >= 12) {
      return(list(is_frameshift = TRUE, frameshift_start_mutation = mutation_data$mutation[i]))
    }
  }
  return(list(is_frameshift = FALSE, frameshift_start_mutation = NA_character_))
}


# --- 4. Load and Prepare the Data ---
cat("Loading and preparing data from:", input_file_path, "\n")
raw_data <- fread(input_file_path, colClasses = "character", fill = TRUE)

# Safely handle locus vs Locus capitalization
if("Locus" %in% names(raw_data) && !"locus" %in% names(raw_data)) {
  setnames(raw_data, "Locus", "locus")
}

# --- 4b. Data Prep (Imported from Heatmap Script) ---
# Add potentially missing columns for safety
for (orf in c("gag", "pro", "pol", "env", "np9", "rec")) {
  if (!orf %in% names(raw_data)) raw_data[, (orf) := NA_character_]
}

# Normalize "N" placeholder to NA
raw_data[tolower(np9) == "n" | tolower(np9) == "na" | np9 == "", np9 := NA_character_]
raw_data[tolower(rec) == "n" | tolower(rec) == "na" | rec == "",  rec  := NA_character_]

# Derive provirus_type
raw_data[, provirus_type := fcase(
  !is.na(np9) & is.na(rec), "type1",
  is.na(np9) & !is.na(rec), "type2",
  default = "unknown"
)]

# Convert to tibble for downstream dplyr use and filter to proviruses
cleaned_data <- as_tibble(raw_data) %>%
  filter(!is.na(locus) & locus != "") %>%
  filter(tolower(Structure) %like% "provirus")

deduplicated_data <- cleaned_data %>% distinct(locus, ID, Haplotype, .keep_all = TRUE)
cat("Filtered down to", nrow(deduplicated_data), "unique provirus observations.\n")


# --- 5. Identify All Disabling Mutations ---
cat("--- Starting Dual-Layer Mutation Analysis ---\n")
# Helper: drop 'Undetermined:...' masked-region markers from a comma-separated
# missense list. These are NOT mutations (they mark divergent regions too far from
# the KCON consensus to call), so they must not be counted or fed to the frameshift
# density heuristic. A token is a mask marker iff it starts with 'undetermined'.
strip_undetermined_tokens <- function(mutation_list_str) {
  vapply(mutation_list_str, function(s) {
    if (is.na(s)) return(NA_character_)
    toks <- str_split(s, ",\\s*")[[1]]
    toks <- toks[!str_starts(tolower(str_trim(toks)), "undetermined")]
    paste(toks, collapse = ", ")
  }, character(1), USE.NAMES = FALSE)
}

mutations_long <- deduplicated_data %>%
  select(locus, ID, Haplotype, provirus_type, starts_with("missense_")) %>%
  pivot_longer(
    cols = starts_with("missense_"),
    names_to = "orf",
    values_to = "mutation_list",
    names_prefix = "missense_"
  ) %>%
  filter(!is.na(mutation_list) & mutation_list != "" & mutation_list != "ComparisonError_ProteinMissing") %>%
  # Skip masked-region 'Undetermined:...' markers before any counting/classification.
  mutate(mutation_list = strip_undetermined_tokens(mutation_list)) %>%
  filter(mutation_list != "") %>%
  mutate(mutation_count = str_count(mutation_list, ",") + 1)

# Step A: Run the Backup Heuristic on highly mutated sequences
data_to_scan <- mutations_long %>% filter(mutation_count >= 12)

frameshift_results <- data_to_scan %>%
  mutate(
    frameshift_info  = map(mutation_list, detect_frameshift_sequential_heuristic),
    is_heuristic_fs  = map_lgl(frameshift_info, "is_frameshift"),
    heuristic_fs_mut = map_chr(frameshift_info, "frameshift_start_mutation")
  ) %>%
  select(locus, ID, Haplotype, orf, is_heuristic_fs, heuristic_fs_mut)

# Step B: Combine Explicit Upstream labels with the Backup Heuristic
mutations_analyzed <- mutations_long %>%
  left_join(frameshift_results, by = c("locus", "ID", "Haplotype", "orf")) %>%
  mutate(is_heuristic_fs = replace_na(is_heuristic_fs, FALSE))

all_classified_mutations <- mutations_analyzed %>%
  separate_rows(mutation_list, sep = ",\\s*") %>%
  rename(mutation = mutation_list) %>%
  filter(mutation != "" & tolower(mutation) != "na" & tolower(mutation) != "n") %>%
  mutate(
    position = parse_mutation_position(mutation),
    mutation_type = case_when(
      str_detect(mutation, fixed("*")) ~ "Nonsense",
      str_detect(tolower(mutation), "frameshift") ~ "Frameshift", # Primary: Explicit upstream tag
      is_heuristic_fs & mutation == heuristic_fs_mut ~ "Frameshift", # Backup: Heuristic trigger
      TRUE ~ "Missense"
    )
  ) %>%
  filter(!is.na(position))

disabling_mutations <- all_classified_mutations %>%
  filter(mutation_type %in% c("Nonsense", "Frameshift"))

cat("Found", nrow(disabling_mutations), "disabling mutations using combined logic.\n")


# --- 6. Process and Summarize Disabling Events ---
cat("Summarizing defects for unspliced (Gag-Pro-Pol) and spliced transcripts...\n")

# --- For Gag-Pro-Pol (Unspliced) ---
gpp_orfs <- c("gag", "pro", "pol")
gpp_primary_defects <- disabling_mutations %>%
  filter(orf %in% gpp_orfs) %>%
  mutate(orf = factor(orf, levels = gpp_orfs)) %>%
  arrange(locus, ID, Haplotype, orf, position) %>%
  group_by(locus, ID, Haplotype) %>%
  slice_head(n = 1) %>%
  ungroup()

summary_gpp <- deduplicated_data %>%
  select(locus, ID, Haplotype) %>%
  left_join(gpp_primary_defects, by = c("locus", "ID", "Haplotype")) %>%
  mutate(
    defect_category = if_else(
      is.na(mutation_type),
      "Intact (Gag-Pro-Pol)",
      str_to_title(paste(orf, mutation_type))
    )
  ) %>%
  count(locus, defect_category)

# --- For Spliced ORFs (Env, Np9/K-rev) ---
# IMPORTANT FILTER: Only look at np9 defects for Type 1, and rec defects for Type 2
spliced_orfs <- c("env", "np9", "rec")
spliced_primary_defects <- disabling_mutations %>%
  filter(orf %in% spliced_orfs) %>%
  filter(
    (orf == "env") |
      (orf == "np9" & provirus_type == "type1") |
      (orf == "rec" & provirus_type == "type2")
  ) %>%
  mutate(orf = factor(orf, levels = spliced_orfs)) %>%
  arrange(locus, ID, Haplotype, orf, position) %>%
  group_by(locus, ID, Haplotype) %>%
  slice_head(n = 1) %>%
  ungroup()

summary_spliced <- deduplicated_data %>%
  select(locus, ID, Haplotype) %>%
  left_join(spliced_primary_defects, by = c("locus", "ID", "Haplotype")) %>%
  mutate(
    defect_category = if_else(
      is.na(mutation_type),
      "Intact (Spliced)",
      str_to_title(paste(orf, mutation_type))
    )
  ) %>%
  count(locus, defect_category)


# --- 7. Generate and Save the Plot ---
cat("Generating two-panel stacked bar chart...\n")

# Order loci by total GPP defects, but include EVERY locus (loci whose proviruses are
# all intact have 0 defects and sort last). The old code filtered to defect rows first,
# so all-intact loci fell outside the factor levels and collapsed into a spurious "NA"
# row. Build the order over the union of GPP + spliced loci so both panels share it and
# neither produces an NA bar.
locus_order <- dplyr::bind_rows(summary_gpp, summary_spliced) %>%
  group_by(locus) %>%
  summarise(total_defects = sum(n[!str_starts(defect_category, "Intact")]), .groups = "drop") %>%
  arrange(desc(total_defects)) %>%
  filter(!is.na(locus) & locus != "") %>%
  pull(locus)

# Create a master color palette for consistency.
# The old palette (RColorBrewer "Paired"/"Set2") encoded defects with red-vs-green pairs,
# which fails colour-blind readers. Replace with an Okabe-Ito-derived scheme: each ORF gets
# a distinct colour-blind-safe HUE, and Nonsense vs Frameshift are distinguished by SHADE
# (strong hue = Nonsense, lightened = Frameshift) rather than by red/green. "Intact" = grey.
all_categories <- c("Intact (Gag-Pro-Pol)", "Gag Nonsense", "Gag Frameshift", "Pro Nonsense", "Pro Frameshift", "Pol Nonsense", "Pol Frameshift",
                    "Intact (Spliced)", "Env Nonsense", "Env Frameshift", "Np9 Nonsense", "Np9 Frameshift", "Rec Nonsense", "Rec Frameshift")

colors <- c(
  "Intact (Gag-Pro-Pol)" = "#BDBDBD",
  "Gag Nonsense"  = "#E69F00", "Gag Frameshift" = "#F6CB6B",   # orange / light orange
  "Pro Nonsense"  = "#0072B2", "Pro Frameshift" = "#7FB8DC",   # blue / light blue
  "Pol Nonsense"  = "#009E73", "Pol Frameshift" = "#74CDAF",   # green / light green
  "Intact (Spliced)" = "#636363",
  "Env Nonsense"  = "#D55E00", "Env Frameshift" = "#F0A06B",   # vermillion / light vermillion
  "Np9 Nonsense"  = "#56B4E9", "Np9 Frameshift" = "#A9DBF4",   # sky blue / light sky blue
  "Rec Nonsense"  = "#CC79A7", "Rec Frameshift" = "#E3B4CD"    # reddish purple / light
)
colors <- colors[all_categories]

# Plot for GPP
plot_gpp <- summary_gpp %>%
  mutate(locus = factor(locus, levels = locus_order)) %>%
  ggplot(aes(x = locus, y = n, fill = defect_category)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(name = "Defect Status", values = colors, drop = FALSE) +
  scale_x_discrete(labels = hml2_short_locus) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Unspliced (Gag-Pro-Pol)", x = NULL, y = "Proportion of Proviruses") +
  theme_pub(base_size = 20) +
  coord_flip() +
  theme(axis.text.y  = element_text(size = 17),
        axis.text.x  = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        legend.text  = element_text(size = 17),
        legend.title = element_text(size = 19))

# Plot for Spliced
plot_spliced <- summary_spliced %>%
  mutate(locus = factor(locus, levels = locus_order)) %>%
  ggplot(aes(x = locus, y = n, fill = defect_category)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(name = "Defect Status", values = colors, drop = FALSE) +
  scale_x_discrete(labels = hml2_short_locus) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Spliced (Env, Np9/Rec)", x = NULL, y = "Proportion of Proviruses") +
  theme_pub(base_size = 20) +
  coord_flip() +
  theme(axis.text.y  = element_text(size = 17),
        axis.text.x  = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        legend.text  = element_text(size = 17),
        legend.title = element_text(size = 19))

# Combine with patchwork. tag_levels = "A" stamps uniform A/B panel letters that were
# previously missing on the two sub-panels (flagged in QC).
combined_plot <- plot_gpp + plot_spliced +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.title = element_text(size = rel(1.6), hjust = 0.5),
    plot.tag = element_text(size = 30, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save the plot. ~70 loci on each y-axis -> tall canvas so the (now larger) locus labels
# have vertical room; this figure ships as a full-page standalone supplement.
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
save_fig(combined_plot, "hml2_disabling_mutations_by_transcript", width = 16, height = 19, dir = output_dir)
cat("--- Script Finished ---\n")
