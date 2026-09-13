source(Sys.getenv("HML2_CONFIG", file.path("manuscript_figures", "R", "config.R")))
# >>> Generates: Figure 3A (hml2_cnv_distribution), S15D (hml2_signature_sharing_heatmap) <<<

# HML-2 Copy Number Variation & Recombination Analyzer
#
# This script analyzes HML-2 copy number variation (CNV) per haplotype,
# detects inter-chromosomal mapping (evidence of structural recombination),
# and tracks the sharing of massive mutation signatures between different loci.

# --- 1. Load Necessary Libraries ---
# Ensure installed: install.packages(c("tidyverse", "data.table", "viridis", "patchwork", "ggrepel"))
library(tidyverse)
library(data.table)
library(viridis)
library(patchwork)
library(ggrepel)

# --- 2. User-Defined Parameters ---
input_file_path <- HML2_ORF_TABLE
output_dir <- HML2_FIG_DIR

# --- 3. Load and Prepare the Data ---
cat("Loading and preparing data from:", input_file_path, "\n")
raw_data <- fread(input_file_path, colClasses = "character", fill = TRUE)

if("Locus" %in% names(raw_data) && !"locus" %in% names(raw_data)) {
  setnames(raw_data, "Locus", "locus")
}

cleaned_data <- as_tibble(raw_data) %>%
  filter(!is.na(locus) & locus != "")

# Collapse arm-confirmed (MATCH) PRESENT acrocentric per-arm rows into their type family, so
# every present acrocentric element counts under acro_type1/acro_type2 in the CNV (Sec 4) and
# mega-signature (Sec 6) analyses -- otherwise the 81 arm-confirmed rows still share signatures
# with the family rows and spuriously read as inter-locus recombination (e.g. 21p13 <-> acro_type2,
# which are the SAME satellite family). Absent per-arm rows are left as-is (they must not inflate
# family copy_number). The per-arm labels are preserved in orig_Locus and drive the S6A scatter (Sec 5).
acro_fam_map <- c("HML-2_13p13" = "HML-2_acro_type1", "HML-2_15p13b" = "HML-2_acro_type1",
                  "HML-2_15p13a" = "HML-2_acro_type2", "HML-2_21p13" = "HML-2_acro_type2",
                  "HML-2_22p13"  = "HML-2_acro_type2")
cleaned_data <- cleaned_data %>%
  mutate(locus = if_else(locus %in% names(acro_fam_map) &
                           str_detect(Structure, "^(Provirus|Solo-LTR|Fragment)"),
                         unname(acro_fam_map[locus]), locus))

cat("Loaded", nrow(cleaned_data), "total observations.\n")


# --- 4. Copy Number Variation (CNV) Analysis ---
cat("\n--- Analyzing Copy Number Variation ---\n")

# Calculate copy number per locus per haplotype
cnv_data <- cleaned_data %>%
  group_by(locus, ID, Haplotype) %>%
  summarise(copy_number = n(), .groups = "drop")

# Identify loci that actually exhibit CNV (max copy number > 1). The main panel
# shows the 20 strongest loci; a 50-60-locus x-axis is not legible at journal width.
N_CNV_MAIN <- 20
cnv_loci_all <- cnv_data %>%
  group_by(locus) %>%
  summarise(
    max_cn = max(copy_number),
    mean_cn = mean(copy_number),
    variance = var(copy_number)
  ) %>%
  filter(max_cn > 1) %>%
  arrange(desc(max_cn), desc(mean_cn))
cnv_loci <- cnv_loci_all %>% slice_head(n = N_CNV_MAIN)

cat("Found", nrow(cnv_loci_all), "CNV loci; showing the top", nrow(cnv_loci),
    "in the main panel.\n")

# Plot 1: CNV Distribution
plot_cnv <- cnv_data %>%
  filter(locus %in% cnv_loci$locus) %>%
  mutate(locus = factor(locus, levels = cnv_loci$locus)) %>%
  ggplot(aes(x = locus, y = copy_number)) +
  geom_boxplot(outlier.shape = NA, fill = "grey88", color = "grey40") +
  geom_jitter(width = 0.2, height = 0.1, alpha = 0.4, color = okabe_ito[["blue"]], size = 1.5) +
  scale_x_discrete(labels = hml2_short_locus) +
  scale_y_continuous(breaks = seq(0, max(cnv_data$copy_number), by = 1)) +
  theme_pub(base_size = 22) +
  labs(
    x = "HML-2 locus",
    y = "Copy number"
  ) +
  # Vertical labels remain unambiguous and leave room for 20 named loci.
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 19, face = "plain"),
    axis.text.y = element_text(size = 22),
    axis.title  = element_text(size = 24),
    panel.grid.minor.y = element_blank()
  )


# --- 5. Acrocentric Inter-chromosomal Mapping Analysis ---
cat("\n--- Analyzing Acrocentric Cross-Mapping ---\n")

acrocentric_loci <- c("HML-2_13p13", "HML-2_14p13", "HML-2_15p13a", "HML-2_15p13b", "HML-2_21p13", "HML-2_22p13")

# Load the NCBI assembly contig -> chromosome map (CM/JB/JA accessions -> chrN).
# Applied ONLY at plot time here; we do NOT rewrite the source TSV (cf. fast_replace.py).
# Two maps: the supplied HPRC map + an NCBI-eutils extension (Tree_Plot/contig
# resolution) covering CM chromosome accessions the first map missed.
ref_dir <- HML2_REF_DIR
contig_map <- do.call(c, lapply(
  file.path(ref_dir, c("hprc_contig_to_chr.txt", "contig_to_chr_ncbi_extra.txt")),
  function(f) if (file.exists(f)) {
    m <- fread(f, header = FALSE, col.names = c("contig", "chr"), colClasses = "character")
    setNames(m$chr, m$contig)
  } else character(0)))
contig_map <- contig_map[!duplicated(names(contig_map))]
if (!length(contig_map)) warning("no contig->chr map found; falling back to raw chr tokens.")

# Extract the mapped chromosome from the Source_Identifier (PanSN sample#hap#contig:coords).
# The contig token (last '#' segment before ':') is an NCBI accession for many assemblies;
# resolve it through the NCBI map, else fall back to any literal chrN already in the string.
mapping_data <- cleaned_data %>%
  # Only PRESENT elements have a real mapping; drop absence placeholders
  # (Placeholder / Validator_Script / Absent / Insertion_Absent) -- they are not
  # "Unknown chromosomes", they're sites where the provirus simply isn't there.
  # Key off orig_Locus (the pre-re-annotation per-arm label): S6A deliberately shows the
  # per-arm -> chromosome SCATTER, which is the evidence that the arm labels are unreliable
  # and motivated collapsing them into the acro_type1/acro_type2 families (Section 4/6 use
  # the family Locus; only this mapping panel keeps the per-arm view).
  filter(orig_Locus %in% acrocentric_loci,
         str_detect(Structure, "^(Provirus|Solo-LTR|Fragment)")) %>%
  mutate(
    contig_token = sub(":.*", "", sub(".*#", "", Source_Identifier)),
    ncbi_chr     = unname(contig_map[contig_token]),                 # NA if not an accession
    raw_chr      = str_match(Source_Identifier, "chr([A-Za-z0-9_]+):")[, 2],
    mapped_chr_clean = case_when(
      # chrUn and unplaced scaffolds are the SAME thing (sequence not assigned to a
      # chromosome), so collapse them into one "Unplaced" category.
      !is.na(ncbi_chr) & ncbi_chr == "chrUn"      ~ "Unplaced",
      !is.na(ncbi_chr)                            ~ ncbi_chr,        # resolved via NCBI map
      !is.na(raw_chr) & str_starts(raw_chr, "Un") ~ "Unplaced",
      !is.na(raw_chr) & !str_detect(raw_chr, "_") ~ paste0("chr", raw_chr),
      TRUE                                        ~ "Unplaced"
    )
  )

mapping_summary <- mapping_data %>%
  count(orig_Locus, mapped_chr_clean) %>%
  group_by(orig_Locus) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Distinguishable, colour-blind-safe chromosome palette. The previous unname(okabe_ito)
# led with pure BLACK (chr13) and re-used near-confusable hues; here every chromosome
# gets a distinct strong Okabe-Ito hue (no pure black, no low-contrast yellow), and the
# "Unplaced" catch-all is a neutral grey so it never competes with a real chromosome.
chrom_palette <- c(
  "chr4"     = okabe_ito[["orange"]],        # #E69F00
  "chr13"    = okabe_ito[["sky_blue"]],      # #56B4E9
  "chr14"    = okabe_ito[["bluish_green"]],  # #009E73
  "chr15"    = okabe_ito[["blue"]],          # #0072B2
  "chr21"    = okabe_ito[["vermillion"]],    # #D55E00
  "chr22"    = okabe_ito[["reddish_purple"]],# #CC79A7
  "Unplaced" = "#999999"                     # neutral grey for the catch-all
)

# Plot 2: Cross-Mapping Bar Chart
plot_mapping <- mapping_summary %>%
  ggplot(aes(x = orig_Locus, y = prop, fill = mapped_chr_clean)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(name = "Mapped to", values = chrom_palette) +
    scale_x_discrete(labels = hml2_short_locus) +
  theme_pub(base_size = 18) +
  labs(
    x = "Per-assembly acrocentric arm label",
    y = "Proportion of alleles"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "plain", size = 18),
    axis.text.y = element_text(size = 18),
    axis.title  = element_text(size = 20),
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 18),
    legend.position = "right"
  )


# --- 6. Mega-Signature Sharing (Evidence of Genetic Recombination) ---
cat("\n--- Analyzing Mega-Signature Sharing ---\n")

# Extract all missense mutations for Gag, Pro, Pol, Env
signatures_long <- cleaned_data %>%
  select(locus, ID_Full, starts_with("missense_")) %>%
  pivot_longer(
    cols = c(missense_gag, missense_pro, missense_pol, missense_env),
    names_to = "orf",
    values_to = "signature",
    names_prefix = "missense_"
  ) %>%
  filter(!is.na(signature) & signature != "" & signature != "ComparisonError_ProteinMissing") %>%
  # Undetermined:START-END tokens mark regions masked by the KCON-misalignment confidence gate:
  # they are NOT mutations. Strip them before counting or exact-matching signatures.
  mutate(signature = gsub("(^,+)|(,+$)", "", gsub(",{2,}", ",", gsub("Undetermined:[^,]*", "", signature)))) %>%
  filter(signature != "") %>%
  # Count mutations to find "Mega-Signatures"
  mutate(mutation_count = str_count(signature, ",") + 1) %>%
  # Filter for massive signatures (>15 mutations) to avoid background noise matches
  filter(mutation_count > 15)

# Find exact signatures that are shared by MULTIPLE DIFFERENT loci
shared_signatures <- signatures_long %>%
  group_by(orf, signature) %>%
  summarise(
    unique_loci_count = n_distinct(locus),
    shared_by = paste(unique(locus), collapse = " <-> "),
    total_occurrences = n(),
    .groups = "drop"
  ) %>%
  filter(unique_loci_count > 1) %>%
  arrange(desc(unique_loci_count), desc(total_occurrences))

cat("Found", nrow(shared_signatures), "unique Mega-Signatures shared across multiple loci.\n")
if(nrow(shared_signatures) > 0) {
  print(shared_signatures %>% select(orf, shared_by, total_occurrences) %>% head(10))
}

# Create a sharing matrix for the heatmap
if(nrow(shared_signatures) > 0) {

  # Get all occurrences of these shared signatures
  sharing_edges <- signatures_long %>%
    inner_join(shared_signatures %>% select(orf, signature), by = c("orf", "signature")) %>%
    select(locus, orf, signature) %>%
    distinct()

  # Create a pairwise co-occurrence matrix of loci sharing signatures
  pairwise_sharing <- sharing_edges %>%
    full_join(sharing_edges, by = c("orf", "signature"), relationship = "many-to-many") %>%
    filter(locus.x != locus.y) %>%
    # Sort pairs alphabetically so A-B is the same as B-A
    mutate(
      loc1 = pmin(locus.x, locus.y),
      loc2 = pmax(locus.x, locus.y)
    ) %>%
    count(loc1, loc2, name = "shared_signatures_count")

  # Plot 3: Signature Sharing Heatmap.
  # Fig 6 panel D sits full-width UNDER the trees in the composite but is only ~13x13
  # tiles, so it is downscaled hard -> push every font well up (cell counts, axis tick
  # locus labels, axis titles, legend) so nothing is tiny on the printed page.
  plot_sharing <- pairwise_sharing %>%
    ggplot(aes(x = loc1, y = loc2, fill = shared_signatures_count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = shared_signatures_count), color = "white", fontface = "plain", size = 9) +
    scale_fill_viridis(name = "Shared\nsignatures", option = HML2_SEQ_OPTION, direction = -1) +
    scale_x_discrete(labels = hml2_short_locus) +
    scale_y_discrete(labels = hml2_short_locus) +
    theme_pub_heatmap(base_size = 26) +
    labs(
      x = "HML-2 locus", y = "HML-2 locus"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "plain", size = 24),
      axis.text.y = element_text(face = "plain", size = 24),
      axis.title  = element_text(size = 28),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 24),
      legend.key.size = unit(1.4, "lines")
    )
} else {
  # Empty placeholder if no sharing is found
  plot_sharing <- ggplot() +
    annotate("text", x=0.5, y=0.5, label="No shared Mega-Signatures detected.") +
    theme_void()
}

# --- 7. Save Plots ---
cat("\nSaving plots to:", output_dir, "\n")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_fig(plot_cnv,     "hml2_cnv_distribution",       width = 13, height = 9, dir = output_dir)
save_fig(plot_mapping, "hml2_acrocentric_mapping",    width = 10, height = 7, dir = output_dir)
if(nrow(shared_signatures) > 0) {
  save_fig(plot_sharing, "hml2_signature_sharing_heatmap", width = 10, height = 8, dir = output_dir)
}

cat("--- Script Finished ---\n")
