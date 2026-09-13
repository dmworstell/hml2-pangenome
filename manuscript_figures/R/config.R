# =============================================================================
# 00_config.R  --  Shared configuration for the HML-2 / HERV-K figure scripts
# -----------------------------------------------------------------------------
# Source this at the top of every plotting script. Direct invocations should
# start in the repository root. HML2_CONFIG can override this configuration.
#
# Provides, in one place:
#   * theme_pub() / theme_pub_heatmap()  -- minimalist high-impact-journal theme
#   * Okabe-Ito colour-blind-safe palettes + scale_*_hml2_* convenience scales
#   * save_fig()  -- writes a PNG *and* a vector PDF (fonts embedded) at once
#   * locus-name cleanup (hml2_clean_locus / hml2_fix_7p22) and HML2_DROP_LOCI
#   * provirus-type derivation + heatmap label helper (shadowtext, with fallback)
#
# Everything degrades gracefully: optional packages (shadowtext, ragg, systemfonts)
# are used when present and silently fall back when not, so sourcing never errors.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
hml2_path <- function(name, default, must_exist = FALSE) {
  value <- Sys.getenv(name, unset = default)
  value <- path.expand(value)
  if (must_exist && !file.exists(value)) {
    stop(sprintf("%s does not exist: %s", name, value), call. = FALSE)
  }
  normalizePath(value, winslash = "/", mustWork = FALSE)
}

# All locations may be overridden without editing source code. Defaults expect
# additional inputs under data/, with generated panels under outputs/panels/
# and assembled manuscript figures under outputs/manuscript/.
HML2_REPO_ROOT       <- hml2_path("HML2_REPO_ROOT", getwd())
# The out-of-repo project resource root (.../HML2_ProjectResources). When set, HML2_DATA_ROOT
# derives from it as ${HML2_PROJECT_RESOURCES}/data -- the canonical layout keeps the data tree one
# level BELOW the project root. An explicitly set HML2_DATA_ROOT still wins; with NEITHER set,
# HML2_DATA_ROOT falls back to <repo>/data. Full R-panel inputs are not bundled.
HML2_PROJECT_RESOURCES <- Sys.getenv("HML2_PROJECT_RESOURCES", unset = "")
.hml2_data_root_default <- if (nzchar(HML2_PROJECT_RESOURCES)) {
  file.path(HML2_PROJECT_RESOURCES, "data")
} else {
  file.path(HML2_REPO_ROOT, "data")
}
HML2_DATA_ROOT       <- hml2_path("HML2_DATA_ROOT", .hml2_data_root_default)
HML2_FIG_DIR         <- hml2_path("HML2_FIG_DIR", file.path(HML2_REPO_ROOT, "outputs", "panels"))
HML2_POP_FIG_DIR     <- hml2_path("HML2_POP_FIG_DIR", file.path(HML2_FIG_DIR, "population"))
HML2_MANUSCRIPT_DIR  <- hml2_path("HML2_MANUSCRIPT_DIR", file.path(HML2_REPO_ROOT, "outputs", "manuscript"))
HML2_STATIC_FIG_DIR  <- hml2_path("HML2_STATIC_FIG_DIR", file.path(HML2_DATA_ROOT, "static_figures"))

HML2_ORF_TABLE       <- hml2_path("HML2_ORF_TABLE", file.path(HML2_DATA_ROOT, "catalog", "combined_hml2_orf_analysis.tsv"))
HML2_SAMPLE_INFO     <- hml2_path("HML2_SAMPLE_INFO", file.path(HML2_DATA_ROOT, "ref", "igsr_samples.tsv"))
HML2_CATALOG_ALIASES <- hml2_path(
  "HML2_CATALOG_ALIASES",
  file.path(HML2_REPO_ROOT, "manuscript_figures", "data", "CATALOG_LOCUS_ALIASES.tsv")
)
HML2_CATALOG_MANIFEST_DIR <- hml2_path(
  "HML2_CATALOG_MANIFEST_DIR",
  file.path(HML2_REPO_ROOT, "manuscript_figures", "data", "catalog_manifest")
)
HML2_REFERENCE_BLIND_SPOT_DIR <- hml2_path(
  "HML2_REFERENCE_BLIND_SPOT_DIR",
  file.path(HML2_REPO_ROOT, "manuscript_figures", "data", "reference_blind_spot")
)
HML2_REF_DIR         <- hml2_path("HML2_REF_DIR", file.path(HML2_DATA_ROOT, "ref"))
HML2_LTR_DIR         <- hml2_path("HML2_LTR_DIR", file.path(HML2_DATA_ROOT, "ltr_sequences"))
HML2_LRSR_DIR        <- hml2_path("HML2_LRSR_DIR", file.path(HML2_DATA_ROOT, "longread_shortread"))
HML2_FIBERSEQ_DIR    <- hml2_path("HML2_FIBERSEQ_DIR", file.path(HML2_DATA_ROOT, "fiberseq"))

# These four haplotypes are reference assemblies, not members of the 292-person
# population panel. Override with a comma-separated HML2_REFERENCE_IDS value if
# a future catalog uses different identifiers.
HML2_REFERENCE_IDS <- strsplit(
  Sys.getenv(
    "HML2_REFERENCE_IDS",
    unset = "GCA,chm13v2.0,hg002v1.1.mat,hg002v1.1.pat.PanSN"
  ),
  ",",
  fixed = TRUE
)[[1]]
HML2_REFERENCE_IDS <- trimws(HML2_REFERENCE_IDS)
HML2_REFERENCE_IDS <- HML2_REFERENCE_IDS[nzchar(HML2_REFERENCE_IDS)]

hml2_is_reference_id <- function(x) as.character(x) %in% HML2_REFERENCE_IDS

hml2_require_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop(sprintf("Missing %s: %s\nSet the corresponding HML2_* environment variable; see DATA_REQUIREMENTS.md.",
                 label, path), call. = FALSE)
  }
  invisible(path)
}

# -----------------------------------------------------------------------------
# Fonts -- pick a clean sans that actually exists on this machine, so neither
# the raster nor the PDF device throws "font family not found" warnings.
# -----------------------------------------------------------------------------
.hml2_pick_font <- function() {
  candidates <- c("Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans")
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    fam <- unique(systemfonts::system_fonts()$family)
    hit <- candidates[candidates %in% fam]
    if (length(hit)) return(hit[[1]])
  }
  "sans"
}
HML2_FONT <- .hml2_pick_font()

# -----------------------------------------------------------------------------
# Colour-blind-safe palettes (Okabe & Ito, 2008)
# -----------------------------------------------------------------------------
okabe_ito <- c(
  black        = "#000000",
  orange       = "#E69F00",
  sky_blue     = "#56B4E9",
  bluish_green = "#009E73",
  yellow       = "#F0E442",
  blue         = "#0072B2",
  vermillion   = "#D55E00",
  reddish_purple = "#CC79A7"
)

# HML-2 internal type: type1 carries np9, type2 carries rec. Aliases included so
# the same vector works whether a script labels them type1/Type I/TypeI etc.
HML2_TYPE_COLORS <- c(
  type1 = "#0072B2", type2 = "#D55E00",
  TypeI = "#0072B2", TypeII = "#D55E00",
  "Type I" = "#0072B2", "Type II" = "#D55E00",
  "Type 1" = "#0072B2", "Type 2" = "#D55E00"
)

# 1000 Genomes super-populations (avoids the low-contrast Okabe-Ito yellow).
HML2_POP_COLORS <- c(
  AFR = "#E69F00",  # orange
  AMR = "#CC79A7",  # reddish purple
  EAS = "#009E73",  # bluish green
  EUR = "#0072B2",  # blue
  SAS = "#D55E00"   # vermillion
)

# Structural classes seen in the ORF table.
HML2_STRUCT_COLORS <- c(
  "Provirus"            = "#0072B2",
  "Provirus_from_Multi" = "#56B4E9",
  "Multi-copy"          = "#D55E00",
  "Solo-LTR"            = "#E69F00",
  "Solo_LTR"            = "#E69F00",
  "SOLO_LTR"            = "#E69F00",
  "Fragment"            = "#999999",
  "Absent"              = "#DDDDDD"
)

# Sequential / diverging continuous fills (both colour-blind-safe).
HML2_SEQ_OPTION  <- "viridis"   # for viridis::scale_*_viridis_c(option = HML2_SEQ_OPTION)
HML2_DIVERGING   <- c("#0072B2", "#F7F7F7", "#D55E00")  # blue-white-vermillion (saturated)
# Lighter diverging -- use when tiles carry BLACK text, so the dark ends don't
# swallow the labels (legibility over punch).
HML2_DIVERGING_LIGHT <- c("#74A9CF", "#FFFFFF", "#FDAE6B")  # light blue-white-light orange

# Population SHAPES (use instead of colour when a plot would otherwise be a
# 5-colour rainbow, e.g. the per-locus frequency grid). Distinct filled glyphs.
HML2_POP_SHAPES <- c(AFR = 16, AMR = 17, EAS = 15, EUR = 18, SAS = 8)
scale_shape_hml2_pop <- function(...) scale_shape_manual(values = HML2_POP_SHAPES, ...)

scale_fill_hml2_type   <- function(...) scale_fill_manual(values = HML2_TYPE_COLORS, ...)
scale_colour_hml2_type <- function(...) scale_colour_manual(values = HML2_TYPE_COLORS, ...)
scale_color_hml2_type  <- scale_colour_hml2_type
scale_fill_hml2_pop    <- function(...) scale_fill_manual(values = HML2_POP_COLORS, ...)
scale_colour_hml2_pop  <- function(...) scale_colour_manual(values = HML2_POP_COLORS, ...)
scale_color_hml2_pop   <- scale_colour_hml2_pop
scale_fill_hml2_struct <- function(...) scale_fill_manual(values = HML2_STRUCT_COLORS, ...)

# -----------------------------------------------------------------------------
# Theme -- the look prestige journals actually print: data-dense and plain, NOT
# decorative. A SHORT, LARGE, NON-BOLD title (journals dislike bold anywhere on
# a figure), large black axis labels, no subtitle (that belongs in the caption),
# both faint gridlines kept (reads data-rich rather than airy), tight margins.
# `titles = FALSE` drops the title for panels that are sub-parts of a montage.
# `grid` chooses which major gridlines survive.
# -----------------------------------------------------------------------------
theme_pub <- function(base_size = 16, base_family = HML2_FONT,
                      grid = c("both", "y", "x", "none"),
                      legend = "right", titles = FALSE) {
  grid <- match.arg(grid)
  th <- theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text             = element_text(colour = "black"),
      axis.title       = element_text(colour = "black", size = rel(1.15)),
      axis.text        = element_text(colour = "grey20", size = rel(1.0)),
      axis.line        = element_line(colour = "grey30", linewidth = 0.4),
      axis.ticks       = element_line(colour = "grey30", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      legend.position  = legend,
      legend.title     = element_text(colour = "black", size = rel(1.0)),
      legend.text      = element_text(size = rel(0.95)),
      legend.key.size  = unit(0.9, "lines"),
      strip.text       = element_text(colour = "black", size = rel(1.05)),
      strip.background = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(8, 8, 6, 6)
    )
  if (titles) {
    # short + large + plain (never bold).
    th <- th + theme(
      plot.title    = element_text(face = "plain", size = rel(2.0), hjust = 0.5,
                                   margin = margin(b = 9)),
      plot.subtitle = element_blank()
    )
  } else {
    th <- th + theme(plot.title = element_blank(), plot.subtitle = element_blank())
  }
  th <- th + theme(plot.caption = element_blank())
  if (grid == "none") th <- th + theme(panel.grid.major = element_blank())
  if (grid == "x")    th <- th + theme(panel.grid.major.y = element_blank())
  if (grid == "y")    th <- th + theme(panel.grid.major.x = element_blank())
  th
}

# Heatmap variant: no grid, no axis lines/ticks, tight, large tick labels.
theme_pub_heatmap <- function(base_size = 16, base_family = HML2_FONT, titles = FALSE) {
  theme_pub(base_size, base_family, grid = "none", titles = titles) +
    theme(
      panel.grid.major = element_blank(),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      panel.border     = element_blank()
    )
}

# Make every plot use the house style by default (scripts can still override).
theme_set(theme_pub())

# -----------------------------------------------------------------------------
# save_fig() -- one call writes a 300-dpi PNG *and* a vector PDF (cairo, so the
# chosen font is embedded). `name` may omit the extension. Returns the paths.
# -----------------------------------------------------------------------------
save_fig <- function(plot, name, width = 8, height = 6, dpi = 300,
                     dir = HML2_FIG_DIR, formats = c("png", "pdf"), bg = "white") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  name  <- sub("\\.(png|pdf)$", "", name)
  paths <- character(0)
  if ("png" %in% formats) {
    p   <- file.path(dir, paste0(name, ".png"))
    dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
    ggsave(p, plot = plot, width = width, height = height, dpi = dpi, bg = bg, device = dev)
    paths <- c(paths, p)
  }
  if ("pdf" %in% formats) {
    p   <- file.path(dir, paste0(name, ".pdf"))
    dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
    ggsave(p, plot = plot, width = width, height = height, device = dev, bg = bg)
    paths <- c(paths, p)
  }
  message("  saved: ", paste(basename(paths), collapse = " + "), "  ->  ", dir)
  invisible(paths)
}

# -----------------------------------------------------------------------------
# Locus-name cleanup
#   * 7p22.1b was a spurious sub-label for what is a single 7p22.1 locus.
#   * HML2_DROP_LOCI: any loci to exclude wholesale (none at present).
# -----------------------------------------------------------------------------
HML2_DROP_LOCI <- character(0)

hml2_fix_7p22 <- function(x) {
  x <- as.character(x)
  gsub("7p22\\.1[ab]", "7p22.1", x, perl = TRUE)
}

hml2_clean_locus <- function(x) {
  hml2_fix_7p22(x)
}

# Strip the "HML-2_" prefix (and any "_hg38" suffix) for axis / tile labels --
# journal figures don't repeat the family name on every tick. Use as a ggplot
# scale labeller, e.g. scale_x_discrete(labels = hml2_short_locus).
hml2_short_locus <- function(x) {
  x <- hml2_clean_locus(x)
  x <- gsub("^HML-?2[_-]", "", x, perl = TRUE)
  gsub("_hg38$", "", x, perl = TRUE)
}

# -----------------------------------------------------------------------------
# Structural-state vocabulary shared by the catalog-wide frequency and
# reference-blind-spot analyses. These labels describe observed catalog rows;
# a missing locus x haplotype row is represented separately as NA and must never
# be passed through this mapper as an inferred Absent call.
# -----------------------------------------------------------------------------
HML2_STRUCTURAL_STATES <- c("Absent", "Solo-LTR", "Fragment", "Other", "Provirus", "Multi-copy")
HML2_STRUCTURAL_PRIORITY <- c("Multi-copy", "Provirus", "Solo-LTR", "Fragment", "Other", "Absent")

hml2_structural_state <- function(x) {
  x <- as.character(x)
  text <- ifelse(is.na(x), "", x)
  state <- rep("Other", length(text))
  state[grepl("Absent", text, ignore.case = TRUE)] <- "Absent"
  state[grepl("Solo[-_]LTR", text, ignore.case = TRUE)] <- "Solo-LTR"
  state[grepl("Fragment", text, ignore.case = TRUE)] <- "Fragment"
  state[grepl("Provirus", text, ignore.case = TRUE)] <- "Provirus"
  state[grepl("Provirus_from_Multi|_part", text)] <- "Multi-copy"
  state
}

hml2_structural_rank <- function(x) match(as.character(x), HML2_STRUCTURAL_PRIORITY)

# -----------------------------------------------------------------------------
# provirus_type derivation -- type1 carries np9, type2 carries rec. Treats
# "", NA, "NA", "Absent", "None" as absent. Vectorised.
# -----------------------------------------------------------------------------
hml2_is_present <- function(v) {
  # The ORF table encodes "absent" several ways, incl. a bare "n".
  !is.na(v) & !(toupper(trimws(as.character(v))) %in% c("", "N", "NA", "ABSENT", "NONE", "."))
}

hml2_derive_provirus_type <- function(np9, rec) {
  ifelse(hml2_is_present(rec), "type2",
  ifelse(hml2_is_present(np9), "type1", NA_character_))
}

# -----------------------------------------------------------------------------
# QC: in the 2026-06 build, 59 graph-derived samples FAILED genotyping (≈all
# Insertion_Absent, ~0 proviruses, vs a real median of ~85). They form a phantom
# cluster that wrecks genotype PCAs / frequency analyses. Return the IDs with a
# plausible provirus count so callers can drop the failures.
# `dt` needs columns ID + Structure.
# -----------------------------------------------------------------------------
hml2_well_genotyped_ids <- function(dt, min_prov = 40) {
  dt <- data.table::as.data.table(dt)
  pc <- dt[, .(np = sum(grepl("^Provirus", Structure))), by = ID]
  pc[np >= min_prov]$ID
}

# -----------------------------------------------------------------------------
# Wilson score confidence interval for a binomial proportion (x of n). Vectorised,
# no extra package. Use for error bars on frequencies / proportions:
#   mutate(lo = hml2_wilson_lo(k, n), hi = hml2_wilson_hi(k, n))
# -----------------------------------------------------------------------------
.hml2_wilson <- function(x, n, conf, side) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- ifelse(n > 0, x / n, NA_real_)
  d <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  if (side == "lo") pmax(0, centre - half) else pmin(1, centre + half)
}
hml2_wilson_lo <- function(x, n, conf = 0.95) .hml2_wilson(x, n, conf, "lo")
hml2_wilson_hi <- function(x, n, conf = 0.95) .hml2_wilson(x, n, conf, "hi")

# -----------------------------------------------------------------------------
# Heatmap tile labels -- plain BLACK text, no halo (white-haloed text reads badly
# in print). Keep it legible by drawing on LIGHT / translucent tiles: pair this
# with scale_fill_viridis(alpha = ...) or the light diverging palette below.
# -----------------------------------------------------------------------------
hml2_tile_label <- function(..., colour = "black", size = 3.4) {
  geom_text(..., colour = colour, size = size)
}

message("[00_config.R] loaded  |  font: ", HML2_FONT,
        "  |  shadowtext: ", requireNamespace("shadowtext", quietly = TRUE),
        "  |  figures -> ", HML2_FIG_DIR)
