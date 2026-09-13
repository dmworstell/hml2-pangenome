source(Sys.getenv("HML2_CONFIG", file.path("R", "config.R")))
# Fusion ORFs — MODEL 1: WITH the canonical HML-2 -1 ribosomal frameshifts.
# gag(fr0) -(-1 FS)-> pro(fr2) -(-1 FS)-> pol(fr1) is the natural Gag-Pro-Pol
# polyprotein, so a provirus can translate it whenever gag, pro AND pol are intact
# (orf_analysis already scores each ORF in its own frame, indels included).
# Classes (with frameshifts): Gag / Gag-Pro / Gag-Pro-Pol / all-ORFs(+env).
# env is normally spliced & out of frame with pol, so a gag-pro-pol-ENV megaprotein
# additionally needs env intact (read-through is examined in the Model-2 scan).
library(tidyverse); library(data.table)
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else
  HML2_ORF_TABLE
output_dir <- HML2_FIG_DIR
d <- fread(input_file, colClasses = "character"); if ("Locus" %in% names(d)) setnames(d, "Locus", "locus")
ok <- function(x) x %in% c("Intact", "Intact_FS_End")
# Undetermined = too divergent from KCON to call: neither an intact ORF nor a disrupted one.
ud <- function(x) tolower(ifelse(is.na(x), "", x)) == "undetermined"
prov <- as_tibble(d) %>% filter(str_detect(Structure, "Provirus")) %>%
  mutate(gi = ok(gag), pri = ok(pro), pli = ok(pol), ei = ok(env),
         uncallable = ud(gag) & ud(pro) & ud(pol) & ud(env),
         fusion = case_when(
           uncallable          ~ "undetermined",
           gi & pri & pli & ei ~ "Gag-Pro-Pol-Env (all ORFs)",
           gi & pri & pli      ~ "Gag-Pro-Pol",
           gi & pri            ~ "Gag-Pro",
           gi                  ~ "Gag only",
           TRUE                ~ "none/disrupted"),
         fusion = factor(fusion, levels = c("undetermined","none/disrupted","Gag only","Gag-Pro",
                                            "Gag-Pro-Pol","Gag-Pro-Pol-Env (all ORFs)")),
         locus = hml2_short_locus(locus))
cat("Overall fusion potential WITH canonical frameshifts (per provirus):\n")
print(prov %>% count(fusion))
cat("\nLoci where the full Gag-Pro-Pol-Env (all-ORF) product is possible in >=1 provirus:\n")
print(prov %>% filter(fusion == "Gag-Pro-Pol-Env (all ORFs)") %>%
        count(locus, name = "n_proviruses", sort = TRUE), n = 30)
cat("\n1q22 breakdown:\n"); print(prov %>% filter(locus == "1q22") %>% count(fusion))

# per-locus stacked fraction
pl <- prov %>% count(locus, fusion) %>% group_by(locus) %>% mutate(frac = n/sum(n)) %>% ungroup() %>%
  group_by(locus) %>% mutate(best = max(as.integer(fusion))) %>% ungroup() %>%
  mutate(locus = fct_reorder(locus, best))
pal <- c("none/disrupted"="grey85","Gag only"="#F0E442","Gag-Pro"="#56B4E9",
         "Gag-Pro-Pol"="#0072B2","Gag-Pro-Pol-Env (all ORFs)"="#D55E00")
p <- ggplot(pl, aes(frac, locus, fill = fusion)) + geom_col(width=0.85) +
  scale_x_continuous(labels = scales::percent, expand=c(0,0)) +
  scale_fill_manual(values = pal, name = "Longest poly-ORF\n(with frameshifts)") +
  labs(title = "HML-2 polyprotein / fusion-ORF potential (canonical -1 frameshifts)",
       x = "Fraction of proviruses", y = NULL) +
  theme_pub(base_size = 13) + theme(axis.text.y = element_text(size=8))
save_fig(p, "hml2_fusion_orf_potential", width = 9, height = 11, dir = output_dir)
write_tsv(prov %>% select(locus, ID, Haplotype, ID_Full, gag, pro, pol, env, fusion),
          file.path(output_dir, "hml2_fusion_orf_calls.tsv"))
cat("\n--- done (Model 1: with frameshifts) ---\n")
