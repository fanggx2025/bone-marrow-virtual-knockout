## =============================================================================
## s08_fig6_consensus.R -- Fig6 three-way consensus core figure
## Unified journal style (Scientific Reports): Arial/sans, white bg, 300 dpi,
## bold lowercase panel letters, no in-figure overall title, English only.
## (reads tier/direction values from the English-translated CSVs)
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(grid); library(gridExtra); library(VennDiagram)
})
ROOT <- "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa"
TBL  <- file.path(ROOT, "results/tables")
FIG  <- file.path(ROOT, "results/figures")

sc   <- read.csv(file.path(TBL, "s04_knk_scores_full.csv"), stringsAsFactors = FALSE)
s06  <- read.csv(file.path(TBL, "s06_consensus_targets.csv"),  stringsAsFactors = FALSE)
s05c <- read.csv(file.path(TBL, "s05_compound_targets.csv"),   stringsAsFactors = FALSE)

KO_set   <- sc %>% filter(direction == "anti-osteogenic" & p_perturb < 0.05) %>% pull(KO) %>% as.character()
## Traj_set: genes whose expression decreases along osteogenic pseudotime (rho<0)
Traj_set <- s06 %>% filter(rho_pseudotime < 0) %>% pull(gene) %>% as.character()
Comp_set <- unique(s05c$gene_symbol) %>% as.character()

n12  <- length(intersect(KO_set, Comp_set))
n13  <- length(intersect(KO_set, Traj_set))
n23  <- length(intersect(Comp_set, Traj_set))
n123 <- length(intersect(KO_set, intersect(Comp_set, Traj_set)))
cat(sprintf("KO=%d Traj=%d Comp=%d | 12=%d 13=%d 23=%d 123=%d\n",
            length(KO_set), length(Traj_set), length(Comp_set), n12, n13, n23, n123))

## ---- Panel a: triple venn --------------------------------------------------
venn <- tryCatch({
  VennDiagram::venn.diagram(
    x = list(`Virtual-KO\ndruggable` = KO_set,
             `Compound\ntargets`      = Comp_set,
             `Osteogenic\ntraj down`  = Traj_set),
    filename = NULL,
    fill = c("#E45756", "#4C78A8", "#72B7B2"),
    cat.cex = 0.95, cex = 0.95, margin = 0.04)
}, error = function(e) NULL)
if (!is.null(venn)) venn <- grid::gTree(children = venn)

## ---- Panel b: tier bar -----------------------------------------------------
tier_short <- c("Tier1 three-way consensus (core targets)" = "Tier1\nconsensus\n(core)",
                "Tier2 KO+compound (trajectory not passed)" = "Tier2\nKO + drug",
                "Tier3 KO+trajectory (mechanistic)" = "Tier3\nKO + traj.",
                "Tier4 KO-only" = "Tier4\nKO only")
tab_tier <- s06 %>% count(tier) %>% mutate(tier_s = tier_short[tier])
tab_tier$tier_s <- factor(tab_tier$tier_s, levels = unname(tier_short[unique(s06$tier)]))
pB <- ggplot(tab_tier, aes(x = tier_s, y = n, fill = tier_s)) +
  geom_bar(stat = "identity", width = 0.72) +
  scale_fill_manual(values = c("#D62728", "#F58518", "#4C78A8", "#9ecae1"),
                    guide = "none") +
  geom_text(aes(label = n), vjust = -0.4, size = 2.8) +
  labs(title = "Candidate target tiers", y = "n genes", x = NULL) +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
        axis.text.x = element_text(size = 7.5),
        axis.title.y = element_text(size = 8),
        panel.grid.minor = element_blank())

## ---- Tier-1 core targets: exported as Supplementary Table S9 (s06_tier1_core.csv) ----
## (no longer drawn inside Fig. 6; presented as a standalone table)
t1 <- s06 %>% filter(grepl("Tier1", tier)) %>%
  select(gene, herbs, compounds, best_nM, rho_pseudotime, med_w_end)

if (is.null(venn)) venn <- textGrob("venn unavailable", gp = gpar(fontsize = 10))

## ---- panel letters ---------------------------------------------------------
lab <- function(s) textGrob(s, x = 0.02, y = 0.99, hjust = 0, vjust = 1,
                            gp = gpar(fontsize = 10, fontface = "bold"))
pglue <- function(l, g) arrangeGrob(l, g, ncol = 1,
                                    heights = unit.c(unit(0.16, "in"), unit(1, "null")))

fig <- grid.arrange(pglue(lab("a"), venn), pglue(lab("b"), pB),
                    ncol = 2, widths = c(1.15, 1.55))
ggsave(file.path(FIG, "Fig6_consensus.png"), fig, width = 7.2, height = 3.3,
       dpi = 300, bg = "white")
ggsave(sub("\\.png$", ".pdf", file.path(FIG, "Fig6_consensus.png")), fig, width = 7.2, height = 3.3, bg = "white")
cat("OK Fig6_consensus.png\n")
write.csv(t1, file.path(TBL, "s06_tier1_core.csv"), row.names = FALSE)
cat("Tier1:", paste(t1$gene, collapse = ", "), "\n")
