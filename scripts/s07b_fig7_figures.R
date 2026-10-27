## =============================================================================
## s07b_fig7_figures.R -- Fig7 OA external validation, unified journal style
## (Scientific Reports: Arial/sans, white bg, 300 dpi, English, panel letters)
##
## Reuses the ALREADY-COMPUTED differential results in
##   results/tables/s07_oa_validation.csv   (roles translated to English)
## and reloads s02_stromal.rds ONLY for panel-a cell counts and panel-c
## per-cell expression -- no FindMarkers rerun.
## =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(tidyr)
})
ROOT <- "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa"
POP  <- c("CAR_cell", "Stromal_fibro", "Osteoblast")
C_NORM <- "#4C78A8"; C_OA <- "#E45756"

tab <- read.csv(file.path(ROOT, "results/tables/s07_oa_validation.csv"),
                stringsAsFactors = FALSE)

## ---- load s02 (metadata + data layer only usage) ---------------------------
obj <- readRDS(file.path(ROOT, "results/objects/s02_stromal.rds"))
obj$group <- ifelse(grepl("oa|arthr|osteo", tolower(as.character(obj$disease))),
                    "OA", "Normal")
sub <- subset(obj, cells = colnames(obj)[obj$subtype %in% POP &
                                         obj$group %in% c("OA", "Normal")])
cat("cells:", ncol(sub), " (OA", sum(sub$group == "OA"),
    "/ Normal", sum(sub$group == "Normal"), ")\n", sep = "")

## ---- Panel a: cell counts ---------------------------------------------------
cnt <- as.data.frame(table(sub$subtype, sub$group))
names(cnt) <- c("subtype", "group", "n")
cnt$subtype <- factor(cnt$subtype, levels = POP)
pA <- ggplot(cnt, aes(x = subtype, y = n, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.72) +
  scale_fill_manual(values = c(Normal = C_NORM, OA = C_OA), name = "Group") +
  labs(title = "a  Osteogenic-lineage subpopulations: cell counts, OA vs Normal",
       y = "Number of cells", x = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9.5, hjust = 0),
        legend.position = "top",
        axis.text.x = element_text(angle = 15, hjust = 1))

## ---- Panel b: log2FC heatmap (from computed table) --------------------------
hb <- tab %>% select(gene, subtype, avg_log2FC) %>%
  tidyr::pivot_wider(names_from = subtype, values_from = avg_log2FC)
hb_mat <- as.matrix(hb[, -1]); rownames(hb_mat) <- hb$gene
hb_mat[is.na(hb_mat)] <- 0
hbl <- data.frame(gene = rownames(hb_mat), hb_mat, stringsAsFactors = FALSE) %>%
  tidyr::pivot_longer(-gene, names_to = "subtype", values_to = "value")
pB <- ggplot(hbl, aes(x = subtype, y = gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D62728",
                       name = "log2FC\n(OA - Normal)") +
  labs(title = "b  OA vs Normal differential expression (red = up in OA)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9.5, hjust = 0),
        axis.text.x = element_text(angle = 15, hjust = 1, size = 8),
        axis.text.y = element_text(size = 4.5),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))

## ---- Panel c: representative-gene violins -----------------------------------
neg <- tab$gene[tab$role == "negative regulator (druggable)"]
rep_genes <- unique(c(intersect(neg, c("TNFAIP3", "PLIN1", "CEBPA", "TLR4", "NFKB2")),
                      "RUNX2", "SP7", "PPARG"))
dat_l <- list()
for (g in rep_genes) {
  if (!(g %in% rownames(sub))) next
  ex <- as.numeric(GetAssayData(sub, layer = "data")[g, ])
  dat_l[[g]] <- data.frame(gene = g, expr = ex,
                           group = sub$group, subtype = sub$subtype)
}
dl <- bind_rows(dat_l)
pC <- ggplot(dl, aes(x = group, y = expr, fill = group)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, linewidth = 0.3) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = c(Normal = C_NORM, OA = C_OA), name = "Group") +
  labs(title = "c  Representative genes: OA vs Normal",
       y = "Log-normalized expression", x = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9.5, hjust = 0),
        legend.position = "top",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1))

## ---- assemble ----------------------------------------------------------------
fig <- pA / (pB | pC) + plot_layout(heights = c(1, 2.1), widths = c(1, 1.25))
ggsave(file.path(ROOT, "results/figures/Fig7_OA_validation.png"), fig,
       width = 7.2, height = 7.0, dpi = 300, bg = "white")
ggsave(file.path(ROOT, "results/figures/Fig7_OA_validation.pdf"), fig,
       width = 7.2, height = 7.0, bg = "white")
cat("Saved Fig7_OA_validation.png\n")
