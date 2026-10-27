## =============================================================================
## s07_oa_validation.R —— 路线 A 的 OA 外部验证
## =============================================================================
## 【为什么做】
##   路线 A 把 OA 从"主对象"降为"外部验证集"。核心假说：
##     骨修复/成骨命运受阻，部分是因为成骨「负调控因子」被异常强化（刹车踩死）。
##   虚拟敲除筛出的 ★★ 可药靶 = 这类负调控因子。若假说成立，应在独立 OA 数据里看到：
##     · 负调控因子（可药靶）在 OA 的基质/成骨细胞中 **上调**
##     · 成骨必需/标记基因在 OA 中 **下调或持平**
##   这构成「算法筛选 → 独立临床样本验证」的闭环。
##
## 【坑位47】S02 的 disease 字段：normal=9058，osteoarthritis=1148；
##   OA 的 Osteoblast 243 / Stromal_fibro 675 / CAR_cell 86，
##   成骨相关亚群在 OA 中显著富集（符合 OA 病理性异位成骨背景）。
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(tidyr)
})
ROOT <- "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa"
POP  <- c("CAR_cell", "Stromal_fibro", "Osteoblast")

## ---- 1. 取算法筛出的关键基因 ----------------------------------------------
sc <- read.csv(file.path(ROOT, "results/tables/s04_knk_scores_full.csv"),
               stringsAsFactors = FALSE)
neg <- sc %>% filter(direction == "抑制成骨" & p_perturb < 0.05) %>% pull(KO)   # 49 可药靶
pos <- sc %>% filter(direction == "促进成骨" & p_perturb < 0.05) %>% pull(KO)    # 成骨必需
ctrl <- c("RUNX2", "SP7", "PPARG", "CEBPB", "CXCL12", "LPL", "KITLG")
markers_bone <- c("IBSP", "BGLAP", "SPP1", "IFITM5", "OMD", "WIF1", "SFRP4", "CHAD")
genes <- unique(c(neg, pos, ctrl, markers_bone))
cat("待验证基因：", length(genes), "（负调控", length(neg), "/ 正调控", length(pos), ")\n")

## ---- 2. 加载 S02，标注 OA / Normal ---------------------------------------
obj <- readRDS(file.path(ROOT, "results/objects/s02_stromal.rds"))
obj$group <- ifelse(grepl("oa|arthr|osteo", tolower(as.character(obj$disease))),
                    "OA", "Normal")
sub <- subset(obj, cells = colnames(obj)[obj$subtype %in% POP &
                                         obj$group %in% c("OA", "Normal")])
cat("验证集细胞：", ncol(sub), "（OA", sum(sub$group == "OA"),
    "/ Normal", sum(sub$group == "Normal"), "）\n", sep = "")

## ---- 3. 逐群 OA vs Normal 差异（FindMarkers：OA - Normal）----------------
res <- list()
for (st in POP) {
  s <- subset(sub, subtype == st)
  if (ncol(s) < 20) { cat("  !", st, "细胞太少，跳过\n"); next }
  Idents(s) <- factor(s$group, levels = c("Normal", "OA"))
  mk <- tryCatch(FindMarkers(s, ident.1 = "OA", ident.2 = "Normal",
                             test.use = "wilcox", min.pct = 0.05, logfc.threshold = 0),
                 error = function(e) { cat("  !", st, "FindMarkers 失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(mk)) next
  mk$gene <- rownames(mk); mk$subtype <- st
  res[[st]] <- mk
}
res <- bind_rows(res)
cat("差异结果行数：", nrow(res), "\n", sep = "")

## ---- 4. 只挑我们关心的基因 + 方向标注 ------------------------------------
tab <- res %>% filter(gene %in% genes) %>%
  mutate(role = case_when(
    gene %in% neg      ~ "负调控(可药靶)",
    gene %in% pos      ~ "成骨必需",
    gene %in% ctrl     ~ "对照",
    TRUE               ~ "成骨标记"),
    sig = p_val_adj < 0.05) %>%
  select(gene, subtype, avg_log2FC, p_val, p_val_adj, role, sig) %>%
  arrange(subtype, role, desc(abs(avg_log2FC)))

## ---- 5. 出图 --------------------------------------------------------------
## Fig7A：三群 OA/Normal 细胞数对比
cnt <- as.data.frame(table(sub$subtype, sub$group)); colnames(cnt) <- c("subtype","group","n")
pA <- ggplot(cnt, aes(x = subtype, y = n, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(Normal = "#4C78A8", OA = "#E45756")) +
  labs(title = "A. 成骨相关亚群细胞数：OA vs Normal", y = "细胞数") +
  theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))

## Fig7B：关心基因的 avg_log2FC 热图（红=OA 高 / 蓝=OA 低）
hb <- tab %>% select(gene, subtype, avg_log2FC) %>%
  tidyr::pivot_wider(names_from = subtype, values_from = avg_log2FC)
hb_mat <- as.matrix(hb[, -1]); rownames(hb_mat) <- hb$gene
hb_mat[is.na(hb_mat)] <- 0
hbl <- data.frame(gene = rownames(hb_mat), hb_mat, stringsAsFactors = FALSE) %>%
  tidyr::pivot_longer(-gene, names_to = "subtype", values_to = "value")
pB <- ggplot(hbl, aes(x = subtype, y = gene, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#D62728",
                                     name = "log2FC\n(OA-Normal)") +
  labs(title = "B. OA vs Normal 差异表达（红=OA 中上调）", x = "", y = "") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(size = 7))

## Fig7C：代表性基因分组小提琴（负调控可药靶 + RUNX2/SP7）
rep_genes <- unique(c(intersect(neg, c("TNFAIP3","PLIN1","CEBPA","TLR4","NFKB2")),
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
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.3) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = c(Normal = "#4C78A8", OA = "#E45756")) +
  labs(title = "C. 代表性基因表达：OA vs Normal", y = "log-normalized expr") +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 20))

fig <- pA / (pB + pC)
ggsave(file.path(ROOT, "results/figures/Fig7_OA_validation.png"), fig,
       width = 11, height = 12, dpi = 150)

## ---- 6. 汇总输出 ----------------------------------------------------------
out <- file.path(ROOT, "results/tables/s07_oa_validation.csv")
write.csv(tab, out, row.names = FALSE)
cat("\n✓ 已保存：", out, "\n", sep = "")

## 方向自检：负调控因子在 OA 中是上调还是下调？
neg_tab <- tab %>% filter(role == "负调控(可药靶)")
cat("\n── 负调控(可药靶) 在 OA 中的方向 ──\n")
cat("  上调(log2FC>0) 显著:", sum(neg_tab$sig & neg_tab$avg_log2FC > 0, na.rm = TRUE),
    "/ 共", nrow(neg_tab), "\n")
cat("  下调(log2FC<0) 显著:", sum(neg_tab$sig & neg_tab$avg_log2FC < 0, na.rm = TRUE), "\n")
cat("\n── 成骨必需/标记 在 OA 中的方向 ──\n")
pos_tab <- tab %>% filter(role %in% c("成骨必需", "成骨标记"))
cat("  下调(成骨受阻) 显著:", sum(pos_tab$sig & pos_tab$avg_log2FC < 0, na.rm = TRUE),
    "/ 共", nrow(pos_tab), "\n")
