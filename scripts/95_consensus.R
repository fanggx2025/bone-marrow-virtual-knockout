## =============================================================================
## 95_consensus.R —— S06 三方共识：虚拟敲除 ∩ 轨迹 ∩ 网络药理学
## =============================================================================
## 【这是全文的落点（Fig6）】
##   立论是"用算法交叉筛选，找出哪些基因被虚拟敲除后能恢复成骨分化命运"。
##   单一来源的证据都不够硬，三方交叉才是卖点：
##
##     ① 虚拟敲除（scTenifoldKnk）：敲掉 g 后成骨模块显著被扰动，
##        且网络方向显示 g 抑制成骨 → g 是负调控因子 ★★
##     ② 轨迹（S03）：g 的表达随拟时序**下降**（rho < 0）
##        → 它在分化的路上被关掉，与"它是刹车"自洽
##     ③ 网络药理学：g 是刺五加/五味子某个成分的已知/预测靶点
##        → 药物有路可走
##
##   ① ∩ ② ∩ ③ = 核心靶点（药→靶→命运 换向的完整闭环）
##   ① ∩ ②     = 算法共识靶点（无药，但可作机制论述）
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix)
})
source("scripts/99_start.R")

POP <- c("CAR_cell", "Stromal_fibro", "Osteoblast")

## ---- 1. 虚拟敲除结果 -------------------------------------------------------
sc <- read.csv(file.path(PROJ$tables, "s04_knk_scores_full.csv"), stringsAsFactors = FALSE)
ko <- sc %>% filter(direction == "抑制成骨" & p_perturb < 0.05) %>%
  arrange(p_perturb)
cat("虚拟敲除 ★★ 可药靶：", nrow(ko), " 个\n", sep = "")

## ---- 2. 轨迹方向（对每个 KO 命中基因重算 rho，而不是只看候选池那 40 个）----
obj <- readRDS(file.path(PROJ$objects, "s02_stromal.rds"))
obj <- subset(obj, cells = colnames(obj)[obj$subtype %in% POP])

ptab <- read.csv(file.path(PROJ$tables, "s03_pseudotime_cells.csv"),
                 stringsAsFactors = FALSE)
keep <- intersect(ptab$cell, colnames(obj))
pt   <- ptab$pseudotime[match(keep, ptab$cell)]

dat <- GetAssayData(obj, layer = "data")[, keep, drop = FALSE]
gs  <- intersect(ko$KO, rownames(dat))
cat("可算轨迹相关的 KO 命中基因：", length(gs), " / ", nrow(ko), "\n", sep = "")

rho <- sapply(gs, function(g) {
  x <- as.numeric(dat[g, ])
  if (sd(x) == 0) return(NA_real_)
  suppressWarnings(cor(x, pt, method = "spearman"))
})
cnt <- Matrix::rowSums(GetAssayData(obj, layer = "counts")[gs, ] > 0) / ncol(dat)

traj <- data.frame(gene = gs, rho_pseudotime = as.numeric(rho),
                   pct_detected = as.numeric(cnt[gs]),
                   stringsAsFactors = FALSE) %>%
  mutate(轨迹下降 = !is.na(rho_pseudotime) & rho_pseudotime < 0)

cat("其中随拟时序下降（rho<0）：", sum(traj$轨迹下降, na.rm = TRUE), " 个\n", sep = "")

## ---- 3. 网络药理学靶点 -----------------------------------------------------
f_tg <- file.path(PROJ$tables, "s05_compound_targets.csv")
tg <- if (file.exists(f_tg)) read.csv(f_tg, stringsAsFactors = FALSE) else NULL
if (!is.null(tg) && nrow(tg)) {
  cat("成分-靶点对：", nrow(tg), " | 去重靶点：", length(unique(tg$gene_symbol)), "\n", sep = "")
} else {
  cat("⚠ 没有成分靶点表（s05_compound_targets.csv），三方退化为两方\n")
}

## ---- 4. 拼装 ---------------------------------------------------------------
## 【坑位45】ko 表从 s04_knk_scores_full.csv 读进来，基因列叫 KO；
##   traj 表叫 gene。left_join 必须显式 by = c("KO" = "gene")，
##   否则 dplyr 会报 "Join columns must be present in data"。
tab <- ko %>% left_join(traj, by = c("KO" = "gene")) %>% rename(gene = KO)

if (!is.null(tg) && nrow(tg)) {
  hit_tg <- tg %>% filter(gene_symbol %in% tab$gene)
  agg <- hit_tg %>% group_by(gene_symbol) %>%
    summarise(herbs = paste(sort(unique(herb)), collapse = "+"),
              compounds = paste(sort(unique(compound)), collapse = "; "),
              n_compound = n_distinct(compound),
              best_nM = min(best_nM, na.rm = TRUE),
              .groups = "drop")
  tab <- tab %>% left_join(agg, by = c("gene" = "gene_symbol")) %>%
    mutate(有药 = !is.na(herbs))
} else {
  tab <- tab %>% mutate(herbs = NA_character_, compounds = NA_character_,
                        n_compound = NA_integer_, best_nM = NA_real_,
                        有药 = FALSE)
}

tab <- tab %>% mutate(
  n_source = 1 + as.integer(轨迹下降) + as.integer(有药),
  tier = case_when(
    有药 & 轨迹下降 ~ "Tier1 三方共识（核心靶点）",
    有药            ~ "Tier2 药+敲除",
    轨迹下降        ~ "Tier3 敲除+轨迹（机制）",
    TRUE            ~ "Tier4 仅敲除"
  )) %>%
  arrange(tier, p_perturb)

cat("\n════════ 分层结果 ════════\n")
print(table(tab$tier))

cat("\n──────── Tier1 三方共识（核心靶点）────────\n")
t1 <- tab %>% filter(tier == "Tier1 三方共识（核心靶点）")
if (nrow(t1)) {
  print(t1 %>% select(gene, n_compound, best_nM, rho_pseudotime,
                      pct_rank_end, p_perturb, med_w_end, herbs, compounds),
        row.names = FALSE, width = 200)
} else cat("（无）\n")

cat("\n──────── Tier3 敲除+轨迹（算法共识，无药）前 20 ────────\n")
t3 <- tab %>% filter(tier == "Tier3 敲除+轨迹（机制）")
if (nrow(t3)) print(head(t3 %>% select(gene, rho_pseudotime, pct_rank_end,
                                       p_perturb, med_w_end), 20), row.names = FALSE)

## ---- 5. 输出 ---------------------------------------------------------------
f <- file.path(PROJ$tables, "s06_consensus_targets.csv")
write.csv(tab, f, row.names = FALSE)
cat("\n✓ 已保存：", f, "\n", sep = "")

saveRDS(list(tab = tab, traj = traj,
             tg = tg), file.path(PROJ$objects, "s06_consensus.rds"))
