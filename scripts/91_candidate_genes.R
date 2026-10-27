## =============================================================================
## 91_candidate_genes.R —— 构建「虚拟敲除候选基因池」
## =============================================================================
## 【本文立论】不是描述骨不愈中成骨谱系如何停滞（已撞题），而是：
##   **用算法交叉筛选，找出哪些基因被虚拟敲除后能恢复成骨分化命运。**
##
## 候选池的质量决定 Fig5 的成败。三个来源，缺一不可：
##
##   A. 文献先验 —— 已知的成骨负调控因子（有据可循，审稿人认）
##   B. 数据驱动 —— 本数据里「随成骨推进而下降」的基因
##                  （轨迹提供的候选，是本文"交叉筛选"的一半）
##   C. 对照     —— RUNX2/SP7 应为 score<0（成骨必需）
##                  PPARG     应为 score>0（成脂主控，抑制成骨）
##                  若对照组方向跑反，说明整条流程有问题，必须回头排查
##
## 【判读方向 —— 最容易搞反的一步】
##   敲除 X 后成骨标记**上调** → X 是成骨的**负调控因子** → 可药靶 ★
##   敲除 X 后成骨标记**下调** → X 是成骨所**必需**     → 不可药靶
##
## 【为什么要排除 end_markers】
##   成骨终点标记（IBSP/BGLAP/SPP1/...）既是"敲除效果的读出指标"，
##   若又拿它们当敲除对象，等于用同一批基因既当自变量又当因变量 → 循环论证。
##   起点标记（CXCL12/LPL/KITLG）不受此限：敲掉它们是为了把细胞"推出"静息态，
##   这正是我们想验证的假设。故 root_markers 保留、end_markers 排除。
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix)
})
source("scripts/99_start.R")

POP <- c("CAR_cell", "Stromal_fibro", "Osteoblast")

## =============================================================================
## A. 文献先验：成骨负调控因子（人基因符号，按通路分组）
## =============================================================================
PRIOR <- list(
  `WNT抑制`   = c("DKK1","DKK2","DKK3","SFRP1","SFRP2","SOST","GSK3B","APC",
                  "AXIN1","AXIN2","CTNNBIP1","TCF7L1","WIF1","NOTUM","ZNRF3","RNF43"),
  `Notch`     = c("HES1","HES4","HES5","HEY1","HEY2","NOTCH1","NOTCH2","JAG1",
                  "DTX1","LFNG","MAML1","RBPJ"),
  `成脂分化`   = c("PPARG","CEBPA","CEBPB","CEBPD","FABP4","ADIPOQ","PLIN1","PLIN2",
                  "LPL","ZNF423","KLF15","SREBF1","NR1H3","RETN","CIDEC"),
  `炎症NFkB`   = c("NFKB1","NFKB2","RELA","RELB","TNF","IL1B","IL6","PTGS2",
                  "TNFAIP3","IKBKB","CHUK","NFKBIA","NFKBIZ","TLR4","MYD88"),
  `软骨决定`   = c("SOX9","SOX5","SOX6","COL2A1","ACAN","COL9A1","COL11A1","MATN3"),
  `抑制型TF`   = c("TWIST1","TWIST2","SNAI1","SNAI2","ZEB1","ZEB2","ID1","ID2",
                  "ID3","ID4","JUNB","FOSL1","FOSL2","ATF3","EGR1","MAFB","IRF8",
                  "HAND2","EN1","PAX1","FOXC1","FOXC2","PRDM16","GAS1","BCL6","MEF2C"),
  `BMPSMAD抑制`= c("NOG","GREM1","GREM2","CHRD","BMPER","FST","INHBA","CER1",
                  "SMAD6","SMAD7","SKI","SKIL","TGIF1","LTBP2","TGFB1","TGFB3","TGFBR3"),
  `衰老SASP`   = c("CDKN1A","CDKN2A","TP53","MDM2","SERPINE1","IGFBP7","IGFBP3",
                  "LMNB1","GDF15","CCND1","MKI67"),
  `基质降解OA` = c("MMP1","MMP2","MMP3","MMP13","MMP14","ADAMTS4","ADAMTS5",
                  "ADAMTS1","TIMP1","TIMP2","SERPINA1"),
  `表观`       = c("EZH2","HDAC4","HDAC5","HDAC9","DNMT3A","KDM1A","EHMT2","SUV39H1"),
  `其他`       = c("VHL","HIF1A","EPAS1","MYC","FGF2","PTHLH","TNFRSF11B","LRP5",
                  "SOST","MEPE","PHEX","DMP1","ANKH","ENPP1")
)

## 对照：方向已知的基因，用来验证流程没跑反
CONTROLS <- c(
  "RUNX2",   # 成骨主控 TF → 敲掉应使成骨标记**下调** → score < 0（阳性对照）
  "SP7",     # Osterix，同上 → score < 0
  "PPARG",   # 成脂主控，与成骨互斥 → 敲掉应使成骨标记**上调** → score > 0
  "CEBPB",   # 成脂，同上预期 > 0
  "CXCL12",  # 静息态 BMSC 标记 → 敲掉把细胞推出静息 → 预期 > 0
  "LPL",     # 起点标记 + 成脂 → 预期 > 0
  "KITLG"    # 起点标记 → 预期 > 0
)

## =============================================================================
## B. 数据驱动：随拟时序显著下降的基因
## =============================================================================
## 【坑位36】s03_trajectory.rds 是 monocle 的 CellDataSet，不是 Seurat 对象，
##   用 obj$pseudotime / GetAssayData 都会报错（第一版因此静默返回 0 个基因）。
##   改用 S03 落盘的 s03_pseudotime_cells.csv（cell → pseudotime）与 S02 表达矩阵 join，
##   绕开 monocle 的 API，也更容易复核。
data_driven <- function(obj_seu, n_top = 40L, min_pct = 0.05, max_rho = -0.12) {
  cat("\n→ 数据驱动：与拟时序负相关的基因\n")

  f <- file.path(PROJ$tables, "s03_pseudotime_cells.csv")
  if (!file.exists(f)) { cat("  ⚠ 找不到 ", f, "\n", sep = ""); return(NULL) }

  ptab <- read.csv(f, stringsAsFactors = FALSE)
  ptab <- ptab[!is.na(ptab$pseudotime), ]
  keep <- intersect(ptab$cell, colnames(obj_seu))
  cat("  轨迹细胞：", nrow(ptab), " → 与 S02 交集 ", length(keep), "\n")
  if (length(keep) < 100) return(NULL)

  pt <- ptab$pseudotime[match(keep, ptab$cell)]
  d  <- GetAssayData(obj_seu, layer = "data")[, keep, drop = FALSE]

  n   <- ncol(d)
  pct <- Matrix::rowSums(d > 0) / n
  use <- names(pct)[pct >= min_pct & pct <= 0.98]
  cat("  可测基因（检出率 5%~98%）：", length(use), "\n")

  ## 逐基因 Spearman（秩相关，避免表达量分布拖累）
  d2  <- as.matrix(d[use, ])
  rho <- apply(d2, 1, function(x) {
    if (stats::sd(x) == 0) return(NA_real_)
    suppressWarnings(stats::cor(x, pt, method = "spearman"))
  })

  res <- data.frame(gene = use, rho = as.numeric(rho),
                    pct = as.numeric(pct[use]), stringsAsFactors = FALSE) %>%
    filter(!is.na(rho)) %>% arrange(rho)

  cat("  rho 最小 12 个：",
      paste(sprintf("%s(%.3f)", head(res$gene, 12), head(res$rho, 12)), collapse = " "), "\n")

  hit <- res %>% filter(rho <= max_rho) %>% head(n_top)
  cat("  入选（rho ≤ ", max_rho, "）：", nrow(hit), " 个\n")
  hit
}

## =============================================================================
## 主流程
## =============================================================================
build_candidates <- function(n_top = 40L, save = TRUE) {
  message("\n########## 构建虚拟敲除候选基因池 ##########")

  ## 读 S02 基质对象（用于检出率过滤）
  obj <- readRDS(file.path(PROJ$objects, "s02_stromal.rds"))
  keep <- obj$subtype %in% POP
  objp <- subset(obj, cells = colnames(obj)[keep])
  cat("用于评估检出率的细胞：", ncol(objp), "\n", sep = "")

  ## 检出率（在原始计数上算，衡量"有多少细胞能测到"）
  cnt <- GetAssayData(objp, slot = "counts")
  n   <- ncol(cnt)
  pct <- Matrix::rowSums(cnt > 0) / n

  ## 数据驱动候选（用 S03 落盘的拟时序表 + S02 表达矩阵，绕开 monocle API）
  dd <- character(0); dd_df <- NULL
  dd_df <- tryCatch(data_driven(objp, n_top = n_top), error = function(e) {
    cat("  ⚠ 数据驱动部分失败：", conditionMessage(e), "\n"); NULL })
  if (!is.null(dd_df)) dd <- dd_df$gene

  prior <- unique(unlist(PRIOR))
  ctrl  <- CONTROLS

  ## ---- 组装 ----
  pool <- unique(c(ctrl, prior, dd))

  ## 排除成骨终点标记（避免循环论证）
  excl <- TRAJECTORY$end_markers
  pool <- setdiff(pool, excl)

  ## 只保留数据中存在的基因
  in_data <- intersect(pool, rownames(objp))
  cat("\n候选池汇总：\n")
  cat("  文献先验：", length(prior), " → 存在于数据 ", length(intersect(prior, rownames(objp))), "\n")
  cat("  对照组  ：", length(ctrl),  " → 存在于数据 ", length(intersect(ctrl, rownames(objp))), "\n")
  cat("  数据驱动：", length(dd),    " → 存在于数据 ", length(intersect(dd, rownames(objp))), "\n")
  cat("  排除终点标记：", paste(intersect(pool, excl), collapse = ", "), "\n")
  cat("  最终候选：", length(in_data), " 个\n")

  tab <- data.frame(
    gene   = in_data,
    source = sapply(in_data, function(g) {
      s <- c(if (g %in% ctrl) "对照",
             if (g %in% prior) "文献先验",
             if (g %in% dd)    "数据驱动")
      paste(s, collapse = "+")
    }),
    pct_detected = as.numeric(pct[in_data]),
    mean_count   = as.numeric(Matrix::rowMeans(cnt[in_data, ])),
    stringsAsFactors = FALSE
  )
  if (!is.null(dd_df)) {
    tab$rho_pseudotime <- as.numeric(dd_df$rho[match(tab$gene, dd_df$gene)])
  }

  ## 关键预警：检出率低于 5% 的基因会被 scQC 滤掉 / 建网不可靠
  tab$qc_flag <- ifelse(tab$pct_detected < 0.05, "⚠低检出(<5%)",
                 ifelse(tab$pct_detected < 0.10, "偏低(5-10%)", "OK"))

  cat("\n检出率分布：\n")
  print(table(tab$qc_flag))

  ## 强制纳入的基因（成骨/轨迹标记必须留在网络里当"读出指标"）
  force <- unique(c(TRAJECTORY$end_markers, TRAJECTORY$root_markers,
                    TRAJECTORY$mid_markers, in_data))

  if (save) {
    dir.create(PROJ$tables, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(PROJ$tables, "s04_candidate_genes.csv")
    write.csv(tab, f, row.names = FALSE)
    cat("\n✓ 已保存：", f, "\n")
    saveRDS(list(tab = tab, force = force,
                 prior = PRIOR, controls = CONTROLS),
            file.path(PROJ$objects, "s04_candidates.rds"))
    cat("✓ 已保存候选对象\n")
  }
  invisible(list(tab = tab, force = force))
}

res <- build_candidates()
cat("\n──────── 候选池（前 40，按来源/检出率） ────────\n")
print(head(res$tab %>% arrange(-pct_detected), 40), row.names = FALSE)
