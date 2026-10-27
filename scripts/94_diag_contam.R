## =============================================================================
## 94_diag_contam.R —— 污染阈值诊断（不修改任何数据，只读+打印）
##
## 目的：S02 的 deplete_contamination 一刀切 0.25 阈值误杀了 71.6% 细胞，
##       其中 Erythroid 一项就杀了 13115 个。骨髓样本普遍含游离血红蛋白 mRNA，
##       基质细胞吸附后 HBB/HBA 呈中低水平阳性 —— 这是**背景**而非红细胞。
##       本脚本用数据确定"背景 vs 真污染"的分界。
##
## 判定逻辑：真污染细胞应当「谱系标记高」**且**「基质标记低」。
##           若一个细胞 HBB 高但 COL1A1 也高 → 是带背景的基质细胞，必须保留。
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr)
})

source("scripts/99_start.R")

## ---- 载入 ----
f <- file.path(PROJ$objects, "s01_h5ad.rds")
message("读 ", f)
seu <- readRDS(f)

STROMAL_TYPES <- c("stromal cell of bone marrow", "bone cell", "smooth muscle cell")
sub <- subset(seu, subset = cell_type %in% STROMAL_TYPES)
sub[["percent.mt"]] <- PercentageFeatureSet(sub, pattern = "^MT-")
sub <- subset(sub, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                            nCount_RNA >= 500 & percent.mt <= 15)
sub <- NormalizeData(sub, verbose = FALSE)
message("基质细胞：", ncol(sub))

d  <- GetAssayData(sub, slot = "data")
gn <- rownames(d)

mean_exp <- function(genes) {
  g2 <- intersect(genes, gn)
  if (!length(g2)) return(setNames(rep(0, ncol(d)), colnames(d)))
  if (length(g2) == 1L) setNames(as.numeric(d[g2, ]), colnames(d))
  else setNames(Matrix::colMeans(d[g2, , drop = FALSE]), colnames(d))
}

## ---- 各类污染标记的平均表达 ----
CM <- CONTAMINATION
sc <- as.matrix(vapply(CM, mean_exp, numeric(ncol(d))))
rownames(sc) <- colnames(d)

## ---- 基质身份标记（真正的基质/间质细胞必须表达这些）----
STROMA_ID <- c("COL1A1","COL1A2","COL3A1","DCN","LUM","CXCL12","PDGFRA",
               "THY1","ENG","NT5E","ACTA2","MYH11","SP7","RUNX2","ALPL")
STROMA_ID <- intersect(STROMA_ID, gn)
message("可用基质身份标记：", paste(STROMA_ID, collapse=", "))
stroma <- mean_exp(STROMA_ID)

## ---- 逐类诊断：被判污染的细胞，基质身份标记表达如何？----
message("\n================ 逐类阈值诊断 ================")
message("对每一类：看「该类标记高表达」的细胞里，基质身份标记是否也高。")
message("若基质标记也高 → 是带背景的基质细胞（误杀）；若基质标记极低 → 真污染。\n")

thr_grid <- c(0.25, 0.5, 1, 2, 3, 5, 8)

for (tp in names(CM)) {
  v <- sc[, tp]
  message("── ", tp, " ──")
  message(sprintf("   %8s %8s %10s %10s %10s",
                  "阈值", "命中数", "占比%", "基质标记均", "判定"))
  for (t in thr_grid) {
    hit <- v > t
    n <- sum(hit)
    if (n == 0) { message(sprintf("   %8.2f %8d %10.2f %10s %10s", t, 0, 0, "-", "-")); next }
    ms <- mean(stroma[hit])
    ## 对照：全体细胞的基质标记均值
    base <- mean(stroma)
    ratio <- ms / base
    verdict <- if (ratio > 0.6) "保留(基质样)" else if (ratio > 0.3) "可疑" else "剔除(真污染)"
    message(sprintf("   %8.2f %8d %10.2f %10.3f   %s (相对全体 %.2f×)",
                    t, n, 100 * n / length(v), ms, verdict, ratio))
  }
  message("")
}

## ---- 关键交叉表：Erythroid 高 vs 基质身份 ----
message("\n================ 重点：Erythroid（血红蛋白背景）================")
hb <- sc[, "Erythroid"]
st <- stroma
q  <- cut(hb, breaks = c(-1, 0.01, 0.25, 0.5, 1, 2, 4, 8, Inf),
          labels = c("0", "0-0.25", "0.25-0.5", "0.5-1", "1-2", "2-4", "4-8", ">8"))
tab <- data.frame(HB = q, stroma = st, nFeature = sub$nFeature_RNA,
                  nCount = sub$nCount_RNA)
agg <- tab %>% group_by(HB) %>%
  summarise(n = n(), stroma_mean = mean(stroma),
            nFeature_med = median(nFeature), nCount_med = median(nCount),
            .groups = "drop")
print(agg)

message("\n解读要点：")
message("  · 红细胞本尊：nFeature 很低(<500)、nCount 低、基质标记≈0")
message("  · 带血红蛋白背景的基质细胞：nFeature 正常(>1000)、基质标记明显>0")
message("  → 只有当 HB 极高且基质标记≈0 时才判红细胞")

## ---- 同样看 panImmune ----
message("\n================ 重点：panImmune（PTPRC/CD45 轴）================")
pi <- sc[, "panImmune"]
q2 <- cut(pi, breaks = c(-1, 0.01, 0.25, 0.5, 1, 2, 4, Inf),
          labels = c("0", "0-0.25", "0.25-0.5", "0.5-1", "1-2", "2-4", ">4"))
tab2 <- data.frame(PI = q2, stroma = st, nFeature = sub$nFeature_RNA)
print(tab2 %>% group_by(PI) %>%
        summarise(n = n(), stroma_mean = mean(stroma),
                  nFeature_med = median(nFeature), .groups = "drop"))

## ---- 建议阈值输出 ----
message("\n================ 建议 ===========================================")
message("把阈值改为「分级 + 双条件」：")
message("  免疫谱系(T/B/Plasma/Myeloid/Neutrophil)：均表达 > 0.25  且 基质标记 < 阈值")
message("  panImmune(PTPRC/CD45 轴)              ：> 0.5           且 基质标记 < 阈值")
message("  Erythroid(血红蛋白)                   ：> 3.0（骨髓基质有背景，须大幅放宽）")
message("  Platelet(PF4/PPBP)                    ：> 1.0")
message("\n具体分界点请以本脚本打印的交叉表为准（看 stroma_mean 掉到接近 0 的那一档）。")
