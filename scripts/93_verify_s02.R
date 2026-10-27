## =============================================================================
## 93_verify_s02.R —— S02 结果验证：不看注释标签，只看 marker 与关键基因表达
##
## 【为什么必须这一步】
## 第一轮 S02 的注释里 C3 被打成 "Plasma_cell"，看起来"注释成功"，
## 但一看 marker（IGHGP/MZB1）才发现 29% 是污染。
## → 结论：**注释标签不可信，marker 才可信。**
##
## 本脚本回答三个决定路线 A 生死的问题：
##   Q1. 各亚群到底是什么？（看 top marker）
##   Q2. 成骨谱系存不存在？（BGLAP/SPP1/RUNX2/SP7/ALPL 的表达梯度）
##   Q3. 成脂谱系存不存在？（ADIPOQ/PPARG/FABP4/CEBPA/LEP）
##   Q4. 供体与疾病构成如何？（路线 A 要求 normal 供体足够多）
## =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
source("scripts/99_start.R")

f <- file.path(PROJ$objects, "s02_stromal.rds")
message("读 ", f)
s <- readRDS(f)
message("细胞数：", ncol(s), "  基因数：", nrow(s))

## ---------------- Q1. 各亚群 top marker ----------------
message("\n################ Q1. 各亚群 TOP marker ################")
mk <- read.csv(file.path(PROJ$tables, "s02_subcluster_markers.csv"),
               stringsAsFactors = FALSE)
## 【坑位19】seurat_clusters 是 factor，as.integer(factor) 拿到的是**水平序号**
##   而非数值本身 → 会整体偏移、还会跳过 cluster 0。
##   必须用 as.numeric(as.character(...))。
cls <- sort(unique(as.numeric(as.character(s$seurat_clusters))))
cls <- cls[!is.na(cls)]
for (cl in cls) {
  s2 <- mk[mk$cluster == cl, ]
  s2 <- s2[order(-s2$avg_log2FC), ]
  lab <- if ("subtype" %in% colnames(s@meta.data)) {
    u <- unique(s$subtype[s$seurat_clusters == cl]); paste(u, collapse = "/")
  } else "-"
  n <- sum(s$seurat_clusters == cl)
  cat(sprintf("C%-3d n=%-5d [%s]\n     %s\n", cl, n, lab,
              paste(head(s2$gene, 14), collapse = ", ")))
}

## ---------------- Q2/Q3. 关键谱系基因的表达梯度 ----------------
message("\n################ Q2/Q3. 成骨 / 成脂关键基因表达 ################")
PANELS <- list(
  "成骨-早期祖" = c("PRRX1","TWIST1","PDGFRA","CXCL12","NT5E","THY1","ENG"),
  "成骨-承诺"   = c("RUNX2","SP7","ALPL","COL1A1","IBSP","SPP1"),
  "成骨-成熟"   = c("BGLAP","IFITM5","DMP1","PHEX","SOST","MEPE"),
  "成脂"       = c("ADIPOQ","PPARG","CEBPA","FABP4","PLIN1","LEP","LPL"),
  "软骨"       = c("SOX9","ACAN","COL2A1","COL11A1","MATN3"),
  "周细胞/SMC" = c("ACTA2","TAGLN","MYH11","PDGFRB","CNN1","RGS5"),
  "内皮"       = c("PECAM1","VWF","CDH5","KDR"),
  "残留污染"   = c("PTPRC","CD3E","MS4A1","MZB1","LYZ","C1QA","MPO","PF4")
)

d  <- GetAssayData(s, slot = "data")
gn <- rownames(d)
cl <- as.character(s$seurat_clusters)

for (pn in names(PANELS)) {
  g <- intersect(PANELS[[pn]], gn)
  if (!length(g)) { message("\n[", pn, "] 无可用基因"); next }
  ## 每个亚群内该组基因的平均表达（先按细胞求均值，再按群求均值）
  cellmean <- if (length(g) == 1L) as.numeric(d[g, ]) else Matrix::colMeans(d[g, , drop = FALSE])
  bycl <- tapply(cellmean, cl, mean)
  ## 阳性细胞比例
  pct <- tapply(cellmean, cl, function(x) mean(x > 0))
  message("\n[", pn, "]  （mean 表达 | 阳性细胞%）")
  for (i in seq_along(bycl))
    message(sprintf("   C%-3s %6.3f  %5.1f%%", names(bycl)[i], bycl[i], 100 * pct[i]))
  ## 逐基因明细（只在小面板时打）
  if (length(g) <= 7) {
    m2 <- as.matrix(t(apply(as.matrix(d[g, , drop = FALSE]), 1,
                            function(x) tapply(x, cl, mean))))
    message("   逐基因：")
    for (gi in rownames(m2))
      message(sprintf("     %-8s %s", gi,
                      paste(sprintf("%s:%.2f", colnames(m2), m2[gi, ]), collapse = " ")))
  }
}

## ---------------- Q4. 供体 / 疾病构成 ----------------
message("\n################ Q4. 供体与疾病构成 ################")
if ("disease" %in% colnames(s@meta.data)) {
  message("\n各亚群的 disease 构成（行=亚群，列=disease）：")
  print(table(cl, s$disease))
  message("\n总体 disease："); print(table(s$disease))
}
if ("donor_id" %in% colnames(s@meta.data)) {
  nd <- length(unique(s$donor_id))
  message("\n供体总数：", nd)
  if ("disease" %in% colnames(s@meta.data)) {
    dn <- s@meta.data %>% distinct(donor_id, disease)
    message("按疾病分的供体数："); print(table(dn$disease))
    message("\n各供体细胞数：")
    print(sort(table(s$donor_id), decreasing = TRUE))
  }
}

## ---------------- 轨迹可行性判据 ----------------
message("\n################ 轨迹可行性判据 ################")
osteo <- intersect(c("RUNX2","SP7","ALPL","COL1A1","BGLAP","IBSP","SPP1","IFITM5"), gn)
adipo <- intersect(c("ADIPOQ","PPARG","CEBPA","FABP4","PLIN1"), gn)
root  <- intersect(c("PRRX1","CXCL12","PDGFRA","NT5E","THY1"), gn)

sc_osteo <- if (length(osteo) == 1L) as.numeric(d[osteo, ]) else Matrix::colMeans(d[osteo, , drop=FALSE])
sc_adipo <- if (length(adipo) == 1L) as.numeric(d[adipo, ]) else Matrix::colMeans(d[adipo, , drop=FALSE])
sc_root  <- if (length(root)  == 1L) as.numeric(d[root, ])  else Matrix::colMeans(d[root,  , drop=FALSE])

res <- data.frame(
  cluster = names(tapply(sc_osteo, cl, mean)),
  root    = as.numeric(tapply(sc_root,  cl, mean)),
  osteo   = as.numeric(tapply(sc_osteo, cl, mean)),
  adipo   = as.numeric(tapply(sc_adipo, cl, mean))
)
res$osteo_vs_adipo <- res$osteo - res$adipo
res <- res[order(-res$osteo_vs_adipo), ]
message("成骨轴 − 成脂轴 排序（正=偏成骨，负=偏成脂）：")
print(res, digits = 3)

span_osteo <- diff(range(res$osteo))
span_adipo <- diff(range(res$adipo))
message("\n成骨评分极差：", round(span_osteo, 3), "   成脂评分极差：", round(span_adipo, 3))
message("\n判据：")
message("  成骨极差 > 0.3 → 成骨分化轴可建（Monocle DDRTree 能找到方向）")
if (span_osteo > 0.3) message("  ✓ 成骨轴 OK") else message("  ✗ 成骨梯度太弱，轨迹建不起来")
message("  成脂极差 > 0.2 → 成脂端可区分")
if (span_adipo > 0.2) message("  ✓ 成脂轴 OK") else message("  ⚠ 成脂梯度弱：骨髓抽吸物中成熟脂肪细胞本就稀少")
