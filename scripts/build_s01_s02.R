## =============================================================================
## build_s01_s02.R —— 从 CELLxGene atlas 重建 S01（全量 Seurat）与 S02（基质细胞）
##
## 【为什么有这个脚本】
##   原 S01->S02 的构建（QC、基质门控、Harmony 整合、聚类、5 亚群注释）当年是在
##   RStudio 里交互式完成的，未落成脚本。本脚本按稿件 Methods 文档的流程【复现】：
##     1) 读入 atlas h5ad（复用 04_load_h5ad.R 的 load_h5ad，已处理 4 个 HDF5 坑）
##     2) QC：线粒体比例 < 20%
##     3) 基质门控：config.R 的 MARKERS_HUMAN，基质谱系阳性模块分 > 非基质/免疫/
##         红细胞阴性模块分（绝对表达硬过滤 + 身份双条件）
##     4) Seurat v5 标准流程 + Harmony 按 gse 批次校正（参数 reduction.use）
##     5) 分辨率 0.5 聚类 + UMAP
##     6) 按 5 个骨髓基质亚群标记（CAR/Pericyte_SMC/Stromal_fibro/Osteoblast/
##        Chondro_prog）的模块分注释 subtype
##     7) 输出 s01_h5ad.rds / s02_stromal.rds
##
## 【与原始交互式构建的差异声明】
##   注释阈值等细节无法从文档完全恢复，边缘细胞的归类可能与原始构建有出入；
##   稿件全部表格所用的确切 S02 已通过 Zenodo 衍生数据集存档（DOI 见 README）。
##   本脚本用于复现主线流程，不作为逐细胞一致的证明。
##
## 【用法】在仓库根目录：
##   Rscript scripts/build_s01_s02.R
## =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(dplyr)
})

## ---- 0. 路径与参数（复刻 config.R 的用法） ----
e <- new.env()
sys.source("config/config.R", envir = e)
PROJ    <- get("PROJ", envir = e)
MARKERS <- get("MARKERS", envir = e)   # = MARKERS_HUMAN
set.seed(get("SEED", envir = e))

H5AD <- file.path(PROJ$data_raw, "bone_marrow_ref.h5ad")
stopifnot("atlas h5ad not found - see README data provenance" = file.exists(H5AD))

## ---- 1. 读入（04_load_h5ad.R 提供 load_h5ad，处理分类列/基因名等坑） ----
source("scripts/04_load_h5ad.R")
s01 <- load_h5ad(H5AD, "BoneMarrow")
cat("S01 loaded:", ncol(s01), "cells x", nrow(s01), "genes\n")

## ---- 2. QC：线粒体 < 20% ----
s01[["percent.mt"]] <- PercentageFeatureSet(s01, pattern = "^MT-")
s01 <- subset(s01, percent.mt < 20)
cat("after mito QC:", ncol(s01), "cells\n")

## ---- 3. 基质门控：基质阳性 vs 非基质（免疫/红细胞/内皮）双条件 ----
stromal_pos <- unique(unlist(MARKERS[c("MSC_generic","Stromal_fibro","CAR_cell",
                                       "Adipo_prog","Osteo_prog","Osteoblast",
                                       "Chondro_prog","Pericyte_SMC")]))
nonstromal   <- unique(unlist(MARKERS[c("Hematopoietic","Monocyte_Mac","T_cell","B_cell",
                                        "Plasma_cell","NK_cell","Neutrophil",
                                        "Erythroid","Endothelial")]))
stromal_pos <- intersect(stromal_pos, rownames(s01))
nonstromal  <- intersect(nonstromal,  rownames(s01))
s01 <- AddModuleScore(s01, features = list(stromal_pos), name = "stroma_pos")
s01 <- AddModuleScore(s01, features = list(nonstromal),  name = "nonstroma")
s02 <- subset(s01, stroma_pos1 > nonstroma1)
cat("after stromal gate:", ncol(s02), "cells\n")

## ---- 4. 整合：Normalize -> HVG -> Scale -> PCA -> Harmony(按 gse) ----
s02 <- NormalizeData(s02)
s02 <- FindVariableFeatures(s02, selection.method = "vst", nfeatures = 2000)
s02 <- ScaleData(s02, vars.to.regress = "percent.mt")
s02 <- RunPCA(s02, features = VariableFeatures(s02))
## 注意：harmony::RunHarmony 的参数名是 reduction.use（不是 reduction）
s02 <- RunHarmony(s02, group.by.vars = "gse", reduction.use = "pca", dims.use = 1:30)

## ---- 5. 聚类 + UMAP ----
s02 <- RunUMAP(s02, reduction = "harmony", dims = 1:30)
s02 <- FindNeighbors(s02, reduction = "harmony", dims = 1:30)
s02 <- FindClusters(s02, resolution = 0.5)

## ---- 6. 五亚群注释（模块分 argmax，低于门槛归 Stromal_fibro 兜底） ----
SUB <- list(
  CAR_cell      = MARKERS$CAR_cell,
  Pericyte_SMC  = MARKERS$Pericyte_SMC,
  Stromal_fibro = MARKERS$Stromal_fibro,
  Osteoblast    = MARKERS$Osteoblast,
  Chondro_prog  = MARKERS$Chondro_prog
)
for (nm in names(SUB)) {
  genes <- intersect(SUB[[nm]], rownames(s02))
  s02 <- AddModuleScore(s02, features = list(genes), name = paste0("sub_", nm))
}
score_cols <- paste0("sub_", names(SUB), "1")
sc <- as.data.frame(s02@meta.data[, score_cols, drop = FALSE])
colnames(sc) <- names(SUB)
s02$subtype <- ifelse(apply(sc, 1, max) < 0.05, "Stromal_fibro",
                      names(SUB)[apply(sc, 1, which.max)])
cat("subtype counts:\n"); print(table(s02$subtype))

## ---- 7. 输出 ----
dir.create(PROJ$objects, showWarnings = FALSE, recursive = TRUE)
saveRDS(s01, file.path(PROJ$objects, "s01_h5ad.rds"))
saveRDS(s02, file.path(PROJ$objects, "s02_stromal.rds"))
cat("saved: s01_h5ad.rds / s02_stromal.rds ->", PROJ$objects, "\n")
