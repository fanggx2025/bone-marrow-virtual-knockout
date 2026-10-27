## =============================================================================
## 90_calibrate_knk.R —— S04 虚拟敲除「耗时标定」（第二版）
## =============================================================================
## 【坑位31（本版修复，非常关键）】
##   第一版直接把 scTenifoldKnk 的 scQC 套在**子集矩阵**上，结果：
##     Library size filter (min=1000): 548/7880 cells retained   ← 只剩 7%！
##     Outlier removal: 502 cells retained
##   根因：scQC 的 minLibSize=1000 是照着**全转录组**定的阈值。
##   我们为了控制网络规模只取了 600~1500 个高变基因，
##   于是每个细胞的"文库大小"只剩这几百个基因的 UMI 之和，
##   大部分细胞自然达不到 1000 → 被当成低质量细胞砍掉。
##   同理 minPCT=0.05 的基因过滤也会砍掉一大片（445/600 保留）。
##
##   ⇒ 修法：**跳过 scQC**。细胞 QC 已在 S02 做过，这里只做：
##        ① 基因表达广度过滤（≥ min_cells 个细胞）
##        ② CPM 归一化（与官方 cpmNormalization 完全一致：t(t(X)/colSums(X))*1e6）
##
## 【坑位32】cpmNormalization 在 scTenifoldKnk 里没导出，
##   实际定义在 scTenifoldNet 命名空间。用 scTenifoldNet::cpmNormalization。
##
## 【核心策略：建网一次，复用 WT 张量】
##   读官方源码可知，7 步里只有 Step3(makeNetworks) / Step4(tensorDecomposition)
##   昂贵，且**与被敲的基因无关**；Step5-7 只是：
##       KO <- WT; KO[g, ] <- 0; MA <- manifoldAlignment(WT, KO); DR <- dRegulation(MA)
##   官方的 transcriptomeWide=TRUE 就是这么做的，但它只保留 distance、
##   **丢掉了 FC / Z / p.adj**，而我们判读"成骨标记上调还是下调"必须有 FC。
##   ⇒ 于是自己复现同一条流水线：建网一次 → 存 WT → 逐基因复用。
##   算法步骤与随机种子(set.seed(1))与官方逐基因调用完全一致。
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(scTenifoldKnk); library(scTenifoldNet); library(Matrix)
})
source("scripts/99_start.R")

POP <- c("CAR_cell", "Stromal_fibro", "Osteoblast")

obj <- readRDS(file.path(PROJ$objects, "s02_stromal.rds"))
obj <- subset(obj, cells = colnames(obj)[obj$subtype %in% POP])
cat("细胞群：", paste(POP, collapse = " / "), " → ", ncol(obj), " 细胞\n", sep = "")

build_mat <- function(obj, hvg_n = 1500L, force = NULL, min_cells = 20L) {
  cnt <- GetAssayData(obj, layer = "counts")
  if (!inherits(cnt, "dgCMatrix")) cnt <- as(as.matrix(cnt), "dgCMatrix")
  n_expr <- Matrix::rowSums(cnt > 0)
  cnt <- cnt[n_expr >= min_cells, , drop = FALSE]
  hv <- VariableFeatures(obj)
  sel <- intersect(hv, rownames(cnt))
  if (!is.null(force)) sel <- union(sel, intersect(force, rownames(cnt)))
  cnt <- cnt[sel, , drop = FALSE]
  m <- as.matrix(cnt); storage.mode(m) <- "double"; m
}

FORCE <- unique(c("RUNX2","SP7","PPARG","CXCL12",
                  TRAJECTORY$end_markers, TRAJECTORY$root_markers,
                  TRAJECTORY$mid_markers))

tmr <- function(label, expr) {
  t0 <- Sys.time(); v <- force(expr)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  [%7.1f s] %s\n", el, label))
  invisible(list(value = v, secs = el))
}

NCELLS <- 500L
for (P in c(600L, 1200L, 2000L)) {
  cat("\n══════════════ 基因数 p =", P, "══════════════\n")
  m <- build_mat(obj, hvg_n = P, force = FORCE)
  if (nrow(m) > P) m <- m[1:P, ]
  cat("  矩阵：", nrow(m), " 基因 × ", ncol(m), " 细胞  (",
      round(as.numeric(object.size(m))/1e6, 1), " MB)\n", sep = "")
  cat("  文库中位数（子集基因）：", median(Matrix::colSums(m)), "\n", sep = "")

  r2 <- tmr("cpmNormalization", scTenifoldNet::cpmNormalization(m))
  mc <- r2$value

  r3 <- tmr("makeNetworks(nNet=1)",
            scTenifoldNet::makeNetworks(X = mc, q = 0.9, nNet = 1L,
                                        nCells = NCELLS, nComp = 3L,
                                        scaleScores = TRUE, symmetric = FALSE,
                                        nCores = 8L))
  cat("     → 推算 nNet=10：", round(r3$secs * 10 / 60, 1), " 分钟\n", sep = "")

  r4 <- tmr("makeNetworks(nNet=3)",
            scTenifoldNet::makeNetworks(X = mc, q = 0.9, nNet = 3L,
                                        nCells = NCELLS, nComp = 3L,
                                        scaleScores = TRUE, symmetric = FALSE,
                                        nCores = 8L))
  r5 <- tmr("tensorDecomposition(K=3)",
            scTenifoldNet::tensorDecomposition(xList = r4$value, K = 3L,
                                               maxError = 1e-5, maxIter = 1000L,
                                               nDecimal = 3L))
  WT <- r5$value$X
  WT <- as.matrix(scTenifoldNet::strictDirection(WT, lambda = 0))
  diag(WT) <- 0; WT <- t(WT)
  cat("  WT 网络：", nrow(WT), " × ", ncol(WT),
      " | RUNX2 在内？", "RUNX2" %in% rownames(WT), "\n", sep = "")

  g <- intersect(c("RUNX2","COL1A1","CXCL12"), rownames(WT))[1]
  r6 <- tmr(paste0("单次敲除 ", g, " (MA + dRegulation)"), {
    KO <- WT; KO[g, ] <- 0
    set.seed(1)
    MA <- scTenifoldNet::manifoldAlignment(WT, KO, d = 2L, nCores = 8L)
    DR <- scTenifoldKnk::dRegulation(MA, empiricalNull = FALSE)
    DR
  })
  cat("     diffRegulation 列：", paste(colnames(r6$value), collapse = ", "), "\n")
  cat("     显著基因(adj<0.05)：", sum(r6$value$p.adj < 0.05), "/", nrow(r6$value), "\n")
  self <- r6$value[r6$value$gene == g, ]
  if (nrow(self)) cat("     自身：FC=", round(self$FC, 1), " p.adj=", signif(self$p.adj, 3), "\n")

  cat("\n  ── 推算（建网一次性 + 逐基因敲除）──\n")
  net10 <- r3$secs * 10
  cat("   建网(nNet=10)：", round(net10/60, 1), " 分钟（一次性）\n", sep = "")
  for (N in c(30, 60, 120)) {
    cat("   ", N, " 个基因串行：", round((net10 + r6$secs*N)/60, 1), " 分钟",
        " | 8 路并行约 ", round((net10 + r6$secs*N/6)/60, 1), " 分钟\n", sep = "")
  }
}
