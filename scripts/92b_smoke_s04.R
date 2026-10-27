## =============================================================================
## 92b_smoke_s04.R —— S04 冒烟测试（秒级，用极小矩阵跑通全链路）
## =============================================================================
## 【为什么需要】
##   建 WT 网络实测 18 分钟。任何一处小错都要赔上 18 分钟。
##   所以正式跑之前，先用 600 细胞 × 200 基因 × 2 张网把整条链路跑一遍，
##   确认 prepare_counts → build_wt → knock_many → score_ko 全都通。
## =============================================================================

source("scripts/99_start.R")
source("R/04_virtual_knockout.R")

set.seed(42)
obj <- readRDS(file.path(PROJ$objects, "s02_stromal.rds"))
obj <- subset(obj, cells = sample(colnames(obj), 600))
cat("冒烟用细胞：", ncol(obj), "\n")

force <- unique(c("RUNX2","SP7","PPARG","CEBPB","CXCL12","LPL","KITLG",
                  TRAJECTORY$end_markers, TRAJECTORY$root_markers,
                  TRAJECTORY$mid_markers))

m  <- prepare_counts(obj, hvg_n = 200L, force_genes = force)
WT <- build_wt(m, n_net = 2L, n_cells = 200L, n_workers = 2L,
               cache_nets = file.path(tempdir(), "smoke_nets.rds"))
cat("\nWT dim:", dim(WT), "\n")
cat("无边基因数：", length(attr(WT, "no_edge_genes")), "\n")

gs <- intersect(c("RUNX2","PPARG","CXCL12","LPL"), rownames(WT))
cat("可敲基因：", paste(gs, collapse = ", "), "\n")
dr <- knock_many(WT, gs, n_workers = 3L)
sc <- score_ko(WT, dr)

cat("\n──────── 冒烟结果 ────────\n")
print(sc, row.names = FALSE)

cat("\n──────── 方向核对 ────────\n")
for (g in gs) {
  em <- intersect(TRAJECTORY$end_markers, rownames(WT))
  w  <- WT[g, em]
  cat(sprintf("%-8s → 成骨标记边：%s\n", g,
              paste(sprintf("%s=%+.4f", names(w), w), collapse = " ")))
}
cat("\n✓ 冒烟测试通过（若上面没报错）\n")
