## 95_audit_markers.R
## 验证两件事：
##   1. 数据集的基因命名是人类还是小鼠（决定 MARKERS 能不能直接用）
##   2. config.R 的 MARKERS（小鼠命名）与 HUMAN_MARKERS 各自的命中率
Sys.setlocale("LC_CTYPE", "English_United States.utf8")
options(encoding = "UTF-8")
suppressMessages(suppressWarnings({ library(rhdf5) }))

f <- H5Fopen("data/raw/bone_marrow_ref.h5ad")
gs <- as.character(h5read(f, "/raw/var/gene_symbols"))
H5Fclose(f)

cat("基因总数:", length(gs), "\n")
cat("样例(1:15):", paste(head(gs, 15), collapse = ", "), "\n\n")

## 大小写统计：人类基因符号几乎全大写，小鼠是首字母大写
upper_only <- grepl("^[A-Z0-9.\\-]+$", gs)
cat("全大写（人类风格）比例:", sprintf("%.1f%%", 100 * mean(upper_only)), "\n")
cat("含小写（小鼠风格）比例  :", sprintf("%.1f%%", 100 * mean(!upper_only)), "\n\n")

## 用几组明确的同源基因对来判定物种
pairs <- list(
  c("RUNX2", "Runx2"), c("COL1A1", "Col1a1"), c("CD3E", "Cd3e"),
  c("PTPRC", "Ptprc"), c("SOX9", "Sox9"), c("SP7", "Sp7"),
  c("ALPL", "Alpl"), c("MMP13", "Mmp13"), c("CXCL12", "Cxcl12")
)
cat("===== 人鼠命名判定 =====\n")
cat(sprintf("%-10s %-10s %-8s %-8s\n", "人类符号", "小鼠符号", "人类在?", "小鼠在?"))
h_cnt <- 0; m_cnt <- 0
for (p in pairs) {
  h <- p[1] %in% gs; m <- p[2] %in% gs
  if (h) h_cnt <- h_cnt + 1
  if (m) m_cnt <- m_cnt + 1
  cat(sprintf("%-10s %-10s %-8s %-8s\n", p[1], p[2],
              ifelse(h, "✓", "-"), ifelse(m, "✓", "-")))
}
cat(sprintf("\n人类符号命中 %d/%d，小鼠符号命中 %d/%d\n",
            h_cnt, length(pairs), m_cnt, length(pairs)))
cat("结论:",
    if (h_cnt > m_cnt) "★ 人类命名 —— config.R 的 MARKERS（小鼠）不能直接用！"
    else "小鼠命名", "\n\n")

## 检查 config.R MARKERS 的实际命中率
cat("===== config.R MARKERS（小鼠命名）在人类数据中的命中率 =====\n")
sys.source("config/config.R", envir = e <- new.env())
ML <- get("MARKERS", envir = e)
for (ct in names(ML)) {
  g <- ML[[ct]]
  hit <- intersect(g, gs)
  pct <- length(hit) / length(g)
  cat(sprintf("  %-16s %2d/%2d (%.0f%%) %s\n", ct, length(hit), length(g), pct*100,
              if (pct == 0) "← 完全落空" else if (pct < 0.5) "← 命中过低" else ""))
}
tot <- unlist(ML); cat(sprintf("\n总体命中: %d/%d = %.1f%%\n",
    sum(tot %in% gs), length(tot), 100*mean(tot %in% gs)))
