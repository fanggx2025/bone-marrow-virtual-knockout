## 用已缓存的 WT 网络与敲除结果重新评分（改判据后无需重建网络/重跑敲除）
source("scripts/99_start.R")
source("R/04_virtual_knockout.R")
k <- readRDS(file.path(PROJ$objects, "s04_knockout.rds"))
sc <- score_ko(k$WT, k$dr)
cat("\n════════ 重新评分（秩和富集版）════════\n")
print(sc, row.names = FALSE)
write.csv(sc, file.path(PROJ$tables, "s04_knk_scores.csv"), row.names = FALSE)
saveRDS(list(WT = k$WT, dr = k$dr, scores = sc),
        file.path(PROJ$objects, "s04_knockout.rds"))
cat("\n✓ 已覆盖保存\n")
