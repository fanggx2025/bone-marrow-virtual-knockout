## =============================================================================
## 94_run_s04_full.R —— 全量候选基因虚拟敲除（194 个）
## =============================================================================
## 【坑位43：WT 网络的基因集必须包含全部候选基因】
##   试点阶段建的 WT 只塞了 7 个对照 + 成骨/轨迹标记（18 个强制基因），
##   194 个候选基因根本不在那张网络里 —— 直接复用缓存会导致
##   `intersect(ko_genes, rownames(WT))` 把它们全滤掉，只剩 20 来个。
##   ⇒ 全量跑之前必须**按含全部候选的基因集重建一次网络**（约 20 分钟，一次性）。
##
## 【规模】p = 1000 = 强制基因（194 候选 + 16 读出标记）+ 高变基因补满
##   单基因敲除约 1 s，194 个基因 10 进程并行 → 1 分钟内。
## =============================================================================

source("scripts/99_start.R")
source("R/04_virtual_knockout.R")

cd <- readRDS(file.path(PROJ$objects, "s04_candidates.rds"))
cat("候选基因：", nrow(cd$tab), " 个\n")

## WT 网络已按「含全部候选基因」的基因集建好并缓存（s04_wt_network.rds），
## 这里直接复用，只跑敲除（约 1 分钟）。
res <- run_s04(ko_genes = cd$tab$gene, hvg_n = 1000L,
               n_workers_ko = 10L,
               rebuild = FALSE, rebuild_nets = FALSE)

sc <- res$scores
cat("\n\n════════════════════════════════════════════\n")
cat("全量结果\n")
cat("════════════════════════════════════════════\n")
cat("判读分布：\n"); print(table(sc$direction))
cat("\nverdict 分布：\n"); print(table(sc$verdict))

hit <- sc %>% filter(direction == "抑制成骨" & p_perturb < 0.05) %>% arrange(p_perturb)
cat("\n──────── ★★ 负调控因子（可药靶，按模块扰动显著性排序）────────\n")
print(hit %>% select(KO, n_top10, pct_rank_end, p_perturb,
                     med_w_end, p_repress, med_w_root, verdict),
      row.names = FALSE)
cat("\n共", nrow(hit), "个 ★★\n")

req <- sc %>% filter(direction == "促进成骨" & p_perturb < 0.05) %>% arrange(p_perturb)
cat("\n──────── 成骨所需（敲除显著阻断成骨，可作反向对照）────────\n")
print(req %>% select(KO, n_top10, pct_rank_end, p_perturb, med_w_end, p_activate),
      row.names = FALSE)

cat("\n──────── 敲除无效（共识网络中无边，需剔除）────────\n")
ne <- sc %>% filter(!has_edges)
cat(nrow(ne), " 个：", paste(ne$KO, collapse = ", "), "\n")

write.csv(sc, file.path(PROJ$tables, "s04_knk_scores_full.csv"), row.names = FALSE)
cat("\n✓ 已保存：results/tables/s04_knk_scores_full.csv\n")
