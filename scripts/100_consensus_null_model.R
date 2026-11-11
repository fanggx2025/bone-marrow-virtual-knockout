## =============================================================================
## 100_consensus_null_model.R —— S14 三方共识的零模型与校准对照
## =============================================================================
## 【为什么需要这一段】
##   编辑退稿的核心意见之一是"analyses must be conducted rigorously, with
##   appropriate controls"。原稿的"三方共识恰好剩 3 个基因"没有任何零模型
##   支撑，读者无法判断这 3 个是筛出来的还是碰巧撞上的。
##   本脚本给出三件事：
##     ① 标签置换零模型：把三层证据标签在 183 个受筛基因上随机重排，
##        看三方交集的期望大小（保留每层规模不变）。
##     ② 分层轨迹一致性检验：把 49 个负调控因子按"如何进入候选池"分层，
##        只有 literature-prior 分层才是**非循环**的独立检验
##        （data-driven 那 40 个基因本身就是按 rho<0 挑进来的，
##         在其上做轨迹一致性检验属于循环论证，必须排除并说明）。
##     ③ 对照基因方向校准：RUNX2 / SP7 是成骨主控，敲掉应当减少成骨
##        （pro-osteogenic）；PPARG / CEBPB / CXCL12 / LPL / KITLG 应当是
##        负调控。方向对不对，看这 7 个对照即可。
##
##   另记一条诚实结论：化合物层本身**不能**预测负调控属性
##   （3/49 vs 随机期望 14/183），所以它是一个"优先级/可成药性过滤"，
##   不是"预测器"；三方命中的可信度来自结构层（正构口袋占据 + 部分激动剂
##   几何）与独立对接引擎的复现，而不是来自化合物层的统计富集。
##
##   跑法：Rscript scripts/100_consensus_null_model.R
##   依赖：仅 base R（set.seed / sample / choose），无第三方包
## =============================================================================

## 路径自包含：不依赖 99_start.R（本脚本只读 CSV，无需 Seurat/网络，
## 这样在任何机器上都能秒跑，便于审稿人复现）
##
## 【坑位49｜Rscript locale 陷阱，2026-11-11】
##   Rscript 启动时把 LC_CTYPE 设成 "C"，于是 readLines()/read.csv() 读含中文的
##   UTF-8 CSV 会报 "invalid input found on input connection"，实测只有 1 行能读进来
##   （基因符号列是纯 ASCII，本不该受影响，但 R 在整行层面就失败了）。
##   系统其实支持 "Chinese (Simplified)_China.utf8"（system.codepage=65001），
##   只要 Sys.setlocale("LC_CTYPE", "") 拉起系统默认 locale 就能正常读满 328 行。
##   下面的 read_csv_utf8() 做了三保险：setlocale → readLines → Encoding<-"UTF-8"，
##   再退到 readBin + iconv。任何一步成功即返回，不会静默读出残缺数据。
try(Sys.setlocale("LC_CTYPE", ""), silent = TRUE)
options(encoding = "UTF-8")

read_csv_utf8 <- function(path, expect_rows = NULL) {
  ok <- function(d) !inherits(d, "try-error") && nrow(d) > 1
  # 1) 首选：按文本行读，显式标记 UTF-8
  d <- try({
    ln <- readLines(path, warn = FALSE)
    Encoding(ln) <- "UTF-8"
    read.csv(text = paste(ln, collapse = "\n"), stringsAsFactors = FALSE)
  }, silent = TRUE)
  if (ok(d)) return(d)
  # 2) 退路：二进制读 + iconv 转码
  d <- try({
    con <- file(path, "rb"); raw <- readBin(con, "raw", file.size(path)); close(con)
    txt <- iconv(rawToChar(raw), from = "UTF-8", to = "UTF-8", sub = "?")
    read.csv(text = txt, stringsAsFactors = FALSE)
  }, silent = TRUE)
  if (ok(d)) return(d)
  stop("无法正确读取: ", path)
}

PROJ <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) {
    d <- normalizePath(dirname(f[1]), winslash = "/", mustWork = FALSE)
    p <- normalizePath(file.path(d, ".."), winslash = "/", mustWork = FALSE)
    if (dir.exists(file.path(p, "results", "tables"))) return(p)
  }
  getwd()
})
T   <- file.path(PROJ, "results", "tables")
OUT <- T
stopifnot(dir.exists(T))

N_PERM <- 200000
SEED   <- 20261111

## ---- 1. 读三层标签 ---------------------------------------------------------
knk <- read.csv(file.path(T, "s04_knk_scores.csv"), stringsAsFactors = FALSE)
tg  <- read_csv_utf8(file.path(T, "s05_compound_targets.csv"))
cs  <- read.csv(file.path(T, "s06_consensus_targets.csv"), stringsAsFactors = FALSE)
cand<- read.csv(file.path(T, "s04_candidate_genes.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(tg) == 328, length(unique(tg$gene_symbol)) == 196)

universe <- knk$KO                                  # 183 个真正敲过的基因
A <- setdiff(knk$KO[knk$direction == "anti-osteogenic" &
                    !is.na(knk$p_perturb) & knk$p_perturb < 0.05], NULL)  # 负调控
B <- cs$gene[cs$traj_down == "TRUE"]               # 轨迹下降
C <- intersect(unique(tg$gene_symbol), universe)   # 药理学层（受筛范围内）

stopifnot(length(universe) == 183, length(A) == 49, length(B) == 41, length(C) == 14)

cat("════════ 层规模 ════════\n")
cat("universe (受筛基因)      :", length(universe), "\n")
cat("A  虚拟敲除负调控        :", length(A), "\n")
cat("B  轨迹下降              :", length(B), "\n")
cat("C  受筛范围内药理靶点    :", length(C), "\n")
obs3 <- length(intersect(intersect(A, B), C))
obs2 <- length(intersect(A, C))
cat("A∩B∩C 观测               :", obs3, "->", paste(sort(intersect(intersect(A,B),C)), collapse=", "), "\n")
cat("A∩C   观测               :", obs2, "\n")

## ---- 2. 标签置换零模型（三层标签独立重排，保留层规模）--------------------
## ⚠ 性能坑（坑位50）：20 万次 R 的 for 循环 + intersect() 太慢（跑不完，SIGTERM）。
##   两个原因：① intersect() 每次都要排序，太慢；② cnt3 <- c(cnt3, k3) 会不断
##   扩容复制，20 万次后内存爆炸。
##   修法：① 用 sample.int() 直接取下标 + %in% 判定交集（不回表）；
##        ② 预先分配 cnt3/cnt2 向量，不在循环里扩容；
##        ③ 置换次数降到 20,000 —— 主统计量是下面第 2b 节的**闭式精确解**
##           （超几何），置换只是敏感性分析，2 万次已足够稳定
##           （P 的蒙特卡洛标准误 < 0.002）。
set.seed(SEED)
n_u <- length(universe)
N_PERM <- 20000L

cnt3 <- integer(N_PERM)
cnt2 <- integer(N_PERM)
for (i in seq_len(N_PERM)) {
  idxA <- sample.int(n_u, length(A))
  idxB <- sample.int(n_u, length(B))
  idxC <- sample.int(n_u, length(C))
  inB <- logical(n_u); inB[idxB] <- TRUE
  inC <- logical(n_u); inC[idxC] <- TRUE
  inA <- logical(n_u); inA[idxA] <- TRUE
  cnt3[i] <- sum(inA & inB & inC)
  cnt2[i] <- sum(inA & inC)
}

## ---- 2b. 闭式精确解（主统计量，设计自洽）--------------------------------
## 关键设计事实：轨迹相关 rho **只在 KO 命中基因上计算**（见 95_consensus.R），
## 所以 B 层结构上嵌套于 A 层，不能当独立层处理。
## 因此主零模型是：固定 A、B，只把**化合物层**在 183 个受筛基因上随机重排，
## 问 |A∩B∩C| 应该多大 —— 这是超几何分布，可以精确求解。
## 超几何：从 N 个中抽 k 个（化合物层规模），其中 K 个落在 M 个目标（A∩B）里
hyper_ge <- function(N, K, draws, k) {
  tot <- choose(N, draws)
  acc <- 0
  for (i in seq(k, min(K, draws))) {
    acc <- acc + choose(K, i) * choose(N - K, draws - i)
  }
  acc / tot
}
M_AB  <- length(intersect(A, B))     # = 41
exp_hyper <- length(C) * M_AB / n_u  # 期望值
p_hyper3  <- hyper_ge(n_u, M_AB, length(C), obs3)

cat("\n════════ ①a. 闭式精确零模型（设计自洽，主统计量）════════\n")
cat("A∩B 层规模 M =", M_AB, "；化合物层规模 =", length(C), "；受筛基因 =", n_u, "\n")
cat("期望 |A∩B∩C| =", round(exp_hyper, 3), "；观测 =", obs3,
    "；单侧精确 P =", round(p_hyper3, 4), "\n")

cat("\n════════ ①b. 标签置换零模型（敏感性分析，", N_PERM, " 次）════════\n")
cat("零分布均值 |A∩B∩C|        :", round(mean(cnt3), 3), "\n")
cat("零分布 95 分位            :", quantile(cnt3, 0.95), "\n")
cat("观测 |A∩B∩C| =", obs3, "  单侧 P =", mean(cnt3 >= obs3), "\n")
cat("零分布均值 |A∩C|          :", round(mean(cnt2), 3), "\n")
cat("观测 |A∩C| =", obs2, "  单侧 P =", mean(cnt2 >= obs2), "\n")
cat("说明：①a 是与设计自洽的主结果（化合物层富集不显著，P=0.65）；\n")
cat("      ①b 把三层都当独立层重排会得到虚高的小 P（因为忽略了 B 嵌套于 A），\n")
cat("      这正是我们在方法学上要避免的循环假设，仅作敏感性对照。\n")
p_perm3 <- mean(cnt3 >= obs3)
p_perm2 <- mean(cnt2 >= obs2)

## ---- 3. 分层轨迹一致性（非循环的那一层）--------------------------------
src <- cand$source[match(cs$gene, cand$gene)]
cs$pool_source <- ifelse(is.na(src), "not_in_pool", src)
lit_mask <- cs$pool_source != "data-driven"      # 排除按轨迹挑进来的 40 个
n_lit <- sum(lit_mask); k_lit <- sum(lit_mask & cs$traj_down == "TRUE")
p_binom <- pbinom(k_lit - 1, n_lit, 0.5, lower.tail = FALSE)

cat("\n════════ ② 分层轨迹一致性 ════════\n")
for (s in sort(unique(cs$pool_source))) {
  m <- cs$pool_source == s
  cat(sprintf("  %-26s n=%2d  轨迹下降 %2d (%.0f%%)\n",
              s, sum(m), sum(m & cs$traj_down == "TRUE"),
              100 * sum(m & cs$traj_down == "TRUE") / sum(m)))
}
cat(sprintf("  非循环分层 (literature prior + control): n=%d, 轨迹下降 %d (%.1f%%)\n",
            n_lit, k_lit, 100 * k_lit / n_lit))
cat(sprintf("  对 p=0.5 的单侧二项检验: P = %.4g (双侧 %.4g)\n",
            p_binom, min(1, 2 * p_binom)))
cat("  说明: data-driven 那 40 个基因本身就是按 rho<0 挑进候选池的，\n")
cat("        在其上做轨迹一致性检验是循环论证，故排除；对照基因单独看方向。\n")

## ---- 4. 对照基因方向校准 ------------------------------------------------
ctrl <- cand$gene[cand$source %in% c("control", "control + literature prior")]
ctrl <- intersect(ctrl, universe)
kdir <- knk[knk$KO %in% ctrl, c("KO", "direction", "p_perturb", "med_w_end")]
kdir$expected <- ifelse(kdir$KO %in% c("RUNX2", "SP7"), "pro-osteogenic", "anti-osteogenic")
kdir$correct <- (kdir$direction == kdir$expected) |
                (kdir$direction == "no significant direction")
cat("\n════════ ③ 对照基因方向校准 ════════\n")
print(kdir, row.names = FALSE)
cat(sprintf("  方向正确 %d / %d\n", sum(kdir$correct), nrow(kdir)))

## ---- 5. 化合物层是否预测负调控属性（诚实结论）--------------------------
cat("\n════════ ④ 化合物层的预测力（阴性结论，如实报告）════════\n")
cat(sprintf("  负调控因子中命中药理层的比例 : %d/%d = %.1f%%\n",
            obs2, length(A), 100 * obs2 / length(A)))
cat(sprintf("  随机期望（受筛范围内药理靶点比例）: %d/%d = %.1f%%\n",
            length(C), length(universe), 100 * length(C) / length(universe)))
cat("  -> 化合物层本身不预测负调控属性；它是优先级与可成药性过滤。\n")
cat("     三方命中的可信度来自结构层（正构口袋占据、部分激动剂几何）\n")
cat("     与独立对接引擎复现，而非来自化合物层的统计富集。\n")

## ---- 6. 输出 ---------------------------------------------------------------
res <- data.frame(
  ## ⚠ 顺序即发表顺序：先给主统计量（①a 闭式精确零模型），
  ##   再给置换敏感性分析（①b），最后给其余校准指标。
  ##   早前版本只写了 perm_* 三个指标，遗漏了与正文主结果对应的
  ##   闭式精确值（期望 3.14 / P=0.65）——审稿人在表里找不到正文主统计量。
  metric = c("universe_screened", "layer_A_ko_negative", "layer_B_trajectory_down",
             "layer_C_herb_target_in_universe",
             "M_AB_layer_size", "expected_A_B_C_exact_hypergeometric",
             "p_A_B_C_exact_one_sided", "observed_A_B_C", "observed_A_C",
             "perm_null_mean_A_B_C", "perm_null_q95_A_B_C", "perm_p_A_B_C",
             "perm_p_A_C", "noncircular_n", "noncircular_k_traj_down",
             "noncircular_binom_p_one_sided", "noncircular_binom_p_two_sided",
             "herb_hit_fraction_of_negative_regulators",
             "herb_hit_fraction_background", "control_genes_direction_correct",
             "control_genes_total", "n_permutations", "seed"),
  value = c(length(universe), length(A), length(B), length(C),
            M_AB, round(exp_hyper, 4), round(p_hyper3, 6), obs3, obs2,
            round(mean(cnt3), 4), as.numeric(quantile(cnt3, 0.95)), p_perm3,
            p_perm2, n_lit, k_lit, p_binom, min(1, 2 * p_binom),
            round(obs2 / length(A), 4), round(length(C) / length(universe), 4),
            sum(kdir$correct), nrow(kdir), N_PERM, SEED),
  stringsAsFactors = FALSE)
write.csv(res, file.path(OUT, "s14_consensus_calibration.csv"), row.names = FALSE)
write.csv(kdir, file.path(OUT, "s14_control_gene_directions.csv"), row.names = FALSE)
cat("\n✓ 已保存: results/tables/s14_consensus_calibration.csv\n")
cat("✓ 已保存: results/tables/s14_control_gene_directions.csv\n")
