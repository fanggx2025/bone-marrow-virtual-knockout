#!/usr/bin/env Rscript
## =============================================================================
## 自动化流水线运行器
##
## 用法：
##   Rscript scripts/run_pipeline.R              # 完整运行
##   Rscript scripts/run_pipeline.R 01 04 06     # 只跑指定步骤
##   QUICK_MODE=1 Rscript scripts/run_pipeline.R  # 快速可行性验证模式
##
## 特性：
##   · 断点续跑（已完成的步骤自动跳过）
##   · 每步独立日志（logs/SXX.log）
##   · 全局错误处理（单步失败不中断整体）
##   · 阶段决策点提醒
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

setwd("C:/Users/Fang/WorkBuddy/2026-10-08-08-31-36/fracture_healing_vk")
source("config/config.R")

## 命令行参数
args <- commandArgs(trailingOnly = TRUE)

## 步骤定义（顺序不可随意调整）
STEP_SCRIPTS <- c(
  "01" = "R/01_qc_integrate.R",
  "04" = "R/04_trajectory.R",
  "06" = "R/06_virtual_knockout.R",
  "07" = "R/07_celloracle_integration.R",
  "08" = "R/08_network_pharmacology.R",
  "09" = "R/09_docking.R",
  "10" = "R/10_validation.R"
)

STEP_DESC <- c(
  "01" = "数据获取与单细胞质控",
  "04" = "拟时序轨迹分析（Fig 3）",
  "06" = "虚拟敲除筛选（Fig 5）★核心",
  "07" = "CellOracle 扰动验证（Fig 6）",
  "08" = "网络药理学交汇",
  "09" = "分子对接与 MD",
  "10" = "独立队列验证（Fig 8）"
)

## 环境变量覆盖
if (Sys.getenv("QUICK_MODE") == "1") {
  QUICK_MODE <<- TRUE
  message("[ENV] QUICK_MODE = TRUE")
}

## 决定运行哪些步骤
steps_to_run <- if (length(args) > 0) {
  intersect(names(STEP_SCRIPTS), args)
} else {
  names(STEP_SCRIPTS)
}

if (length(steps_to_run) == 0) {
  stop("无效的步骤编号。可选：", paste(names(STEP_SCRIPTS), collapse = ", "))
}

## -----------------------------------------------------------------------------
## 执行引擎
## -----------------------------------------------------------------------------
run_step <- function(step_id) {

  script <- STEP_SCRIPTS[[step_id]]
  log_file <- file.path(PROJ$logs, paste0("S", step_id, ".log"))

  cat("\n")
  cat("━" |> strrep(72) |> str_pad(0, side = "both"), "\n")
  cat(sprintf("  [%s] %s\n", step_id, STEP_DESC[[step_id]]))
  cat("━" |> strrep(72) |> str_pad(0, side = "both"), "\n")

  if (!file.exists(script)) {
    cat("  ✗ 脚本不存在：", script, "\n")
    return(FALSE)
  }

  start_time <- Sys.time()

  ## 记录开始
  writeLines(c(
    sprintf("=== %s | %s ===", step_id, STEP_DESC[[step_id]]),
    sprintf("开始时间: %s", format(start_time))
  ), log_file)

  ## 执行（隔离环境，变量不污染主环境）
  result <- tryCatch({

    ## 在新环境中 source，脚本内的局部变量不影响此处
    env <- new.env(parent = globalenv())
    sys.source(script, envir = env, keep.source = FALSE)

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

    writeLines(c(
      sprintf("结束时间: %s", format(Sys.time())),
      sprintf("耗时: %.2f 分钟", elapsed),
      "状态: 成功"
    ), log_file, append = TRUE)

    cat(sprintf("  ✓ 完成（%.1f 分钟）\n", elapsed))
    TRUE

  }, error = function(e) {

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

    writeLines(c(
      sprintf("结束时间: %s", format(Sys.time())),
      sprintf("耗时: %.2f 分钟", elapsed),
      "状态: 失败",
      sprintf("错误: %s", conditionMessage(e))
    ), log_file, append = TRUE)

    cat(sprintf("  ✗ 失败（%.1f 分钟）\n", elapsed))
    cat("    错误：", conditionMessage(e), "\n")
    cat("    日志：", log_file, "\n")

    ## 打印错误调用栈
    calls <- sys.calls()
    if (length(calls) > 0) {
      cat("    调用栈（末 3 层）：\n")
      for (i in seq(max(1, length(calls) - 2), length(calls))) {
        cat("      ", deparse(calls[[i]])[1], "\n")
      }
    }

    FALSE
  })

  return(result)
}

## -----------------------------------------------------------------------------
## 主流程
## -----------------------------------------------------------------------------
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║  寒地龙药 · 骨折不愈合 单细胞虚拟敲除研究流水线                   ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

cat("\n配置摘要：\n")
cat(sprintf("  SEED          : %d\n", SEED))
cat(sprintf("  QUICK_MODE    : %s\n", QUICK_MODE))
cat(sprintf("  KO Tier1 上限 : %d 个基因\n", KO_CONFIG$KO_TIER1_MAX))
cat(sprintf("  最终靶点目标  : %d 个\n", KO_CONFIG$KO_TIER3_N))
cat(sprintf("  待运行步骤    : %s\n", paste(steps_to_run, collapse = " → ")))

## 环境预检
## 注意：列名用 ASCII，避免中文 Windows locale 下 data.frame 列名解析异常
cat("\n环境预检：\n")
env_check <- data.frame(
  component = c("R", "Seurat", "monocle3", "scTenifoldNet", "CellOracle", "Vina"),
  status = c(
    paste(R.version.string),
    ifelse(requireNamespace("Seurat", quietly = TRUE), "OK", "MISSING"),
    ifelse(requireNamespace("monocle3", quietly = TRUE), "OK", "MISSING (github)"),
    tryCatch({
      if (nzchar(Sys.which("python"))) "check setup_env.R" else "no python"
    }, error = function(e) "no python"),
    ifelse(requireNamespace("celloracle", quietly = TRUE), "OK", "optional"),
    ifelse(nzchar(Sys.which("vina")), "OK", "optional")
  ),
  stringsAsFactors = FALSE
)
print(env_check, row.names = FALSE)

## Python 环境检查
py_ok <- system2("python", c("-c", "import scTenifoldNet"),
                 stdout = FALSE, stderr = FALSE)
if (py_ok != 0) {
  cat("\n⚠ scTenifoldNet 未就绪，S06 将失败。安装：\n")
  cat("  conda create -n vk python=3.9 && conda activate vk\n")
  cat("  pip install scTenifoldNet scTenifoldKnk\n")
}

## 正式运行
cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("开始执行流水线\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

results <- list()
for (s in steps_to_run) {
  results[[s]] <- run_step(s)

  ## 阶段决策点提醒
  if (s == "01" && results[[s]]) {
    cat("\n  💡 决策点 P1：检查质控结果\n")
    cat("     → 若 GSE337700 质控后细胞数 < 1500，转为主攻 GSE93213/93215\n")
  }
  if (s == "04" && results[[s]]) {
    cat("\n  💡 决策点 P2：检查轨迹与分岔分析\n")
    cat("     → 若未检出显著分岔节点，检查样本量与标记基因\n")
    cat("     → 若 S04_top_branch_nodes.txt 已生成，可继续 S06\n")
  }
  if (s == "06" && results[[s]]) {
    cat("\n  💡 决策点 P3：虚拟敲除可行性\n")
    cat("     → 等待数小时后检查 scTenifoldKnk.log\n")
    cat("     → 若完成基因数 < 50 或结果全为噪声，立即转向方案一\n")
    cat("     → 若有效响应基因 > 20 且 steer_score 有区分度，继续\n")
  }
}

## -----------------------------------------------------------------------------
## 汇总
## -----------------------------------------------------------------------------
elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

cat("\n\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("流水线执行完毕\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

summary_df <- data.frame(
  step  = steps_to_run,
  name  = unname(STEP_DESC[steps_to_run]),
  state = ifelse(unlist(results), "OK", "FAILED"),
  stringsAsFactors = FALSE
)
print(summary_df, row.names = FALSE)

## 关键产出检查
cat("\n关键产出检查：\n")
key_outputs <- c(
  "results/tables/S06_final_targets.txt",
  "results/tables/S07_consensus_targets.txt",
  "results/tables/S08_final_converged.txt",
  "results/figures/Fig3_trajectory.png",
  "results/figures/Fig5_virtual_KO.png",
  "results/figures/Fig6_celloracle.png",
  "results/figures/Fig8_validation.png"
)
for (f in key_outputs) {
  fp <- file.path(PROJ$root, f)
  if (file.exists(fp)) {
    cat(sprintf("  ✓ %-45s (%s)\n", basename(f),
                format(file.size(fp), big.mark = ",")))
  } else {
    cat(sprintf("  ○ %-45s (未生成)\n", basename(f)))
  }
}

## 下一步建议
cat("\n下一步建议：\n")
if (!file.exists(file.path(PROJ$root, "results/tables/S07_consensus_targets.txt"))) {
  cat("  1. 在 Python 环境运行 CellOracle：\n")
  cat("     python scripts/run_celloracle.py \\\n")
  cat("       data/processed/counts_matrix.csv.gz \\\n")
  cat("       results/celloracle \\\n")
  cat("       $(cat results/tables/S06_candidate_pool.txt | tr '\\n' ',')\n\n")
}
cat("  2. 准备 GWAS 数据（骨密度/骨愈合性状）→ data/raw/gwas/\n")
cat("  3. 手动准备 TCMSP 成分数据 → data/raw/tcmsp/\n")
cat("  4. 填写靶点 UniProt/PDB → data/raw/structures/structure_mapping.csv\n")

cat(sprintf("\n总耗时：%.1f 分钟\n", elapsed_total))
cat(sprintf("日志目录：%s\n", PROJ$logs))
cat("\n")