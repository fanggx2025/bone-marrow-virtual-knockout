## =============================================================================
## 项目路径配置向导（在你的机器上跑一次即可）
##
## 作用：自动找到你的项目目录，把 config/config.R 里的路径改对，
##       并创建所有需要的目录结构。
##
## 用法（在你自己的 RStudio 里）：
##   source("C:/Users/Fang/WorkBuddy/2026-10-08-08-31-36/fracture_healing_vk/scripts/00_setup_paths.R")
##
## 注意：本脚本只创建 config/ R/ scripts/ data/ results/ logs/ 这些目录，
##       绝不会碰你已有的 Fig2 / Fig3 / Fig4 / Fig5 / Fig6。
## =============================================================================

Sys.setenv(LANG = "en_US.UTF-8", LC_ALL = "en_US.UTF-8")

cat("==================================================\n")
cat(" 项目路径配置向导\n")
cat("==================================================\n\n")

## ============================================================================
## 1. 自动探测项目根目录
## ============================================================================
cat("[1/4] 自动探测项目目录…\n")

## 关键教训：用户的「S3:」不是网络盘，而是桌面上的 S3 文件夹。
## R 会把 "S3:/xxx" 当成盘符 S3: 去解析，导致 setwd/file.copy 全部失败。
## 所以这里一律用「完整绝对路径」，不用盘符简写。
user <- Sys.getenv("USERPROFILE")

candidates <- c(
  ## 最可能：桌面 S3 文件夹（本次实际路径）
  file.path(user, "Desktop/S3/ciwujia_wuweizi_oa"),
  ## 常见网络映射盘
  "S:/ciwujia_wuweizi_oa",
  "S3:/ciwujia_wuweizi_oa",
  ## 桌面直接放
  file.path(user, "Desktop/ciwujia_wuweizi_oa"),
  file.path(user, "ciwujia_wuweizi_oa"),
  ## 文档
  file.path(user, "Documents/ciwujia_wuweizi_oa"),
  ## 工作区副本
  "C:/Users/Fang/WorkBuddy/2026-10-08-08-31-36/fracture_healing_vk"
)

cat("  检测候选路径：\n")
found <- NULL
for (p in candidates) {
  ok <- dir.exists(p)
  cat(sprintf("    %-52s %s\n", p, if (ok) "✓ 存在" else "✗"))
  if (ok && is.null(found)) found <- p
}

if (is.null(found)) {
  cat("\n✗ 所有候选路径都不存在。\n")
  cat("  请手动输入项目根目录的【完整绝对路径】，例如：\n")
  cat("    C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa\n")
  cat("  提示：在资源管理器里打开项目文件夹，点地址栏复制即可。\n\n")
  ans <- readline("  路径：")
  found <- trimws(ans)
  if (!nzchar(found) || !dir.exists(found)) {
    cat("✗ 路径无效：", found, "\n")
    quit(status = 1)
  }
}

##规范化反斜杠
PROJ_ROOT <- gsub("\\\\", "/", found)
if (!dir.exists(PROJ_ROOT)) {
  cat("✗ 路径不存在：", PROJ_ROOT, "\n")
  quit(status = 1)
}

cat("\n  ✓ 项目根目录：", PROJ_ROOT, "\n")

## ---- 安全检查：绝不覆盖用户已有内容 ----
cat("\n  检查已有内容（只读，不会修改）：\n")
existing_figs <- character(0)
for (d in c("Fig2", "Fig3", "Fig4", "Fig5", "Fig6")) {
  p <- file.path(PROJ_ROOT, d)
  ok <- dir.exists(p)
  cat(sprintf("    %-6s %s\n", d, if (ok) "✓ 存在（本脚本不会碰它）" else "○ 不存在"))
  if (ok) existing_figs <- c(existing_figs, d)
}
if (length(existing_figs) > 0) {
  cat("\n  ⚠ 检测到已有 Figure 目录：", paste(existing_figs, collapse = ", "), "\n")
  cat("    本脚本【不会读取、修改或覆盖】它们。\n")
  cat("    新代码框架会放在 ", PROJ_ROOT, "/ 的 config/ R/ scripts/ 下，\n", sep = "")
  cat("    与Figure 目录相互独立。\n")
}

## ============================================================================
## 2. 复制代码框架
## ============================================================================
cat("\n[2/4] 准备代码框架…\n")

WS <- "C:/Users/Fang/WorkBuddy/2026-10-08-08-31-36/fracture_healing_vk"

if (normalizePath(WS, winslash = "/") == normalizePath(PROJ_ROOT, winslash = "/")) {
  cat("  ✓ 代码框架已在项目目录中，无需复制\n")
  FRAMEWORK_READY <- TRUE
} else if (dir.exists(file.path(PROJ_ROOT, "R")) &&
           file.exists(file.path(PROJ_ROOT, "config", "config.R"))) {
  cat("  ✓ 代码框架已存在，跳过复制\n")
  FRAMEWORK_READY <- TRUE
} else {
  FRAMEWORK_READY <- FALSE
  if (!dir.exists(WS)) {
    cat("  ✗ 工作区代码框架不存在：", WS, "\n")
    quit(status = 1)
  }
  cat("  从工作区复制：", WS, "\n")
  cat("  → ", PROJ_ROOT, "\n\n")
  ans <- readline("  现在复制吗？(y/n) ")
  if (tolower(ans) != "y") {
    cat("  已取消。若要手动复制，把工作区的 config/ R/ scripts/ 拖到项目目录即可。\n")
    quit(status = 1)
  }
  for (d in c("config", "R", "scripts")) {
    src <- file.path(WS, d)
    dst <- file.path(PROJ_ROOT, d)
    if (!dir.exists(src)) next
    dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    fs <- list.files(src, full.names = TRUE)
    if (length(fs) > 0) {
      file.copy(fs, dst, overwrite = TRUE)
      cat("    ✓ ", d, "/ (", length(fs), " 个文件)\n", sep = "")
    }
  }
  for (f in c("README.md", "START_HERE.md", "RSTUDIO_GUIDE.md")) {
    s <- file.path(WS, f)
    if (file.exists(s)) file.copy(s, PROJ_ROOT, overwrite = TRUE)
  }
  cat("  ✓ 复制完成\n")
  FRAMEWORK_READY <- TRUE
}

## ============================================================================
## 3. 创建标准目录结构
## ============================================================================
cat("\n[3/4] 创建目录结构…\n")

dirs <- c("data/raw", "data/raw/GSE337700", "data/raw/ebi_fastq",
          "data/raw/tcmsp", "data/raw/targets", "data/raw/structures",
          "data/raw/cellranger",
          "data/processed", "results/figures", "results/tables",
          "results/objects", "results/ko_scTenifoldKnk",
          "logs", "notebooks")

n_created <- 0
for (d in dirs) {
  p <- file.path(PROJ_ROOT, d)
  if (!dir.exists(p)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    n_created <- n_created + 1
  }
}
cat("  ✓ 目录就绪（新建 ", n_created, " 个，共 ", length(dirs), " 个）\n", sep = "")

## ============================================================================
## 4. 写回 config.R
## ============================================================================
cat("\n[4/4] 更新 config/config.R…\n")

cfg <- file.path(PROJ_ROOT, "config", "config.R")
if (!file.exists(cfg)) {
  cat("  ✗ config.R 不存在。\n")
  quit(status = 1)
}

txt <- suppressWarnings(readLines(cfg, encoding = "UTF-8", warn = FALSE))

new_proj <- c(
  'PROJ <- list(',
  paste0('  root           = "', PROJ_ROOT, '",'),
  paste0('  data_raw       = "', PROJ_ROOT, '/data/raw",'),
  paste0('  data_processed = "', PROJ_ROOT, '/data/processed",'),
  paste0('  results        = "', PROJ_ROOT, '/results",'),
  paste0('  figures        = "', PROJ_ROOT, '/results/figures",'),
  paste0('  tables         = "', PROJ_ROOT, '/results/tables",'),
  paste0('  objects        = "', PROJ_ROOT, '/results/objects",'),
  paste0('  logs           = "', PROJ_ROOT, '/logs"'),
  ')'
)

start <- grep("^PROJ\\s*<-\\s*list\\(", txt)
if (length(start) == 0) {
  cat("  ✗ 未在 config.R 中找到 PROJ 定义\n")
  quit(status = 1)
}
start <- start[1]

end <- start
for (i in (start + 1):min(start + 15, length(txt))) {
  if (grepl("^\\s*\\)\\s*$", txt[i]) || grepl("^\\s*\\}\\s*$", txt[i])) {
    end <- i
    break
  }
}

txt_new <- c(
  txt[1:(start - 1)],
  new_proj,
  if (end < length(txt)) txt[(end + 1):length(txt)] else character(0)
)

## 备份后写入
## 注意：writeLines() 在 R 4.3 中不接受 encoding 参数（会报 unused argument）
backup <- paste0(cfg, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
file.copy(cfg, backup)

con <- file(cfg, open = "w", encoding = "UTF-8")
writeLines(txt_new, con, useBytes = FALSE)
close(con)

cat("  ✓ 已更新（原文件备份为 ", basename(backup), "）\n", sep = "")

## ============================================================================
## 验证
## ============================================================================
cat("\n==================================================\n")
cat(" 验证\n")
cat("==================================================\n\n")

## 在干净环境里source，确保能跑
rm(list = ls())
source(cfg)

cat("当前配置：\n")
cat("  root        : ", PROJ$root, "\n", sep = "")
cat("  data_raw    : ", PROJ$data_raw, "\n", sep = "")
cat("  results     : ", PROJ$results, "\n", sep = "")
cat("  SEED        : ", SEED, "\n", sep = "")
cat("  撞题回避基因: ", length(AVOID$genes), " 个\n", sep = "")
cat("  root_markers: ", paste(TRAJECTORY$root_markers, collapse = ", "), "\n", sep = "")

cat("\n写入测试：\n")
for (sub in c("logs", "data/raw")) {
  tf <- file.path(PROJ_ROOT, sub, "_test_write.txt")
  if (file.exists(tf)) file.remove(tf)
  ok <- tryCatch({ writeLines("ok", tf); TRUE }, error = function(e) FALSE,
                 warning = function(w) FALSE)
  if (ok && file.exists(tf)) {
    file.remove(tf)
    cat("  ✓ ", sub, "/ 可写\n", sep = "")
  } else {
    cat("  ✗ ", sub, "/ 无法写入\n", sep = "")
  }
}

cat("\n==================================================\n")
cat(" 配置完成\n")
cat("==================================================\n\n")
cat("下一步：\n")
cat("  1) 设工作目录：setwd(\"", PROJ_ROOT, "\")\n", sep = "")
cat("  2) 把已下载的数据解压到：\n     ", PROJ$data_raw, "/ebi_fastq\n", sep = "")
cat("     来源：C:/Users/Fang/WorkBuddy/2026-10-08-08-31-36/fracture_healing_vk/data/raw/ebi_fastq_backup.tar.gz\n")
cat("  3) 装 CellRanger 并跑 count（6-12 小时）\n")
cat("\n详细图文步骤见：RSTUDIO_GUIDE.md\n")
cat("==================================================\n")
