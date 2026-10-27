## =============================================================================
## P0: 环境体检脚本（先跑这个，再跑别的）
##
## 用途：一次性检查 R / Python / 磁盘 / 内存是否满足项目要求
## 运行：Rscript scripts/00_check_env.R
## 退出码：0 = 可继续；1 = 有硬性缺失
## =============================================================================

message("=================== P0 环境体检 ===================\n")

problems <- character(0)
warnings <- character(0)

## ---- 1. R 版本 ----
message("\n[1/7] R 环境")
r_version <- getRversion()
message("  R 版本：", getRversion())
if (r_version < "4.2.0") {
  problems <- c(problems, "R 版本过低（需 >= 4.2.0），Bioconductor 依赖装不上")
} else {
  message("  ✓ R 版本满足要求")
}

## ---- 1b. 库目录遮挡检测 ----
## 背景：包安装中途失败会留下 00LOCK-xxx 锁目录，或空包目录。
##       .libPaths() 里用户库排在前，空目录会把后面主库里完好的同名包遮住，
##       表现为 requireNamespace 报 FALSE / packageVersion 报"无此包"，
##       但 library(pkg) 反而能载入——极具迷惑性。这里主动检出。
message("\n[1b/7] 库目录健康检查")
libs <- .libPaths()
shadows <- character(0)
stale_locks <- character(0)

for (lb in libs) {
  if (!dir.exists(lb)) next
  entries <- list.dirs(lb, full.names = FALSE, recursive = FALSE)
  for (e in entries) {
    if (grepl("^00LOCK", e)) {
      stale_locks <- c(stale_locks, file.path(lb, e))
      next
    }
    ## 认定为可疑：目录存在但既无 DESCRIPTION 也无 Meta
    has_desc <- file.exists(file.path(lb, e, "DESCRIPTION"))
    has_meta <- file.exists(file.path(lb, e, "Meta", "package.rds"))
    if (!has_desc && !has_meta) {
      n <- length(list.files(file.path(lb, e), all.files = TRUE, no.. = TRUE))
      if (n == 0) shadows <- c(shadows, file.path(lb, e))
    }
  }
}

if (length(stale_locks) > 0) {
  message("  ⚠ 发现 ", length(stale_locks), " 个残留锁目录（安装中断遗留）:")
  for (s in stale_locks) message("      ", s)
  message("     修复：直接删除即可（主库里的包不受影响）")
  problems <- c(problems, paste0("残留锁目录需手动删除：", paste(basename(stale_locks), collapse = ", ")))
}

if (length(shadows) > 0) {
  message("  ✗ 发现 ", length(shadows), " 个空包目录，会遮挡后续库中的同名包:")
  for (s in shadows) message("      ", s)
  message("     修复：删除这些空目录")
  problems <- c(problems, paste0("空包目录需删除：", paste(basename(shadows), collapse = ", ")))
} else if (length(stale_locks) == 0) {
  message("  ✓ 库目录干净，无空目录/锁残影")
}

## ---- 2. 核心 R 包 ----
##必需（缺失则无法开始）
message("\n[2/7] 核心 R 包（必需）")
core_pkgs <- c("Seurat", "dplyr", "tidyr", "stringr",
               "ggplot2", "patchwork", "jsonlite")
for (p in core_pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  if (ok) {
    ver <- tryCatch(as.character(utils::packageVersion(p)),
                    error = function(e) "已装(版本读取异常)")
    message(sprintf("  ✓ %-14s %s", p, ver))
  } else {
    message(sprintf("  ✗ %-14s 未安装", p))
    problems <- c(problems, paste0("缺少 R 包：", p))
  }
}

## 可选包（缺失有降级方案，不算硬性问题）
message("\n[2b/7] 可选包（缺失不影响启动）")
opt_core <- c("SeuratDisk" = "h5Seurat 格式支持；本项目用标准 Seurat 流程",
              "GEOquery"  = "自动下载 GEO 数据；也可手动下载后跳过",
              "R.utils"   = "gzip 解压")
for (p in names(opt_core)) {
  ok <- requireNamespace(p, quietly = TRUE)
  message(sprintf("  %s %-14s %s", if (ok) "✓" else "○", p,
                  if (ok) as.character(utils::packageVersion(p))
                  else paste0("可选 — ", opt_core[[p]])))
}

## ---- 3. 分析专用 R 包 ----
## 分两类：必需（阻塞对应步骤）vs 可选（有降级方案）
message("\n[3/7] 分析专用 R 包")

## 可选包：缺失不阻塞，各自有替代方案
opt_pkgs <- c("monocle3", "CellRank", "SeuratDisk",
              "harmony", "limma", "clusterProfiler", "org.Hs.eg.db")
for (p in opt_pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  if (ok) {
    message(sprintf("  ✓ %-16s %s", p, as.character(utils::packageVersion(p))))
  } else {
    note <- switch(p,
      monocle3     = "  → 拟时序必需，见下方说明",
      CellRank     = "  → 可选，steer_score 的 fate 分量可改用轨迹位置代理",
      SeuratDisk   = "  → 可选，CRAN 已撤回，本项目用标准 Seurat 流程",
      "harmony"    = "  → 建议装，多样本整合效果明显优于 PCA 校正",
      "limma"      = "  → Bulk 验证必需",
      "clusterProfiler" = "  → 富集分析必需",
      "org.Hs.eg.db"    = "  → 基因ID 转换必需",
      "")
    message(sprintf("  ○ %-16s 待安装%s", p, note))
  }
}

## monocle3 特殊处理：不在 CRAN，且无 Windows 二进制
if (!requireNamespace("monocle3", quietly = TRUE)) {
  message("")
  message("  [monocle3 安装说明]")
  message("    · 不在 CRAN，须从 GitHub：")
  message("      remotes::install_github('cole-trapnell-lab/monocle3')")
  message("    · 依赖 igraph(已装) + Matrix(已装) + bpmnchess 等，通常可直接编译")
  message("    · 若编译失败：装 Rtools44 后重试")
  message("      https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html")
  message("    · 替代方案：用 DDRTree(已装) 或 Slingshot 走 Seurat wrappers")
}

## ---- 4. Python 环境 ----
message("\n[4/7] Python 环境（虚拟敲除必需）")
py_ok <- FALSE
py_version <- NA_character_
for (py in c("python", "python3")) {
  v <- tryCatch(
    suppressWarnings(system2(py, c("-c", shQuote(
      "import sys;print('.'.join(map(str,sys.version_info[:3])))")),
      stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(v) > 0 && grepl("^[0-9]", v[1])) {
    message("  ✓ ", py, " → Python ", v[1])
    py_version <- v[1]
    py_ok <- TRUE
    break
  }
}
if (!py_ok) {
  message("  ✗ 未找到 Python")
  problems <- c(problems, "未找到 Python（scTenifoldKnk 依赖）")
}

if (py_ok && compareVersion(paste_version <- py_version, "3.9.0") < 0) {
  message("  ⚠ Python ", py_version, " 偏低，scTenifoldNet 建议 3.9")
  warnings <- c(warnings, "Python 版本建议 3.9（当前 " + py_version + "）")
}

## scTenifoldNet 是虚拟敲除的核心包，单独确认
if (py_ok) {
  st_ok <- tryCatch({
    out <- suppressWarnings(system2("python",
      c("-c", shQuote("import scTenifoldNet;print('ok')")),
      stdout = TRUE, stderr = FALSE))
    any(grepl("ok", out))
  }, error = function(e) FALSE)

  if (st_ok) {
    message("  ✓ scTenifoldNet 已安装")
  } else {
    message("  ✗ scTenifoldNet 未安装")
    message("    安装：conda create -n vk python=3.9 && conda activate vk")
    message("           pip install scTenifoldNet")
    warnings <- c(warnings, "scTenifoldNet 未安装（Step 06 前必须装）")
  }
}

## ---- 5. 磁盘空间 ----
message("\n[5/7] 磁盘空间")
## Windows 用 PowerShell 查C 盘（disk.table 在 Windows 上不可靠）
ds <- tryCatch({
  if (.Platform$OS.type == "windows") {
    out <- suppressWarnings(system2("powershell.exe",
      c("-NoProfile", "-Command",
        "(Get-PSDrive C).Free/1GB"),
      stdout = TRUE, stderr = FALSE))
    free_gb <- as.numeric(gsub("[^0-9.]", "", out[1]))
    if (is.na(free_gb)) NA_real_ else free_gb
  } else {
    df <- utils::disk.table()
    sum(df$free * df$incr) / 1024^3
  }
}, error = function(e) NA_real_)

if (!is.na(ds)) {
  message(sprintf("  C盘可用：%.1f GB", ds))
  ## GSE337700约 1-5 GB（fastq 更小），解压 + 中间对象 + 虚拟敲除结果，预留 30GB
  if (ds < 30) {
    message("  ✗ 空间不足 30 GB")
    problems <- c(problems, "C 盘可用空间不足 30 GB")
  } else if (ds < 60) {
    message("  ⚠ 空间偏紧，虚拟敲除结果约需 20-40 GB")
    warnings <- c(warnings, "磁盘空间偏紧（虚拟敲除需 20-40 GB）")
  } else {
    message("  ✓ 空间充足")
  }
} else {
  message("  ○ 无法检测（可手动查看：C 盘剩余空间）")
}

## ---- 6. 内存 ----
message("\n[6/7] 内存")
mem_gb <- tryCatch({
  if (.Platform$OS.type == "windows") {
    out <- suppressWarnings(system2("powershell.exe",
      c("-NoProfile", "-Command",
        "[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1)"),
      stdout = TRUE, stderr = FALSE))
    as.numeric(gsub("[^0-9.]", "", out[1]))
  } else {
    NA_real_
  }
}, error = function(e) NA_real_)

if (!is.na(mem_gb)) {
  message(sprintf("  物理内存：%.1f GB", mem_gb))
  if (mem_gb < 16) {
    message("  ⚠ 内存 < 16 GB，虚拟敲除建议把KO_CONFIG$n_cells_per_sample 降到 300")
    warnings <- c(warnings, "内存 < 16 GB，建议 n_cells_per_sample 降到 300")
  } else {
    message("  ✓ 内存充足")
  }
} else {
  message("  ○ 无法检测")
}

## ---- 汇总 ----
message("\n=================== 体检汇总 ===================")
if (length(problems) == 0) {
  message("✓ 无硬性问题，可以开始 Step 01")
} else {
  message("✗ 存在", length(problems), "个必须先解决的问题：")
  for (i in seq_along(problems)) message("  [", i, "] ", problems[i])
}

if (length(warnings) > 0) {
  message("\n⚠ ", length(warnings), "个提醒：")
  for (i in seq_along(warnings)) message("  [", i, "] ", warnings[i])
}

message("\n下一步：")
if (length(problems) == 0) {
  message("  1) 手动下载 GSE337700 → data/raw/GSE337700/")
  message("  2) Rscript scripts/00_prepare_gse337700.R")
  message("  3) Rscript scripts/run_pipeline.R 01")
} else {
  message("  先解决上述硬性问题。安装命令见 README 第三节")
}
message("==================================================\n")

if (length(problems) > 0) quit(status = 1)
