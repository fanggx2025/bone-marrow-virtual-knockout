## ------------------------------------------------------------------
## 97_verify_h5ad.R  验证 h5ad 文件是否完整（防截断）
## 用法：Rscript scripts/97_verify_h5ad.R
##      或在 RStudio 里 source("scripts/97_verify_h5ad.R")
## ------------------------------------------------------------------

Sys.setlocale("LC_CTYPE", "English_United States.utf8")
options(encoding = "UTF-8")
suppressMessages(suppressWarnings({
  library(rhdf5)
}))

## ---- 1. 定位项目根（ASCII 前缀模糊匹配，避中文路径硬编码） ----
find_root <- function() {
  cands <- c(
    "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa",
    "C:/Users/Fang/WorkBuddy/2026-10-08-08-31-36/fracture_healing_vk"
  )
  for (d in cands) if (dir.exists(d)) return(d)
  ## 兜底：扫 S3 下的中文目录
  base <- "C:/Users/Fang/Desktop/S3"
  if (dir.exists(base)) {
    subs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
    hit  <- subs[grepl("^\u523a", basename(subs))]  # 以"刺"开头
    if (length(hit)) return(hit[1])
  }
  getwd()
}

root <- find_root()
cat("项目根：", root, "\n\n")

fp <- file.path(root, "data", "raw", "bone_marrow_ref.h5ad")

if (!file.exists(fp)) {
  cat("!! 找不到文件：", fp, "\n")
  cat("   data/raw 现有：\n")
  print(list.files(file.path(root, "data", "raw")))
  quit(status = 1)
}

## ---- 2. 字节数对比 ----
sz  <- file.size(fp)
EXP <- 645050714  ## CELLxGENE 官方声明的字节数
cat(sprintf("本地大小  : %s 字节 (%.1f MB)\n", format(sz, big.mark = ","), sz / 1024^2))
cat(sprintf("官方声明  : %s 字节 (%.1f MB)\n", format(EXP, big.mark = ","), EXP / 1024^2))
cat(sprintf("差额      : %s 字节 (%.1f MB)\n\n",
            format(EXP - sz, big.mark = ","), (EXP - sz) / 1024^2))
if (sz < EXP) {
  cat(">> 大小不足，需确认是否截断。下面做结构自洽性检查...\n\n")
}

## ---- 3. 结构自洽性检查（这才是真正的判据） ----
cat("===== 结构自洽性检查 =====\n")
ok_all <- TRUE

f <- tryCatch(H5Fopen(fp), error = function(e) {
  cat("!! H5Fopen 失败：", conditionMessage(e), "\n"); NULL
})
if (is.null(f)) quit(status = 1)

grab <- function(path) tryCatch(h5read(f, path), error = function(e) {
  cat("   !! 读 ", path, " 失败：", conditionMessage(e), "\n", sep = "")
  ok_all <<- FALSE; NULL
})

## 顶层
cat("顶层对象：", paste(h5ls(f)$name[h5ls(f)$group == "/"], collapse = ", "), "\n\n")

## obs / var 行数
obs_idx <- grab("/obs/_index")
var_idx <- grab("/var/_index")
if (!is.null(obs_idx)) cat(sprintf("细胞数(/obs) : %s\n", format(length(obs_idx), big.mark = ",")))
if (!is.null(var_idx)) cat(sprintf("基因数(/var) : %s\n", format(length(var_idx), big.mark = ",")))

## X 的稀疏结构
indptr <- grab("/X/indptr")
if (!is.null(indptr)) {
  nr <- length(indptr) - 1L
  cat(sprintf("/X 推断行数 : %s  (来自 indptr 长度-1)\n", format(nr, big.mark = ",")))
  if (!is.null(obs_idx) && nr != length(obs_idx)) {
    cat("   !! 行数与 /obs 不一致 -> 文件大概率截断\n"); ok_all <- FALSE
  } else {
    cat("   OK 行数与 /obs 一致\n")
  }
  nnz <- indptr[length(indptr)]
  cat(sprintf("/X 非零元素 : %s\n", format(nnz, big.mark = ",")))
}

## 关键：末尾的数据能不能读出来（截断文件通常卡在这里）
data_tail <- tryCatch({
  h5read(f, "/X/data", index = list(NULL))  ## 全读，截断必炸
}, error = function(e) {
  cat("   !! 读 /X/data 失败：", conditionMessage(e), "\n")
  ok_all <<- FALSE; NULL
})
if (!is.null(data_tail)) {
  cat(sprintf("/X/data 长度: %s  (期望 = 非零元素数)\n", format(length(data_tail), big.mark = ",")))
  if (!is.null(indptr) && length(data_tail) != indptr[length(indptr)]) {
    cat("   !! /X/data 长度与 indptr 末尾不符 -> 截断\n"); ok_all <- FALSE
  } else cat("   OK /X/data 完整\n")
}

indices_ok <- tryCatch({
  idx <- h5read(f, "/X/indices")
  cat(sprintf("/X/indices 长度: %s\n", format(length(idx), big.mark = ",")))
  if (!is.null(indptr) && length(idx) != indptr[length(indptr)]) {
    cat("   !! indices 长度不符 -> 截断\n"); ok_all <<- FALSE
  } else cat("   OK /X/indices 完整\n")
  TRUE
}, error = function(e) {
  cat("   !! 读 /X/indices 失败：", conditionMessage(e), "\n"); ok_all <<- FALSE; FALSE
})

H5Fclose(f)

cat("\n===== 结论 =====\n")
if (ok_all) {
  cat("文件结构自洽，可以放心 load_h5ad()。\n")
} else {
  cat("!! 文件不完整，必须用 -C - 续传补全，别直接拿去分析。\n")
  cat("   续传命令：\n")
  cat("   curl -L -C - -o data/raw/bone_marrow_ref.h5ad <URL>\n")
}
