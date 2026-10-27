## =============================================================================
## 03_search_cellxgene.R —— 在 CELLxGENE 里搜索数据集（纯 R，无需 Python）
##
## CELLxGENE 有 2238 个数据集，全部提供【已处理好的】H5AD，
## 不需要 CellRanger / kallisto 等任何定量工具（这些在 Windows 上都不可用）。
##
## 用法（RStudio）：
##   source("scripts/03_search_cellxgene.R")
##
##   # 搜骨相关
##   res <- cxg_search(c("bone","osteoblast","mesenchymal"))
##
##   # 只要人源
##   res_h <- cxg_search(c("bone","osteoblast"), organism="Homo sapiens")
##
##   # 打印可下载清单
##   cxg_print(res_h, top=20)
## =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

CXG_API <- "https://api.cellxgene.cziscience.com/curation/v1/datasets"

## 安全取值：字段不存在时给默认（%||% 只在 rlang 里有，这里自己实现）
get_or <- function(x, default) {
  if (is.null(x)) return(default)
  if (length(x) == 0) return(rep(default, max(1, length(x))))
  x[is.na(x)] <- default
  return(x)
}

## ---- 拉取元数据（会缓存到 data/processed，避免重复请求） ----
cxg_fetch <- function(force = FALSE) {
  cache <- "data/processed/cxg_datasets.rds"
  if (!force && file.exists(cache)) {
    d <- readRDS(cache)
    message("使用缓存（", nrow(d), " 个数据集）。要刷新用 cxg_fetch(force=TRUE)")
    return(d)
  }

  message("拉取 CELLxGENE 数据集清单（约 2238 个，需 10-30 秒）…")
  txt <- tryCatch(
    readLines(CXG_API, warn = FALSE),
    error = function(e) {
      stop("无法访问 CELLxGENE。检查网络后重试。\n  原始错误: ", conditionMessage(e))
    }
  )
  lst <- fromJSON(paste(txt, collapse = ""), flatten = TRUE)

  ## 整理成扁平 data.frame
  pick <- function(x) {
    if (is.data.frame(x) && nrow(x) > 0) paste(x$label, collapse = "; ")
    else if (is.data.frame(x)) ""
    else ""
  }

  n <- max(lengths(lst[c("dataset_id", "title")]), 0L)
  d <- data.frame(
    dataset_id  = get_or(lst$dataset_id, NA_character_),
    title       = get_or(lst$title, NA_character_),
    cell_count  = as.integer(get_or(lst$cell_count, NA_integer_)),
    stringsAsFactors = FALSE
  )

  d$organism <- sapply(lst$organism, pick)
  d$tissue   <- sapply(lst$tissues, pick)
  d$assay    <- sapply(lst$assay, pick)
  d$disease  <- sapply(lst$diseases, pick)

  ## 下载链接
  d$h5ad_url <- sapply(lst$assets, function(a) {
    if (is.data.frame(a) && nrow(a) > 0) {
      u <- a$url[a$filetype == "H5AD"]
      if (length(u) > 0) u[1] else NA_character_
    } else NA_character_
  })

  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(d, cache)
  message("✓ 已缓存 ", nrow(d), " 个数据集 → ", cache)
  return(d)
}

## ---- 搜索 ----
cxg_search <- function(keywords, organism = NULL, min_cells = 0) {
  d <- cxg_fetch()

  ## 拼一个可搜索的文本
  blob <- tolower(paste(d$title, d$tissue, d$disease, d$organism, d$assay))

  hit <- rep(FALSE, nrow(d))
  for (k in keywords) {
    hit <- hit | grepl(tolower(k), blob, fixed = TRUE)
  }
  d <- d[hit, , drop = FALSE]

  if (!is.null(organism)) {
    d <- d[grepl(tolower(organism), tolower(d$organism), fixed = TRUE), , drop = FALSE]
  }
  if (min_cells > 0) {
    d <- d[!is.na(d$cell_count) & d$cell_count >= min_cells, , drop = FALSE]
  }

  d <- d[order(-d$cell_count), , drop = FALSE]
  return(d)
}

## ---- 打印 ----
cxg_print <- function(d, top = 15) {
  if (nrow(d) == 0) {
    cat("没有匹配的数据集\n")
    return(invisible(NULL))
  }
  cat(sprintf("共 %d 个匹配，显示前 %d 个：\n\n", nrow(d), min(top, nrow(d))))
  for (i in seq_len(min(top, nrow(d)))) {
    r <- d[i, ]
    cat(sprintf("[%d] %s\n", i, substr(r$title, 1, 68)))
    cat(sprintf("    物种 %-16s 组织 %-26s 细胞 %s\n",
                substr(r$organism, 1, 16), substr(r$tissue, 1, 26),
                format(r$cell_count, big.mark = ",")))
    cat(sprintf("    id   %s\n", substr(r$dataset_id, 1, 36)))
    cat(sprintf("    h5ad %s\n", substr(r$h5ad_url, 1, 92)))
    cat("\n")
  }
}

## ---- 直接跑时的引导 ----
if (sys.nframe() == 0L) {
  cat("==================================================\n")
  cat(" CELLxGENE 数据集搜索\n")
  cat("==================================================\n\n")
  cat("已内置示例，直接跑会搜索「骨相关 + 人源」：\n\n")

  ok <- tryCatch({
    res <- cxg_search(c("bone", "osteoblast", "mesenchymal", "skeletal"),
                      organism = "Homo sapiens")
    cxg_print(res, top = 12)
    TRUE
  }, error = function(e) {
    cat("搜索失败:", conditionMessage(e), "\n")
    FALSE
  })

  cat("\n--------------------------------------------------\n")
  cat("下载方式（Git Bash / Anaconda Prompt）：\n")
  cat("  bash scripts/02_download_matrix.sh <dataset_id> <输出名>\n\n")
  cat("例：\n")
  cat("  bash scripts/02_download_matrix.sh \\\n")
  cat("       df38503b-8dba-4800-957d-6613747924bb bone_dev\n")
  cat("==================================================\n")
}
