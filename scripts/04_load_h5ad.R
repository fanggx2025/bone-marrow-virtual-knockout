## =============================================================================
## 04_load_h5ad.R —— 把 CELLxGENE 的 H5AD 读成 Seurat 对象（纯 rhdf5，无需 Python）
##
## 为什么不用 zellkonverter：
##   zellkonverter 依赖 basilisk → 需要完整 Python 环境。本机 pip 装 pysam 失败
##   （源码包需 C 编译），conda 也无 win-64 的 kb-python。本机已有 rhdf5 2.46.1，
##   可直接读 HDF5。
##
## 本文件已用真实 CELLxGENE h5ad 实测通过（dataset 7b278368，8,978 细胞）。
## 实测踩到的 4 个坑（代码里都已处理）：
##   坑1  /X 里【没有 shape】，只有 data / indices / indptr
##        → 行数 = length(indptr)-1；列数必须取 nrow(var)，
##          不能取 max(indices)+1（实测 35468 vs 35469，末尾基因全为 0）
##   坑2  分类列（cell_type / tissue / disease / sex …）是 H5I_GROUP，
##        里面是 categories(STRING) + codes(INTEGER, 0-based, -1=缺失)
##        → 必须用 categories[codes + 1] 还原
##   坑3  /var/_index 是 Ensembl ID，SYMBOL 在 /var/gene_symbols（也是 categorical）
##        → 默认用 SYMBOL 做行名，Ensembl ID 另存一份备用
##   坑4  ENUM 类型（is_primary_data / highly_variable）rhdf5 读出来是 factor
##        → 统一转 character
##
## 用法：
##   source("scripts/04_load_h5ad.R")
##   seu <- load_h5ad("data/raw/bone.h5ad", "Bone")
## =============================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(Seurat)
})

## ---------------------------------------------------------------------------
## 工具 1：读「普通 dataset」（含 ENUM → character）
## ---------------------------------------------------------------------------
h5_read_dataset <- function(f, path) {
  v <- tryCatch(h5read(f, path), error = function(e) NULL)
  if (is.null(v)) return(NULL)
  if (is.factor(v)) v <- as.character(v)
  ## 二维以上的（obsm 的 PCA 等）跳过
  if (!is.null(dim(v)) && length(dim(v)) > 1L) return(NULL)
  v <- as.vector(v)
  if (length(v) == 0L) return(NULL)
  return(v)
}

## ---------------------------------------------------------------------------
## 工具 2：读「分类变量」 group = {categories, codes}
## ---------------------------------------------------------------------------
h5_read_categorical <- function(f, path) {
  items <- h5ls(f, datasetinfo = FALSE)
  nm <- items$name[items$group == path]
  if (!all(c("categories", "codes") %in% nm)) return(NULL)

  cats <- tryCatch(as.character(h5read(f, paste0(path, "/categories"))),
                   error = function(e) NULL)
  cods <- tryCatch(as.integer(h5read(f, paste0(path, "/codes"))),
                   error = function(e) NULL)
  if (is.null(cats) || is.null(cods)) return(NULL)

  out <- rep(NA_character_, length(cods))
  ok  <- !is.na(cods) & cods >= 0L & cods < length(cats)
  out[ok] <- cats[cods[ok] + 1L]        # Python 0-based → R 1-based
  return(out)
}

## ---------------------------------------------------------------------------
## 工具 3：取基因 SYMBOL（按可靠性从高到低尝试）
##
## 坑13（2026-10-08 实测，b6037d0f / 117,308 细胞）：
##   CELLxGENE 的 /var/gene_symbols 常常是「压缩 categorical」——
##   categories 只有 1,269 个，而基因有 25,574 个，
##   也就是只有约 5% 的基因有 SYMBOL，其余 code = -1，还原出来大片 NA。
##   表现：前 20 个基因全是 NA，CD3E / PTPRC / CD68 / SOX9 / SP7 全部"缺失"
##         （这些基因绝不可能没有，一看到就该警觉）。
##   正解：完整的基因名在 /raw/var/gene_symbols —— 纯 STRING，25,574 个，0 个 NA。
##   判据：不要只看"能不能读出来"，必须校验 length == 基因数 且 sum(is.na) == 0。
## ---------------------------------------------------------------------------
h5_gene_symbols <- function(f, n_expected, verbose = TRUE) {
  nm_all <- h5ls(f, datasetinfo = FALSE)
  is_ds  <- function(g, n) any(nm_all$group == g & nm_all$name == n &
                               nm_all$otype == "H5I_DATASET")
  good <- function(v) !is.null(v) && length(v) == n_expected && sum(is.na(v)) == 0L

  ## ① /raw/var/gene_symbols —— 完整 STRING，最可靠
  if (is_ds("/raw/var", "gene_symbols")) {
    v <- tryCatch(as.character(h5read(f, "/raw/var/gene_symbols")),
                  error = function(e) NULL)
    if (good(v)) {
      if (verbose) message("        基因名来源: /raw/var/gene_symbols (完整STRING)")
      return(v)
    }
  }
  ## ② /var/gene_symbols 本身是纯 STRING
  if (is_ds("/var", "gene_symbols")) {
    v <- tryCatch(as.character(h5read(f, "/var/gene_symbols")),
                  error = function(e) NULL)
    if (good(v)) {
      if (verbose) message("        基因名来源: /var/gene_symbols (STRING)")
      return(v)
    }
  }
  ## ③ /var/gene_symbols 是 categorical —— 先校验 NA 比例
  if (any(nm_all$group == "/var/gene_symbols" & nm_all$name == "categories")) {
    v <- tryCatch(h5_read_categorical(f, "/var/gene_symbols"),
                  error = function(e) NULL)
    if (good(v)) {
      if (verbose) message("        基因名来源: /var/gene_symbols (categorical)")
      return(v)
    }
    if (!is.null(v) && verbose)
      message("        ⚠ /var/gene_symbols 仅 ", sum(!is.na(v)), "/", n_expected,
              " 个非NA，弃用（坑13：压缩categorical）")
  }
  return(NULL)
}

## ---------------------------------------------------------------------------
## 读 AnnData 的 obs / var（自动判断 dataset 还是 categorical）
## ---------------------------------------------------------------------------
read_anndata_df <- function(f, group, verbose = TRUE) {
  items <- h5ls(f, datasetinfo = FALSE)
  sub <- items[items$group == group, ]
  if (nrow(sub) == 0L) return(NULL)

  out <- list()
  for (i in seq_len(nrow(sub))) {
    nm   <- sub$name[i]
    path <- paste0(group, "/", nm)
    v <- if (sub$otype[i] == "H5I_GROUP") {
      h5_read_categorical(f, path)
    } else {
      h5_read_dataset(f, path)
    }
    if (is.null(v)) next
    out[[nm]] <- v
  }

  if (length(out) == 0L) return(NULL)

  ## 行名
  idx <- NULL
  if ("_index" %in% names(out)) {
    idx <- as.character(out[["_index"]])
    out[["_index"]] <- NULL
  }

  ## 只保留长度等于众数的列（对齐主表）
  lens <- vapply(out, length, integer(1))
  if (length(lens) == 0L) return(NULL)
  target <- as.integer(names(sort(table(lens), decreasing = TRUE))[1])
  out <- out[lens == target]
  if (length(out) == 0L) return(NULL)

  df <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(idx) && nrow(df) == length(idx)) rownames(df) <- idx

  return(df)
}

## ---------------------------------------------------------------------------
## 读 /X（稀疏 CSR → dgRMatrix）
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## 工具 4：CSR 健全性检查（坑14）
##
## 坑14（2026-10-08 实测，b6037d0f / 117,308 细胞）：
##   CELLxGENE 某些 h5ad 的 /X 是「等间距 padding」的伪 CSR ——
##   diff(indptr) 恒等于某个常数（本例 1,269），标准差 = 0，
##   与 obs/n_genes（真实检出基因数）完全不相关。
##   而 /raw/X 的 diff(indptr) 与 obs/n_genes 相关系数高达 0.9986。
##   后果：用 /X 建出来的矩阵每细胞"检出基因数"全部相同，
##         再做 min.features 过滤时行为诡异（本例 117,308 → 只剩 9,801）。
##   判据：cv = sd(diff(indptr)) / mean(diff(indptr))，cv < 0.01 即判为伪 CSR。
##   注意：只看魔数/维度/能不能读出来，全都发现不了，必须看 indptr 的变异。
## ---------------------------------------------------------------------------
h5_csr_sanity <- function(f, group = "/X") {
  ip <- tryCatch(as.numeric(h5read(f, paste0(group, "/indptr"))),
                 error = function(e) NULL)
  if (is.null(ip) || length(ip) < 3L) return(NULL)
  d <- diff(ip)
  if (length(d) == 0L || !is.finite(mean(d)) || mean(d) == 0) return(NULL)
  cv <- sd(d) / mean(d)
  list(n = length(d), mean = mean(d), sd = sd(d), cv = cv, uniform = cv < 0.01)
}

read_anndata_X <- function(f, ncol_hint = NULL, group = "/X", verbose = TRUE) {
  items <- h5ls(f, datasetinfo = FALSE)
  xn <- items$name[items$group == group]

  if (all(c("data", "indices", "indptr") %in% xn)) {
    if (verbose) message("        ", group, " = 稀疏 CSR（无 shape，维度由 indptr + var 推断）")

    data    <- h5read(f, paste0(group, "/data"))
    indices <- as.integer(h5read(f, paste0(group, "/indices")))
    indptr  <- as.numeric(h5read(f, paste0(group, "/indptr")))  # 可能超 2^31，用 double

    nr <- length(indptr) - 1L
    nc <- if (!is.null(ncol_hint)) as.integer(ncol_hint) else (max(indices) + 1L)

    if (verbose) {
      message("        非零元素: ", format(length(data), big.mark = ","),
              "   行=细胞: ", nr, "   列=基因: ", nc)
    }

    i_vec <- rep(seq_len(nr), times = diff(indptr))
    if (length(i_vec) != length(data)) {
      stop("CSR 结构不一致: length(i)=", length(i_vec),
           " 但 length(data)=", length(data))
    }

    m <- sparseMatrix(i = i_vec, j = indices + 1L, x = as.numeric(data),
                      dims = c(nr, nc), repr = "R")   # "R" = 行优先 = CSR，免排序
    rm(i_vec, data, indices, indptr); gc()
    return(m)
  }

  ## 稠密
  if (verbose) message("        ", group, " = 稠密")
  m <- t(h5read(f, group))
  return(as(as.matrix(m), "CsparseMatrix"))
}

## ---------------------------------------------------------------------------
## 主函数
## ---------------------------------------------------------------------------
load_h5ad <- function(h5_path, project = "scDataset",
                      min_cells = 3, min_features = 200,
                      use_raw = FALSE, gene_name = "symbol",
                      verbose = TRUE) {

  if (!file.exists(h5_path)) stop("文件不存在: ", h5_path)
  if (verbose) message("[load_h5ad] ", basename(h5_path),
                       "  (", round(file.size(h5_path) / 1048576, 1), " MB)")

  f <- H5Fopen(h5_path)
  on.exit(try(H5Fclose(f), silent = TRUE), add = TRUE)

  ## use_raw 只影响表达矩阵 X；
  ## var / obs 始终取信息最全的 /var、/obs ——
  ## 因为 /raw/var 只有少数几列，且 /raw/obs 根本不存在（会让元数据全丢）

  ## ---- [1/4] 基因注释 /var -------------------------------------------------
  if (verbose) message("  [1/4] 读基因注释 /var")
  var_df <- read_anndata_df(f, "/var", verbose)
  if (is.null(var_df)) stop("读不到 /var，文件可能损坏或不是 AnnData 格式")

  ## 基因名：优先 SYMBOL（坑13：走 h5_gene_symbols 三级回退）
  n_col <- nrow(var_df)
  gsym  <- if (identical(gene_name, "symbol")) {
    h5_gene_symbols(f, n_expected = n_col, verbose = verbose)
  } else NULL

  genes <- if (!is.null(gsym)) {
    gsym
  } else {
    if (verbose && identical(gene_name, "symbol"))
      message("        ⚠ 取不到 SYMBOL，回退 Ensembl ID（/var/_index）")
    as.character(rownames(var_df))
  }
  ## 去重，避免 Seurat 报重复行名
  if (anyDuplicated(genes)) {
    genes <- make.unique(genes)
    if (verbose) message("        ⚠ 基因名有重复，已 make.unique")
  }
  if (verbose) message("        基因数: ", n_col, "   例: ", paste(head(genes, 5), collapse = ", "))

  ## 基因名自检：这几个基因任何人类 10x 数据都必然存在，
  ## 一旦大量"缺失"就说明 SYMBOL 取错了（坑13 的典型症状）
  if (verbose && !all(grepl("^ENSG", genes))) {
    housekeep <- c("ACTB", "GAPDH", "CD3E", "PTPRC", "CD68", "MS4A1", "LYZ")
    miss <- setdiff(housekeep, genes)
    if (length(miss) > 2) {
      message("        ⚠ 基因名自检未通过，缺失: ", paste(miss, collapse = ", "),
              "\n          强烈建议 gene_name = 'ensembl' 或检查文件")
    } else {
      message("        基因名自检: OK（管家/谱系基因均在）")
    }
  }

  ## ---- [2/4] 细胞注释 /obs -------------------------------------------------
  if (verbose) message("  [2/4] 读细胞注释 /obs")
  obs_df <- read_anndata_df(f, "/obs", verbose)
  if (!is.null(obs_df) && verbose) {
    message("        细胞数: ", nrow(obs_df), "   元数据列: ", ncol(obs_df))
    if ("cell_type" %in% names(obs_df)) {
      message("        cell_type 种类: ", length(unique(obs_df$cell_type)))
    }
  }

  ## ---- [3/4] 表达矩阵 ------------------------------------------------------
  ## 坑14 自动择优：/X 可能是等间距 padding 的伪 CSR，此时改用 /raw/X
  x_group <- if (use_raw) "/raw/X" else "/X"

  if (!use_raw) {
    s <- h5_csr_sanity(f, "/X")
    if (!is.null(s) && s$uniform) {
      sr <- h5_csr_sanity(f, "/raw/X")
      if (!is.null(sr) && !sr$uniform) {
        if (verbose)
          message("        ⚠ /X 的 indptr 完全均匀（每细胞恒为 ", round(s$mean, 1),
                  " 个非零, CV=", round(s$cv, 4), "）→ 判为伪CSR，自动改用 /raw/X（坑14）")
        x_group <- "/raw/X"
      } else if (verbose) {
        message("        ⚠ /X 的 indptr 均匀，但 /raw/X 也不可用，仍用 /X（结果可能失真）")
      }
    } else if (!is.null(s) && verbose) {
      message("        /X 的 indptr CV=", round(s$cv, 3), "（正常非均匀）")
    }
  }
  if (verbose) message("  [3/4] 读表达矩阵 ", x_group)
  m <- read_anndata_X(f, ncol_hint = n_col, group = x_group, verbose = verbose)

  if (!is.null(obs_df) && nrow(obs_df) == nrow(m)) {
    rownames(m) <- rownames(obs_df)
  } else {
    rownames(m) <- paste0("cell_", seq_len(nrow(m)))
  }
  colnames(m) <- genes

  H5Fclose(f)

  ## ---- [4/4] 构建 Seurat ---------------------------------------------------
  if (verbose) message("  [4/4] 构建 Seurat 对象（转置成 基因×细胞）")
  seu <- CreateSeuratObject(
    counts = as(t(m), "CsparseMatrix"),     # dgRMatrix → dgCMatrix
    project = project,
    min.cells = min_cells,
    min.features = min_features
  )
  rm(m); gc()

  ## 挂载元数据
  if (!is.null(obs_df)) {
    common <- intersect(rownames(obs_df), colnames(seu))
    if (length(common) > 0L) {
      md <- obs_df[common, , drop = FALSE]
      md <- md[, vapply(md, function(x) !is.list(x), logical(1)), drop = FALSE]
      ## 已是 character，不再转 factor
      seu <- AddMetaData(seu, md)
      if (verbose) message("        已挂载 ", ncol(md), " 个元数据列")
    }
  }

  if (verbose) {
    message("\n✓ 完成")
    message("  细胞数: ", ncol(seu), "   基因数: ", nrow(seu))
    if ("cell_type" %in% names(seu@meta.data)) {
      message("  细胞类型 top8:")
      tb <- sort(table(seu$cell_type), decreasing = TRUE)
      for (k in seq_len(min(8, length(tb)))) {
        message("    ", sprintf("%-46s %s", names(tb)[k], tb[k]))
      }
    }
  }
  return(seu)
}

## ---------------------------------------------------------------------------
## 只探测结构、不读矩阵（选数据集时用，秒级）
## ---------------------------------------------------------------------------
inspect_h5ad <- function(h5_path) {
  f <- H5Fopen(h5_path)
  on.exit(try(H5Fclose(f), silent = TRUE), add = TRUE)

  items <- h5ls(f, datasetinfo = FALSE)
  cat("文件: ", basename(h5_path), "  ",
      round(file.size(h5_path) / 1048576, 1), " MB\n", sep = "")

  ## 细胞数 / 基因数
  n_cell <- NA; n_gene <- NA
  ip <- tryCatch(h5read(f, "/X/indptr"), error = function(e) NULL)
  if (!is.null(ip)) n_cell <- length(ip) - 1L
  vi <- tryCatch(h5read(f, "/var/_index"), error = function(e) NULL)
  if (!is.null(vi)) n_gene <- length(vi)
  cat("细胞数: ", n_cell, "   基因数: ", n_gene, "\n", sep = "")

  ## obs 列
  obs_nm <- unique(items$name[items$group == "/obs"])
  cat("\nobs 列 (", length(obs_nm), "): ", paste(sort(obs_nm), collapse = ", "), "\n", sep = "")

  ## 分类列取值
  for (k in c("cell_type", "tissue", "disease", "sex", "development_stage", "assay")) {
    v <- if (paste0("/obs/", k) %in% paste(items$group, items$name, sep = "/")) {
      h5_read_categorical(f, paste0("/obs/", k))
    } else {
      tryCatch(as.character(h5read(f, paste0("/obs/", k))), error = function(e) NULL)
    }
    if (is.null(v)) next
    tb <- sort(table(v), decreasing = TRUE)
    cat("\n-- ", k, " (", length(tb), " 类) --\n", sep = "")
    print(head(tb, 12))
  }
  cat("\n顶层: ", paste(h5ls(f, recursive = FALSE)$name, collapse = ", "), "\n", sep = "")
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## 直接跑时的自检
## ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  cat("==================================================\n")
  cat(" h5ad → Seurat 加载器（自检）\n")
  cat("==================================================\n\n")

  args <- commandArgs(trailingOnly = TRUE)
  target <- if (length(args) >= 1) args[1] else
    "C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/data/raw/bone.h5ad"

  if (file.exists(target)) {
    mode <- if (length(args) >= 2) args[2] else "load"
    if (identical(mode, "inspect")) {
      inspect_h5ad(target)
    } else {
      seu <- load_h5ad(target, "scDataset")
      cat("\n元数据列名:\n"); print(names(seu@meta.data))
      cat("\n前 5 个细胞:\n"); print(head(seu@meta.data, 5))
    }
  } else {
    cat("未找到文件:\n  ", target, "\n\n")
    cat("用法:\n")
    cat("  Rscript scripts/04_load_h5ad.R <h5ad路径>            # 载入成 Seurat\n")
    cat("  Rscript scripts/04_load_h5ad.R <h5ad路径> inspect    # 只看结构\n\n")
    cat("在 RStudio 里:\n")
    cat("  source('scripts/04_load_h5ad.R')\n")
    cat("  seu <- load_h5ad('data/raw/bone.h5ad', 'Bone')\n\n")
    cat("先下载:\n")
    cat("  bash scripts/02_download_matrix.sh <dataset_id> bone\n")
  }
  cat("==================================================\n")
}
