## =============================================================================
## 路线C 下载器：从 CELLxGENE 下载【已处理好的】单细胞表达矩阵
##
## 为什么需要这个：
##   Windows 上没有任何可用的 fastq→矩阵 定量工具
##   （CellRanger 仅 Linux；kb-python 无 Windows 二进制；kallisto 无 win-64 conda 包）
##   所以改走 CELLxGENE —— 它直接提供处理好的 H5AD，无需定量。
##
## 实测数据（2026-10-08）：
##   连通性 OK，HTTP 206 支持断点续传，HDF5 魔数 8948 4446 正确
##   小文件 162 MB → 7.2 MB/s；615 MB → 15.7 MB/s（峰值 19 MB/s）
##   9 GB 的大文件会被限速到约 2 MB/s（约 1 小时），夜里挂着下即可
##   若 curl 报连接失败，多半是代理问题，在命令里加 --noproxy '*' 再试
##
## 用法（Git Bash / Anaconda Prompt）：
##   bash 02_download_matrix.sh <数据集ID> <输出名>
##   bash 02_download_matrix.sh df38503b-8dba-4800-957d-6613747924bb bone_dev
## =============================================================================

set -u

OUTDIR="${VK_DATA_DIR:-C:/Users/Fang/Desktop/S3/ciwujia_wuweizi_oa/data/raw}"
mkdir -p "$OUTDIR"

DATASET_ID="${1:-}"
NAME="${2:-dataset}"

if [ -z "$DATASET_ID" ]; then
  echo "用法: bash 02_download_matrix.sh <dataset_id> <输出名>"
  echo ""
  echo "人源骨相关数据集（CELLxGENE，2026-10-08 实测可下载，按推荐度排序）："
  echo "  b6037d0f-c415-45f7-8d01-0e69407fba9c  人骨髓腔参考 117,308 细胞  615 MB  ← 推荐"
  echo "      组织=bone marrow cavity；含骨关节炎 23,307 vs 正常 92,338"
  echo "      细胞类型 12 类：基质细胞 25,498 / 单核 34,156 / 骨细胞 2,812 / 造血前体 22,951"
  echo "      obs 46 列，另有 author_celltype、leiden、predicted.celltype.l1/l2"
  echo ""
  echo "  c54035d5-a021-4bd0-9496-e9ec0cf9dcc8  人间充质细胞   8,030 细胞  317 MB"
  echo " 98ee24fd-6102-4c51-843a-8cf40e6df724  人骨髓类器官  31,040 细胞  359 MB"
  echo "  7b278368-e146-4811-bc55-b2b022586de3  人胎儿骨髓 CITE-seq  8,978 细胞  162 MB"
  echo "  df38503b-8dba-4800-957d-6613747924bb  人骨发育图谱  908,046 细胞  9,000 MB"
  echo ""
  echo "示例："
  echo "  bash scripts/02_download_matrix.sh \\"
  echo "       b6037d0f-c415-45f7-8d01-0e69407fba9c bone_marrow_ref"
  echo ""
  echo "列更多数据集：source('scripts/03_search_cellxgene.R') ; cxg_search('bone')"
  exit 1
fi

URL="https://datasets.cellxgene.cziscience.com/${DATASET_ID}.h5ad"
OUT="$OUTDIR/${NAME}.h5ad"

echo "数据集 : $DATASET_ID"
echo "输出到 : $OUT"
echo ""

## 断点续传 + 重试
attempt=1
while [ "$attempt" -le 8 ]; do
  before=0
  [ -f "$OUT" ] && before=$(stat -c%s "$OUT" 2>/dev/null || echo 0)
  echo "[$(date '+%H:%M:%S')] 第 ${attempt} 次 (已有 $(( before / 1048576 )) MB)"

  curl -L -C - --retry 5 --retry-delay 5 --retry-all-errors \
       --max-time 10800 --speed-time 180 --speed-limit 1024 \
       -o "$OUT" "$URL"

  rc=$?
  after=0
  [ -f "$OUT" ] && after=$(stat -c%s "$OUT" 2>/dev/null || echo 0)

  ## HDF5 魔数校验：89 48 44 46 = .HDF
  magic=$(xxd -p -l 4 "$OUT" 2>/dev/null)

  if [ "$magic" = "89484446" ]; then
    echo "[$(date '+%H:%M:%S')] ✓ 完成，$(( after / 1048576 )) MB，HDF5 格式校验通过"
    break
  elif [ "$after" -gt "$before" ]; then
    echo "[$(date '+%H:%M:%S')] … 续传至 $(( after / 1048576 )) MB，继续"
  else
    echo "[$(date '+%H:%M:%S')] ! 无进展 (curl rc=$rc)，等待重试"
    sleep 12
  fi
  attempt=$(( attempt + 1 ))
done

magic=$(xxd -p -l 4 "$OUT" 2>/dev/null)
if [ "$magic" != "89484446" ]; then
  echo "✗ 下载未完成或文件损坏（魔数=$magic）"
  echo "  重跑本命令即可续传"
  exit 1
fi

echo ""
echo "=========================================================="
echo " 下载完成"
echo "=========================================================="
echo "文件: $OUT ($(( $(stat -c%s "$OUT") / 1048576 )) MB)"
echo ""
echo "下一步（R 里）："
echo "  source('scripts/04_load_h5ad.R')   # h5ad → Seurat 对象"
echo "=========================================================="
