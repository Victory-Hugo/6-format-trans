#!/usr/bin/env bash
set -euo pipefail

# ======== 可调参数 ========
BASE_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta"
ALIGN_DIR="${BASE_DIR}/Core_CDS_Align"
OUT_DIR="${BASE_DIR}/Core_CDS_strain"
FINAL_OUT="${BASE_DIR}/Core_CDS_concat.aln.fasta"
PY_DIR='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-VCF2FASTA/2-FASTA_CDS/src/1-3-merge_optimized.py'
PY_EXE='/home/luolintao/miniconda3/envs/pyg/bin/python3'

# ======== 导出环境变量，供 Python 读取 ========
export ALIGN_DIR OUT_DIR

echo "[INFO] 开始执行优化版本..."

# ======== 执行优化的 Python 脚本（不再需要并行，因为已经在内部优化） ========
${PY_EXE} ${PY_DIR}

# ======== 输出完成信息 ========
echo "[DONE] 第一步执行完成."

# ======== 开始第二步：更高效的文件合并 ========
echo "[INFO] 开始合并所有输出文件..."

# 检查输出目录是否存在文件
if [ ! "$(ls -A ${OUT_DIR}/*.fasta 2>/dev/null)" ]; then
    echo "[ERROR] 在 ${OUT_DIR} 中未找到 .fasta 文件"
    exit 1
fi

# 使用更高效的方式合并文件
# 方法1：直接使用cat一次性合并所有文件（推荐）
cat ${OUT_DIR}/*.fasta > "${FINAL_OUT}"

echo "[DONE] 文件合并完成，输出文件: ${FINAL_OUT}"

# 显示统计信息
echo "[INFO] 统计信息:"
echo "  - 输入文件数: $(ls ${OUT_DIR}/*.fasta | wc -l)"
echo "  - 输出文件大小: $(du -h "${FINAL_OUT}" | cut -f1)"
echo "  - 序列数量: $(grep -c '^>' "${FINAL_OUT}")"

echo "[DONE] 全部执行完成."
