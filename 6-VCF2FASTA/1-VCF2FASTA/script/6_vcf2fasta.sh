#!/bin/bash

# 设置脚本在遇到错误时立即退出，并在使用未定义变量时退出
set -euo pipefail

#########################
# 配置部分
#########################

# 定义参考基因组路径
REFERENCE="/mnt/d/幽门螺旋杆菌/参考序列/NC_000915.fasta"

# 定义 VCF 文件路径（压缩格式）
VCF="/mnt/d/幽门螺旋杆菌/Script/分析结果/merged.vcf.gz"

# 定义输出目录，用于存放生成的 FASTA 序列
OUTPUT_DIR="/mnt/f/Fasta_sequences"
mkdir -p "$OUTPUT_DIR"

# 定义最终的 Multi-FASTA 输出文件
MULTI_FASTA="/mnt/f/all_samples.fasta"

# 定义临时目录，用于存放中间文件
TMP_DIR='/mnt/f/Temp'
mkdir -p "$TMP_DIR"
# 定义过滤后的 VCF 路径
FILTERED_VCF="${TMP_DIR}/filtered.vcf.gz"
bcftools index $FILTERED_VCF
# 函数：在脚本退出时清理临时目录
# cleanup() {
#     rm -rf "$TMP_DIR"
# }
# trap cleanup EXIT

#########################
# 过滤 VCF 文件
#########################

echo "开始过滤 VCF 文件，移除插入变异..."

# 步骤 1：将多等位基因记录拆分为单等位基因记录
bcftools norm -m -any -Oz -o "${TMP_DIR}/biallelic.vcf.gz" "$VCF"

# 步骤 2：索引拆分后的 VCF 文件
bcftools index "${TMP_DIR}/biallelic.vcf.gz"

# 步骤 3：过滤掉插入变异（保留 ALT 长度小于或等于 REF 长度的记录）
bcftools view -i 'strlen(ALT) <= strlen(REF)' "${TMP_DIR}/biallelic.vcf.gz" -Oz -o "$FILTERED_VCF"

# 步骤 4：索引过滤后的 VCF 文件
bcftools index "$FILTERED_VCF"

echo "VCF 文件过滤完成。过滤后的 VCF 路径：$FILTERED_VCF"

#########################
# 生成 Multi-FASTA 文件
#########################

# 提取样本名称并存储到数组中
echo "提取样本名称..."
mapfile -t samples < <(bcftools query -l "$FILTERED_VCF")

# 定义临时 FASTA 目录，用于存放各样本的临时 FASTA 文件
FASTA_TMP_DIR="${TMP_DIR}/fasta_tmp"
mkdir -p "$FASTA_TMP_DIR"

# 导出环境变量以供并行处理函数使用
export REFERENCE FILTERED_VCF OUTPUT_DIR FASTA_TMP_DIR

# 定义处理单个样本的函数
process_sample() {
    local SAMPLE="$1"
    echo "正在处理样本：$SAMPLE"

    # 使用 bcftools 生成单个样本的 FASTA 序列
    if bcftools consensus -f "$REFERENCE" -s "$SAMPLE" "$FILTERED_VCF" > "${FASTA_TMP_DIR}/${SAMPLE}.fasta"; then
        # 修改 FASTA 头部为样本名称，并保存到临时文件
        {
            echo ">${SAMPLE}"
            tail -n +2 "${FASTA_TMP_DIR}/${SAMPLE}.fasta"
        } > "${OUTPUT_DIR}/${SAMPLE}.tmp.fasta"
    else
        echo "处理样本时出错：$SAMPLE" >&2
    fi
}

# 导出处理函数，以供 GNU Parallel 调用
export -f process_sample

# 使用 GNU Parallel 并行处理样本
echo "开始并行生成 FASTA 序列..."
parallel -j 10 --bar process_sample ::: "${samples[@]}"

# 合并所有临时的 FASTA 文件为一个 Multi-FASTA 文件
echo "合并所有样本的 FASTA 序列为一个 Multi-FASTA 文件..."
cat "${OUTPUT_DIR}"/*.tmp.fasta > "$MULTI_FASTA"

# 移除临时的 FASTA 文件
# rm "${OUTPUT_DIR}"/*.tmp.fasta

echo "Multi-FASTA 文件已生成：$MULTI_FASTA"
