#!/bin/bash

# 退出条件：任何命令失败时立即退出
set -euo pipefail

# ================== 配置参数 ==================
VCF_PATH="/mnt/f/0_现代DNA处理流程/output/merge/merged_filtered_去除hots.vcf.gz"
OUT_PREFIX="/mnt/f/0_现代DNA处理流程/output/merge_fasta/merged_filtered_anno_去除hots"
PYTHON_PATH="/home/luolintao/miniconda3/bin/python3"
SCRIPT_PATH="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-VCF2FASTA/1-VCF2FASTA/script/6_vcf2fasta_通过plink.py"
MAPPING_FILE="/home/luolintao/07_20K_CPGDP/1_单倍群分型/data/质量控制_ID_Hap.tsv" #? 可以不添加这个参数，注释rename_fasta_ids函数
THREADS=32


# ================== 函数定义 ==================

function run_plink() {
    echo "[INFO] 正在使用 PLINK 处理 VCF 文件..."
    /mnt/e/Scientifc_software/plink_linux_x86_64_20241022/plink --vcf "$VCF_PATH" \
        --threads "$THREADS" \
        --recode \
        --double-id \
        --out "$OUT_PREFIX"
    echo "[INFO] PLINK 处理完成！"
}

function ped_to_fasta() {
    local ped_file="${OUT_PREFIX}.ped"
    local fasta_file="${OUT_PREFIX}.fasta"

    echo "[INFO] 正在将 PED 文件转为 FASTA 格式..."
    awk '{
        printf(">%s\n", $1);
        seq = "";
        for (i = 7; i <= NF; i += 2) {
            seq = seq $i;
        }
        printf("%s\n", seq);
    }' "$ped_file" > "$fasta_file"
    echo "[INFO] PED 转换为 FASTA 成功！"
}

function rename_fasta_ids() {
    local input_fasta="${OUT_PREFIX}.fasta"
    local output_fasta="${OUT_PREFIX}_rename.fasta"

    echo "[INFO] 正在调用 Python 脚本重命名 FASTA 序列..."
    "$PYTHON_PATH" "$SCRIPT_PATH" \
        --mapping "$MAPPING_FILE" \
        --input_fasta "$input_fasta" \
        --output_fasta "$output_fasta"
    echo "[INFO] FASTA 序列重命名完成！"
}

# ================== 主流程执行 ==================

run_plink
ped_to_fasta
# rename_fasta_ids
rm "$OUT_PREFIX".ped
rm "$OUT_PREFIX".map
rm "$OUT_PREFIX".log
rm "$OUT_PREFIX".nosex
echo "[DONE] 所有步骤执行完毕。"
