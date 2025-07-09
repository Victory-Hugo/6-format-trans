#!/usr/bin/env bash
set -euo pipefail

VCF="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544.vcf.NoN.maf99.WGS.recode.fillN_SNP.vcf.gz"
SAMPLE_LIST="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/conf/7544个样本list.txt"
REF="/mnt/d/幽门螺旋杆菌/参考序列/NC_000915.fasta"
OUTDIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta"
mkdir -p "$OUTDIR"

# # 索引
# bcftools index "$VCF"
# samtools faidx "$REF"

# 过滤掉空行，然后为每个样本生成 fasta
grep -v '^[[:space:]]*$' "$SAMPLE_LIST" | while read -r SAMPLE; do
  echo ">>> processing $SAMPLE"
  bcftools consensus \
    -f "$REF" \
    -s "$SAMPLE" \
    -o "$OUTDIR/${SAMPLE}.fasta" \
    "$VCF"
done

echo "► 所有样本的 fasta 已写入 $OUTDIR/"
