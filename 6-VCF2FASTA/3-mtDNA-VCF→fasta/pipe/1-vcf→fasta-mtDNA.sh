#!/usr/bin/env bash

set -euo pipefail

INPUT_VCF="/home/luolintao/0_现代DNA处理流程/output/Archive/merged_clean.vcf.gz"

BASE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/6-VCF2FASTA/3-mtDNA-VCF→fasta"

OUTPUT_DIR="${BASE_DIR}/output"
REF_FASTA="${BASE_DIR}/data/rcrs.fasta"
SAMPLE_LIST="${BASE_DIR}/data/ID.txt"
VCF2PHYLO="${BASE_DIR}/python/vcf2phylo.py"

MAX_NEXUS_TAXA=20000

PYTHON3="/home/luolintao/miniconda3/envs/BigLin/bin/python3"

mkdir -p "${OUTPUT_DIR}"

"${PYTHON3}" "${VCF2PHYLO}" \
    -v "${INPUT_VCF}" \
    -o "${OUTPUT_DIR}" \
    -r "${REF_FASTA}" \
    -S "${SAMPLE_LIST}" \
    --max-nexus-taxa "${MAX_NEXUS_TAXA}"
#? 后续命令，使用iqtree3进行系统发育分析
iqtree3 \
  -s /mnt/c/Users/Administrator/Desktop/output/alignment.fasta \
  -p /mnt/c/Users/Administrator/Desktop/output/partition.nex \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T AUTO \
  --prefix /mnt/c/Users/Administrator/Desktop/output/alignment_MFP