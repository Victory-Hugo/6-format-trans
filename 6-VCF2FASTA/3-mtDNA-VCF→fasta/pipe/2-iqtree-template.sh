#!/usr/bin/bash
# ====== IQ-TREE3 混合模型与定根：极简命令清单 ======
BASE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/6-VCF2FASTA/3-mtDNA-VCF→fasta" #* 基础目录
ALIGN="$BASE_DIR/output/alignment.fasta" #! 或 concat_genes.aln.fasta
PART="$BASE_DIR/output/partition.nex"     #* 分区文件（可选）
OUTGROUP="RSRS"  #* 外群物种名（可选）
mkdir -p "$BASE_DIR/output/iqtree3"
iqtree3 \
  -s "$ALIGN" \
  -p "$PART" \
  -m MF+MERGE \
  -T 16 \
  --prefix "$BASE_DIR/output/iqtree3/分区-模型选择"

# iqtree3 \
#   -s "$ALIGN" \
#   -p "$PART" \
#   -m MFP \
#   -B 1000 --alrt 1000 \
#   -T AUTO \
#   --prefix /output/iqtree3/分区-模型选择-构树

# #*----------蛋白质模型----------
# iqtree3 -s "$ALIGN" \
#        -redo \
#        -m LG+R4 \
#        -nt 64 \
#        -o ${OUTGROUP} \
#        -pre /output/concat_genes_tree

# # ---------- 混合模型 ----------
# # 手动两模型选择
# iqtree3 -s "$ALIGN" \
#         -m 'MIX{JC,HKY}' \
#         -B 1000 -alrt 1000 -T AUTO

# # 加速率异质性（2模型×4速率=8选择）
# iqtree3 -s "$ALIGN" \
#         -m 'MIX{JC,HKY}+G4' \
#         -B 1000 -alrt 1000 -T AUTO

# # ---------- 混合模型优化 ----------
# # 自动选择模型选择
# iqtree3 -s "$ALIGN" \
#         -m MIX+MF -T AUTO

# # 优化并建树
# iqtree3 -s "$ALIGN" \
#         -m MIX+MFP -B 1000 -alrt 1000 -T AUTO

# # AIC准则 / LRT(p=0.05)
# iqtree3 -s "$ALIGN" \
#         -m MIX+MF -merit AIC -T AUTO
# iqtree3 -s "$ALIGN" \
#         -m MIX+MF -lrt 0.05 -T AUTO

# # 限制模型与速率异质性选项
# iqtree3 -s "$ALIGN" \
#         -m MIX+MF -mset HKY,GTR \
#         -mrate E,I,G,I+G -qmax 5 \
#         -mrate-twice 1 -T AUTO

# # ---------- 分区 + 混合 ----------
# iqtree3 -s "$ALIGN" \
#         -p "$PART" \
#         -m MIX+MFP \
#         -B 1000 -alrt 1000 -T AUTO

# # ===================================================
# #                 系统发育树定根
# # ===================================================

# # 外群定根（推荐）
# iqtree3 -s "$ALIGN" \
#         -p "$PART" \
#         -m MFP \
#         -B 1000 -alrt 1000 \
#         -o "$OUTGROUP" -T AUTO

# # 无外群（非可逆模型定根）
# iqtree3 -s "$ALIGN" \
#         -p "$PART" \
#         -B 1000 -T AUTO \
#         --prefix rev_stage

# iqtree3 -s "$ALIGN" \
#         -p rev_stage.best_scheme.nex \
#         --model-joint UNREST \
#         -B 1000 -T AUTO \
#         --prefix nonrev_stage

