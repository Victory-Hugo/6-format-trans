#!/bin/bash
set -euo pipefail

# 定义临时目录
temp_dir=$(mktemp -d)

# 拷贝所需ctl文件到临时目录（假设你希望在临时目录操作）
cp /mnt/f/0_现代DNA处理流程/script/2_ML_MRCA/conf/baseml2.ctl \
    "$temp_dir/"

# 切换到临时目录
cd "$temp_dir"

# 运行PAML baseml
/mnt/e/Scientifc_software/paml4.10.8/bin/baseml \
    /mnt/f/0_现代DNA处理流程/script/2_ML_MRCA/conf/baseml2.ctl

# 可选：保存输出结果
# cp mlb /mnt/e/Scientifc_software/paml-master/ML_MRCA/output

# 清理临时目录
cd /
rm -rf "$temp_dir"
