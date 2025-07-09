#############################
##第三步，提取每个基因的核心区域#
#############################
#!/bin/bash

# 设置输入和输出目录
INPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/CDS_Align"
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/Core_CDS_Align"
SCRIPT_PATH="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-VCF2FASTA/2-FASTA_CDS/src/1-2_GetCoreCDS.py"

# 设置并行任务数（根据你的CPU核心数调整）
NUM_JOBS=4

# 检查 GNU Parallel 是否已安装
if ! command -v parallel &> /dev/null; then
    echo "错误: GNU Parallel 未安装。请安装后重试。"
    exit 1
fi

# 检查 Python 脚本是否存在
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "错误: Python 脚本 $SCRIPT_PATH 不存在。"
    exit 1
fi

# 创建输出目录（如果不存在）
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "输出目录不存在，正在创建: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
    if [ $? -ne 0 ]; then
        echo "错误: 无法创建输出目录 $OUTPUT_DIR。"
        exit 1
    fi
fi

# 导出变量以供 parallel 使用
export SCRIPT_PATH
export OUTPUT_DIR

# 定义处理单个文件的函数
process_file() {
    INPUT_FILE="$1"
    FILENAME=$(basename "$INPUT_FILE")
    
    # 修改输出文件名：将 '_merged' 替换为 '_core'
    NEW_FILENAME="${FILENAME/_merged/_core}"
    OUTPUT_FILE="$OUTPUT_DIR/$NEW_FILENAME"

    echo "处理文件: $FILENAME -> $NEW_FILENAME"

    /home/luolintao/miniconda3/envs/pyg/bin/python3 "$SCRIPT_PATH" -i "$INPUT_FILE" -o "$OUTPUT_FILE"

    if [ $? -eq 0 ]; then
        echo "成功生成核心对齐文件: $OUTPUT_FILE"
    else
        echo "错误: 处理文件 $FILENAME 失败。"
    fi

    echo "----------------------------------------"
}

export -f process_file

# 使用 GNU Parallel 并行处理所有 .aln 文件
find "$INPUT_DIR" -type f -name "*.aln" | parallel -j "$NUM_JOBS" process_file {}

echo "所有文件处理完成。"

