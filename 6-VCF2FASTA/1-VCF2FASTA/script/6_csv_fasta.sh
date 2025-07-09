#!/bin/bash

# 定义输入目录
INPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/Final_CSV_files/"

# 定义参考序列和输出目录
REFERENCE_FASTA="/mnt/d/幽门螺旋杆菌/参考序列/NC_000915.fasta"
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/不考虑InDel/"

# 定义日志文件路径，用于记录已处理的文件
LOG_FILE="${OUTPUT_DIR}/processed_files.log"

# 检查输出目录是否存在，不存在则创建
mkdir -p "$OUTPUT_DIR"

# 如果日志文件不存在，创建一个空的日志文件
# 并根据已存在的FASTA文件初始化日志
if [ ! -f "$LOG_FILE" ]; then
    touch "$LOG_FILE"
    echo "创建新的处理日志文件并根据现有FASTA文件进行初始化..."
    
    # 查找已经存在的所有FASTA文件
    find "$OUTPUT_DIR" -type f -name "*.fasta" | while read fasta_file; do
        # 提取文件名（不含路径和扩展名）
        basename=$(basename "$fasta_file" .fasta)
        # 找到对应的CSV文件路径
        csv_file=$(find "$INPUT_DIR" -type f -name "${basename}.csv")
        # 如果找到对应的CSV文件，则将其添加到已处理日志中
        if [ -n "$csv_file" ]; then
            echo "$csv_file" >> "$LOG_FILE"
        fi
    done
    
    echo "日志初始化完成，已找到 $(wc -l < "$LOG_FILE") 个已处理的文件"
else
    echo "使用现有日志文件继续处理: $LOG_FILE"
    echo "已处理文件数: $(wc -l < "$LOG_FILE")"
fi

# 创建一个临时文件列表，包含所有需要处理的CSV文件
TEMP_FILE_LIST="${OUTPUT_DIR}/all_csv_files.txt"
find "$INPUT_DIR" -type f -name "*.csv" > "$TEMP_FILE_LIST"
total_files=$(wc -l < "$TEMP_FILE_LIST")
echo "共发现CSV文件: $total_files"

# 创建一个新的临时文件列表，只包含尚未处理的文件
PENDING_FILE_LIST="${OUTPUT_DIR}/pending_csv_files.txt"
cat "$TEMP_FILE_LIST" | grep -v -f "$LOG_FILE" > "$PENDING_FILE_LIST" || true
pending_files=$(wc -l < "$PENDING_FILE_LIST")
echo "待处理CSV文件: $pending_files"
echo "已跳过处理: $((total_files - pending_files)) 个文件"

if [ "$pending_files" -eq 0 ]; then
    echo "所有文件都已处理完毕，无需进一步处理。"
    rm -f "$TEMP_FILE_LIST" "$PENDING_FILE_LIST"
    exit 0
fi

# 使用GNU Parallel并行处理未处理的CSV文件，并记录已处理的文件
cat "$PENDING_FILE_LIST" | parallel -j 4 --joblog "${OUTPUT_DIR}/parallel_joblog.txt" \
    "echo 处理: {}; \
    csv_basename=\$(basename {} .csv); \
    output_fasta=\"${OUTPUT_DIR}/\${csv_basename}.fasta\"; \
    if [ ! -f \"\$output_fasta\" ]; then \
        /home/luolintao/miniconda3/envs/pyg/bin/python3 \
        /mnt/f/OneDrive/文档（科研）/脚本/Download/3-VCF2FASTA/1-VCF2FASTA/script/6_csv_fasta_noInDel.py \
        {} \"$REFERENCE_FASTA\" \"$OUTPUT_DIR\" && \
        echo {} >> \"$LOG_FILE\" && \
        echo \"完成: \${csv_basename}\"; \
    else \
        echo \"跳过: \${csv_basename} (FASTA已存在)\"; \
        echo {} >> \"$LOG_FILE\"; \
    fi"

# 清理临时文件
rm -f "$TEMP_FILE_LIST" "$PENDING_FILE_LIST"

echo "CSV文件处理已完成。"
processed_files=$(wc -l < "$LOG_FILE")
echo "已处理文件总数: $processed_files / $total_files"
