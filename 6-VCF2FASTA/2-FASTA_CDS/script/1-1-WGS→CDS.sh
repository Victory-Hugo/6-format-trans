#!/bin/bash
##################################
##高效并行流式处理：CDS对齐文件####
##跳过中间步骤，减少内存使用######
##使用多进程并行处理加速########
##################################

# 定义变量
BASE_DIR='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-VCF2FASTA/2-FASTA_CDS'
PYTHON3_PATH="/home/luolintao/miniconda3/envs/pyg/bin/python3"
cds_file="${BASE_DIR}/conf/posCDS.csv" # CDS位置文件

# 输入和输出路径
input_folder="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/不考虑InDel/7544_Sample/" # 输入FASTA文件夹
output_folder="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/CDS_Align_new" # 直接输出对齐文件

# 并行进程数 (可根据CPU核心数调整，建议不超过CPU核心数)
processes=8

echo "开始高效并行流式处理CDS提取和对齐..."
echo "输入文件夹: $input_folder"
echo "输出文件夹: $output_folder"
echo "CDS文件: $cds_file"
echo "并行进程数: $processes"
echo "模式: 并行流式处理 (最大化性能，最小化内存使用)"

# 检查输入文件夹是否存在
if [ ! -d "$input_folder" ]; then
    echo "错误: 输入文件夹不存在: $input_folder"
    exit 1
fi

# 检查CDS文件是否存在
if [ ! -f "$cds_file" ]; then
    echo "错误: CDS文件不存在: $cds_file"
    exit 1
fi

# 创建输出目录
mkdir -p "$output_folder"

# 记录开始时间
start_time=$(date +%s)

# 运行并行流式Python脚本
${PYTHON3_PATH} "${BASE_DIR}/src/1-1-WGS→CDS.PY" \
    --input_folder "$input_folder" \
    --cds_file "$cds_file" \
    --output_folder "$output_folder" \
    --processes "$processes"

# 计算运行时间
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo ""
echo "处理完成！"
echo "总运行时间: ${runtime} 秒"
