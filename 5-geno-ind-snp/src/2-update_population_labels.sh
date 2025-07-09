#!/bin/bash

# 更新EIGENSTRAT .ind文件中的群体标签
# 用法: ./update_population_labels.sh <ind_file> <sample_info_file>

set -e

if [ $# -lt 1 ]; then
    echo "用法: $0 <ind_file> [sample_info_file]"
    echo ""
    echo "参数说明:"
    echo "  ind_file         - EIGENSTRAT格式的.ind文件路径"
    echo "  sample_info_file - 可选，包含样本-群体对应关系的文件"
    echo ""
    echo "如果不提供sample_info_file，将进行交互式标签替换"
    exit 1
fi

IND_FILE="$1"
SAMPLE_INFO_FILE="$2"

# 检查输入文件
if [ ! -f "$IND_FILE" ]; then
    echo "错误: .ind文件不存在: $IND_FILE"
    exit 1
fi

# 备份原始文件
backup_file="${IND_FILE}.backup.$(date +%Y%m%d_%H%M%S)"
cp "$IND_FILE" "$backup_file"
echo "已备份原始文件到: $backup_file"

if [ -n "$SAMPLE_INFO_FILE" ] && [ -f "$SAMPLE_INFO_FILE" ]; then
    echo "使用样本信息文件进行批量更新: $SAMPLE_INFO_FILE"
    
    # 假设sample_info_file格式为: SampleID Population
    # 读取样本信息并更新
    while IFS=$'\t' read -r sample_id population; do
        if [ -n "$sample_id" ] && [ -n "$population" ]; then
            sed -i "s/^$sample_id\t\([MF]\)\t?/$sample_id\t\1\t$population/" "$IND_FILE"
        fi
    done < "$SAMPLE_INFO_FILE"
    
    echo "批量更新完成"
else
    echo "交互式群体标签更新"
    echo "当前.ind文件中的唯一标签:"
    awk '{print $3}' "$IND_FILE" | sort | uniq -c
    echo ""
    
    # 获取所有唯一的群体标签
    unique_labels=$(awk '{print $3}' "$IND_FILE" | sort | uniq)
    
    for label in $unique_labels; do
        if [ "$label" = "?" ]; then
            echo "发现占位符'?', 需要替换"
            echo -n "请输入新的群体名称 (或按Enter跳过): "
            read new_population
            
            if [ -n "$new_population" ]; then
                sed -i "s/\t?$/\t$new_population/" "$IND_FILE"
                echo "已将'?'替换为'$new_population'"
            fi
        else
            echo "发现群体标签: $label"
            echo -n "是否需要替换? (输入新名称或按Enter保持不变): "
            read new_population
            
            if [ -n "$new_population" ]; then
                sed -i "s/\t$label$/\t$new_population/g" "$IND_FILE"
                echo "已将'$label'替换为'$new_population'"
            fi
        fi
    done
fi

echo ""
echo "更新后的群体标签统计:"
awk '{print $3}' "$IND_FILE" | sort | uniq -c

echo ""
echo "群体标签更新完成!"
echo "修改后的文件: $IND_FILE"
echo "备份文件: $backup_file"
