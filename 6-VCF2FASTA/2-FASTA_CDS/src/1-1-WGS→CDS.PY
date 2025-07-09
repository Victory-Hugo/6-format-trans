#!/usr/bin/env python3
"""
高效流式处理：使用并行处理和更优化的内存管理
每个CDS使用独立进程处理，减少内存峰值
"""

import argparse
import os
from Bio import SeqIO
import glob
import multiprocessing as mp
from functools import partial

# 定义互补碱基字典
complementary_bases = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '-': '-', 'N': 'N'
}

def get_complementary_sequence(sequence):
    """获取反向互补序列"""
    return ''.join(
        complementary_bases.get(base, 'N')
        for base in reversed(sequence)
    )

def load_cds_info(cds_file_path):
    """加载CDS位置信息"""
    cds_info = []
    with open(cds_file_path, 'r') as csv_file:
        for line in csv_file:
            fields = line.strip().split('\t')
            cds_info.append({
                'start': int(fields[3]) - 1,  # 起始位置 (0基准)
                'end': int(fields[4]) - 1,    # 终止位置 (0基准)
                'strand': fields[6],          # 链方向
                'gene_name': fields[2] if len(fields) > 2 else f"CDS_{len(cds_info) + 1}"
            })
    return cds_info

def extract_single_cds_from_fasta(fasta_file_path, cds):
    """从单个FASTA文件提取指定的CDS序列 - 优化版本"""
    sample_name = os.path.splitext(os.path.basename(fasta_file_path))[0]
    
    start = cds['start']
    end = cds['end']
    strand = cds['strand']
    
    try:
        with open(fasta_file_path, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                sequence = str(record.seq)
                
                # 检查索引是否有效
                if end >= len(sequence):
                    return sample_name, None
                
                # 提取序列
                if strand == '+':  # 正链
                    extracted_sequence = sequence[start:end + 1]
                elif strand == '-':  # 反链
                    extracted_sequence = get_complementary_sequence(sequence[start:end + 1])
                else:
                    return sample_name, None
                
                return sample_name, extracted_sequence
    except Exception as e:
        print(f"Error processing {fasta_file_path}: {e}")
        return sample_name, None
    
    return sample_name, None

def process_single_cds_worker(args):
    """单个CDS处理工作函数，用于并行处理"""
    cds_idx, cds, fasta_files, output_file = args
    
    cds_name = f"CDS_{cds_idx + 1}"
    processed_count = 0
    
    try:
        with open(output_file, 'w') as out_f:
            for i, fasta_file in enumerate(fasta_files):
                sample_name, sequence = extract_single_cds_from_fasta(fasta_file, cds)
                
                if sequence is not None:
                    processed_count += 1
                    out_f.write(f">{sample_name}\n")
                    # 每行60个字符
                    for j in range(0, len(sequence), 60):
                        out_f.write(sequence[j:j+60] + "\n")
                
                # 每处理100个文件打印进度
                if (i + 1) % 100 == 0:
                    print(f"  CDS {cds_idx + 1}: {i + 1}/{len(fasta_files)} files processed...")
        
        return cds_idx, cds_name, processed_count, True
    
    except Exception as e:
        print(f"Error processing CDS {cds_idx + 1}: {e}")
        return cds_idx, cds_name, 0, False

def process_all_fastas_parallel(input_folder, cds_file_path, output_folder, num_processes=None):
    """并行流式处理所有FASTA文件并生成对齐文件"""
    # 加载CDS信息
    cds_info = load_cds_info(cds_file_path)
    print(f"Loaded {len(cds_info)} CDS regions")
    
    # 获取所有FASTA文件
    fasta_files = glob.glob(os.path.join(input_folder, "*.fasta"))
    print(f"Found {len(fasta_files)} FASTA files")
    
    # 确保输出目录存在
    os.makedirs(output_folder, exist_ok=True)
    
    # 设置进程数
    if num_processes is None:
        num_processes = min(mp.cpu_count(), len(cds_info))
    print(f"Using {num_processes} processes")
    
    # 准备任务参数
    tasks = []
    for cds_idx, cds in enumerate(cds_info):
        cds_name = f"CDS_{cds_idx + 1}"
        output_file = os.path.join(output_folder, f"{cds_name}_merged.aln")
        tasks.append((cds_idx, cds, fasta_files, output_file))
    
    # 并行处理
    print("Starting parallel processing...")
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_single_cds_worker, tasks)
    
    # 汇总结果
    successful = 0
    for cds_idx, cds_name, sample_count, success in results:
        if success:
            successful += 1
            print(f"✓ {cds_name}: {sample_count} samples")
        else:
            print(f"✗ {cds_name}: FAILED")
    
    print(f"\nProcessing completed!")
    print(f"Successfully processed: {successful}/{len(cds_info)} CDS regions")

def main():
    parser = argparse.ArgumentParser(description="高效并行流式处理：直接从FASTA文件生成CDS对齐文件")
    parser.add_argument("--input_folder", required=True, help="输入FASTA文件夹路径")
    parser.add_argument("--cds_file", required=True, help="CDS位置文件路径")
    parser.add_argument("--output_folder", required=True, help="输出对齐文件夹路径")
    parser.add_argument("--processes", type=int, default=None, help="并行进程数 (默认自动选择)")
    
    args = parser.parse_args()
    
    process_all_fastas_parallel(args.input_folder, args.cds_file, args.output_folder, args.processes)

if __name__ == "__main__":
    main()
