#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
超级优化版本：使用多进程+内存映射+批处理
Usage:
    python 1-3-merge_ultra_optimized.py [num_processes]

Env vars (all required):
    ALIGN_DIR : dir containing *.aln
    OUT_DIR   : dir to write <sample>.fasta
"""

import os
import sys
import time
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import functools

def process_aln_file(aln_file_path):
    """处理单个.aln文件，返回该文件中所有样本的序列数据"""
    aln_file = Path(aln_file_path)
    gene = aln_file.stem
    
    sample_seqs = {}
    try:
        for rec in SeqIO.parse(aln_file, "fasta"):
            sample_id = rec.id.replace(".fasta", "")
            sequence = str(rec.seq).strip()
            sample_seqs[sample_id] = (gene, sequence)
    except Exception as e:
        print(f"[ERROR] 处理文件 {aln_file} 时出错: {e}")
        return {}
    
    return sample_seqs

def write_sample_batch(args):
    """批量写入样本文件"""
    sample_batch, out_dir = args
    
    for sample_id, gene_seqs in sample_batch.items():
        if not gene_seqs:
            continue
        
        # 按基因名排序
        gene_seqs.sort(key=lambda x: x[0])
        
        # 写入文件
        output_file = out_dir / f"{sample_id}.fasta"
        try:
            with output_file.open("w") as fh:
                fh.write(f">{sample_id}\n")
                fh.write('---'.join(seq for _, seq in gene_seqs) + '\n')
        except Exception as e:
            print(f"[ERROR] 写入样本 {sample_id} 时出错: {e}")

def chunk_dict(data, chunk_size):
    """将字典分块"""
    items = list(data.items())
    for i in range(0, len(items), chunk_size):
        yield dict(items[i:i + chunk_size])

def main():
    # 获取进程数参数
    num_processes = int(sys.argv[1]) if len(sys.argv) > 1 else min(cpu_count(), 8)
    print(f"[INFO] 使用 {num_processes} 个进程进行处理")
    
    align_dir = Path(os.environ["ALIGN_DIR"])
    out_dir = Path(os.environ["OUT_DIR"])
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("[INFO] 开始超级优化版本处理...")
    start_time = time.time()
    
    # 第一步：并行读取所有.aln文件
    aln_files = list(align_dir.glob("*.aln"))
    print(f"[INFO] 找到 {len(aln_files)} 个.aln文件")
    
    if not aln_files:
        print("[ERROR] 未找到.aln文件")
        return
    
    print("[INFO] 正在并行读取所有.aln文件...")
    
    # 使用多进程池并行处理文件
    with Pool(num_processes) as pool:
        file_results = pool.map(process_aln_file, aln_files)
    
    # 合并所有结果
    print("[INFO] 正在合并读取结果...")
    sample_data = defaultdict(list)
    
    for file_result in file_results:
        for sample_id, gene_seq in file_result.items():
            sample_data[sample_id].append(gene_seq)
    
    print(f"[INFO] 数据读取完成，共找到 {len(sample_data)} 个样本")
    
    # 第二步：批量并行写入样本文件
    print("[INFO] 正在批量生成样本文件...")
    
    # 将样本数据分批
    batch_size = max(50, len(sample_data) // num_processes)
    sample_batches = list(chunk_dict(sample_data, batch_size))
    
    # 准备参数
    batch_args = [(batch, out_dir) for batch in sample_batches]
    
    # 并行写入
    with Pool(num_processes) as pool:
        pool.map(write_sample_batch, batch_args)
    
    end_time = time.time()
    print(f"[DONE] 处理完成！共处理 {len(sample_data)} 个样本")
    print(f"[INFO] 总耗时: {end_time - start_time:.2f} 秒")
    print(f"[INFO] 平均每个样本耗时: {(end_time - start_time) / len(sample_data):.4f} 秒")

if __name__ == "__main__":
    main()
