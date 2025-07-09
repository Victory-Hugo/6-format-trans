#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
优化版本：一次性读取所有文件，批量处理所有样本
Usage:
    python 1-3-merge_optimized.py

Env vars (all required):
    ALIGN_DIR : dir containing *.aln
    OUT_DIR   : dir to write <sample>.fasta
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
import time

def main():
    align_dir = Path(os.environ["ALIGN_DIR"])
    out_dir = Path(os.environ["OUT_DIR"])
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("[INFO] 开始优化版本处理...")
    start_time = time.time()
    
    # 第一步：一次性读取所有.aln文件，构建数据结构
    print("[INFO] 正在读取所有.aln文件...")
    
    # 数据结构：sample_data[sample_id] = [(gene, seq), ...]
    sample_data = defaultdict(list)
    
    aln_files = list(align_dir.glob("*.aln"))
    print(f"[INFO] 找到 {len(aln_files)} 个.aln文件")
    
    for i, aln_file in enumerate(aln_files):
        gene = aln_file.stem  # 如: CDS_1_merged
        
        # 读取当前基因的所有样本序列
        for rec in SeqIO.parse(aln_file, "fasta"):
            sample_id = rec.id.replace(".fasta", "")
            sequence = str(rec.seq).strip()
            sample_data[sample_id].append((gene, sequence))
        
        if (i + 1) % 100 == 0:
            print(f"[INFO] 已处理 {i + 1}/{len(aln_files)} 个基因文件")
    
    print(f"[INFO] 数据读取完成，共找到 {len(sample_data)} 个样本")
    
    # 第二步：为每个样本生成输出文件
    print("[INFO] 正在生成样本文件...")
    
    processed_count = 0
    for sample_id, gene_seqs in sample_data.items():
        if not gene_seqs:
            print(f"[WARN] 样本 {sample_id} 没有找到序列")
            continue
        
        # 按基因名排序，确保输出顺序一致
        gene_seqs.sort(key=lambda x: x[0])
        
        # 写入文件
        output_file = out_dir / f"{sample_id}.fasta"
        with output_file.open("w") as fh:
            fh.write(f">{sample_id}\n")
            fh.write('---'.join(seq for _, seq in gene_seqs) + '\n')
        
        processed_count += 1
        if processed_count % 100 == 0:
            print(f"[INFO] 已处理 {processed_count}/{len(sample_data)} 个样本")
    
    end_time = time.time()
    print(f"[DONE] 处理完成！共处理 {processed_count} 个样本")
    print(f"[INFO] 总耗时: {end_time - start_time:.2f} 秒")

if __name__ == "__main__":
    main()
