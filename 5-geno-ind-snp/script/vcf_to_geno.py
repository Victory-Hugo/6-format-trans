#!/usr/bin/env python3
"""
VCF文件转换为geno, ind, snp格式 - Python高效版本
用于处理大型VCF文件，适用于ADMIXTOOLS/f3统计分析
"""

import gzip
import sys
import os
from pathlib import Path

def process_vcf_to_geno_format(input_vcf, output_dir):
    """
    将VCF文件转换为geno、ind、snp格式
    """
    print(f"开始处理VCF文件: {input_vcf}")
    
    # 检查输入文件
    if not os.path.exists(input_vcf):
        raise FileNotFoundError(f"输入VCF文件不存在: {input_vcf}")
    
    # 创建输出目录
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # 输出文件路径
    geno_file = os.path.join(output_dir, "hp_7544.geno")
    ind_file = os.path.join(output_dir, "hp_7544.ind")
    snp_file = os.path.join(output_dir, "hp_7544.snp")
    
    # 读取VCF文件头部
    print("读取VCF头部信息...")
    with gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, 'r') as f:
        sample_names = []
        for line in f:
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                sample_names = headers[9:]  # 样本名从第10列开始
                break
            elif not line.startswith('#'):
                raise ValueError("未找到VCF头部信息")
    
    n_samples = len(sample_names)
    print(f"样本数: {n_samples}")
    
    # 创建ind文件
    print("创建ind文件...")
    with open(ind_file, 'w') as f:
        for sample in sample_names:
            f.write(f"{sample}\tU\t{sample}\n")
    
    # 处理VCF数据并创建geno和snp文件
    print("开始处理VCF数据...")
    
    processed_snps = 0
    
    with (gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, 'r')) as vcf_f, \
         open(geno_file, 'w') as geno_f, \
         open(snp_file, 'w') as snp_f:
        
        for line in vcf_f:
            # 跳过头部行
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            # 提取SNP信息
            chrom = fields[0]
            pos = fields[1]
            ref_allele = fields[3]
            alt_allele = fields[4]
            
            # 提取基因型（从第10列开始）
            genotypes = fields[9:]
            
            # 转换基因型编码
            geno_codes = []
            for gt_info in genotypes:
                gt = gt_info.split(':')[0]  # 提取GT字段
                
                if gt in ['0/0', '0|0']:
                    geno_codes.append('2')  # 纯合子参考
                elif gt in ['0/1', '1/0', '0|1', '1|0']:
                    geno_codes.append('1')  # 杂合子
                elif gt in ['1/1', '1|1']:
                    geno_codes.append('0')  # 纯合子替换
                else:
                    geno_codes.append('9')  # 缺失或其他
            
            # 写入geno文件（每行是一个SNP）
            geno_f.write(''.join(geno_codes) + '\n')
            
            # 写入snp文件
            snp_id = f"snp_{pos}"
            genetic_dist = "0"  # 遗传距离设为0
            snp_f.write(f"{snp_id}\t{chrom}\t{genetic_dist}\t{pos}\t{ref_allele}\t{alt_allele}\n")
            
            processed_snps += 1
            
            # 显示进度
            if processed_snps % 10000 == 0:
                print(f"已处理 {processed_snps} 个SNP")
    
    print(f"\n处理完成! 总共处理了 {processed_snps} 个SNP")
    
    # 验证输出文件
    print("\n验证输出文件:")
    for file_path, file_type in [(geno_file, "geno"), (ind_file, "ind"), (snp_file, "snp")]:
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                lines = sum(1 for _ in f)
            size_mb = os.path.getsize(file_path) / (1024 * 1024)
            print(f"{file_type}文件: {lines} 行, {size_mb:.2f} MB")
    
    print(f"\n输出文件:")
    print(f"- geno文件: {geno_file}")
    print(f"- ind文件: {ind_file}")
    print(f"- snp文件: {snp_file}")
    
    print("\n注意事项:")
    print("1. geno文件编码: 0=纯合子替换, 1=杂合子, 2=纯合子参考, 9=缺失")
    print("2. ind文件中所有样本性别设为'U'(未知)，群体标签为样本名")
    print("3. snp文件中遗传距离设为0，可根据需要修改")
    print("4. 如需修改群体标签，请编辑ind文件的第三列")

if __name__ == "__main__":
    # 输入和输出路径
    input_vcf = "/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544.NoN.SNP.maf99.vcf.recode.vcf.gz"
    output_dir = "/mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/5-f3/script/"
    
    try:
        process_vcf_to_geno_format(input_vcf, output_dir)
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)
