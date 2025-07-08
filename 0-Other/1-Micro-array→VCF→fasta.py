#!/usr/bin/env python3
# /doc: 本脚本将 mtDNA 的 VCF 文件转换为 FASTA 格式，视为单倍体并按照规则最小化 N 的使用

from pathlib import Path
from typing import List, Dict


def get_base_from_genotype(gt: str, ref: str, alt: List[str]) -> str:
    """
    根据基因型字符串返回对应的碱基：
      - 含 “.” 或 未知 则返回 N
      - 纯合参考 (0/0 或 0) 返回参考碱基
      - 纯合替代 (1/1、2/2… 或 单个数字) 返回对应的替代碱基
      - 杂合或多等位基因不一致 返回 N

    :param gt: 基因型字段 (如 "0/0","0/1","1","./.")
    :param ref: 参考碱基
    :param alt: 替代碱基列表 (如 ["A","T"])，若无替代填空列表
    :return: 单个碱基字符
    """
    if "." in gt:
        return "N"
    # 将可能的 | 替换为 / 并拆分等位基因
    alleles = gt.replace("|", "/").split("/")
    # 全部为参考
    if all(a == "0" for a in alleles):
        return ref
    # 杂合则 N
    if len(set(alleles)) > 1:
        return "N"
    # 纯合非零
    try:
        idx = int(alleles[0]) - 1
        return alt[idx]
    except (ValueError, IndexError):
        return "N"


def vcf_to_fasta_mtDNA(vcf_path: Path, output_path: Path) -> None:
    """
    将 VCF 转换为 FASTA 格式：
      - 读取 #CHROM 行获得样本名
      - 逐行解析变异：参考 vs 替代 vs 杂合 vs 缺失
      - 按样本拼接序列并输出 FASTA

    :param vcf_path: 输入 VCF 文件路径
    :param output_path: 输出 FASTA 文件路径
    """
    sample_names: List[str] = []
    sequences: Dict[str, List[str]] = {}

    with vcf_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                # 提取样本名并初始化序列容器
                cols = line.strip().split("\t")
                sample_names = cols[9:]
                for name in sample_names:
                    sequences[name] = [f">{name}\n"]
                continue

            cols = line.strip().split("\t")
            ref_base = cols[3]
            # 如果 ALT 列为 'X' 或 '.' 则当无替代
            alt_bases = [] if cols[4] in ("X", ".") else cols[4].split(",")

            # 遍历每个样本的基因型字段
            for name, field in zip(sample_names, cols[9:]):
                gt = field.split(":", 1)[0]
                base = get_base_from_genotype(gt, ref_base, alt_bases)
                sequences[name].append(base)

    # 打印过滤规则说明
    print("过滤规则：含 '.' 或 杂合 -> N；0/0 -> 参考碱基；纯合替代 -> 替代碱基")

    # 写入 FASTA 文件
    with output_path.open("w", encoding="utf-8") as out_f:
        for name, seq_list in sequences.items():
            out_f.write("".join(seq_list) + "\n")


def main():
    # ——— 请将以下路径替换为您的实际文件路径 ——— #
    vcf_file    = Path(r"C:/Users/victo/Desktop/98pops10220indAffy_mtDNA_modified_AC0.vcf") # todo芯片的VCF文件
    fasta_file  = Path(r"C:/Users/victo/Desktop/Affey_chip_mtDNA.fasta")

    vcf_to_fasta_mtDNA(vcf_file, fasta_file)
    print(f"FASTA 文件已生成：{fasta_file}")


if __name__ == "__main__":
    main()
