#!/usr/bin/env python3
# /doc: 根据分组文件和 FASTA 序列，直接生成包含 Traits 区块的 NEXUS (.nex) 文件，无需中间文件

from pathlib import Path
from typing import Dict
import pandas as pd


def read_fasta(fasta_path: Path) -> Dict[str, str]:
    """
    解析 FASTA 文件，将每个样本的序列合并为单行字符串。

    :param fasta_path: FASTA 文件路径
    :return: 字典，键为样本 ID，值为序列字符串
    """
    seqs: Dict[str, str] = {}
    with fasta_path.open('r', encoding='utf-8') as fh:
        current_id = None
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                current_id = line[1:].strip()
                seqs[current_id] = ''
            elif current_id:
                seqs[current_id] += line
    return seqs


def write_id_haplogroup_fasta(seqs: Dict[str, str],
                              haplo_map: Dict[str, str],
                              out_path: Path) -> None:
    """
    将样本 ID 与单倍群拼接后输出新的 FASTA 文件。

    :param seqs: 原始样本序列字典，键为样本 ID，值为序列
    :param haplo_map: 样本 ID 到单倍群的映射
    :param out_path: 输出的新 FASTA 文件路径
    """
    with out_path.open('w', encoding='utf-8') as fh:
        for sid, seq in seqs.items():
            hap = haplo_map.get(sid, 'Unknown')
            fh.write(f">{sid}_{hap}\n{seq}\n")


def write_nexus_with_traits(fasta_with_hap: Path,
                            metadata_df: pd.DataFrame,
                            nex_path: Path) -> None:
    """
    根据带单倍群的 FASTA 文件和分组 metadata，生成包含 Data 与 Traits 区块的 NEXUS 文件。

    :param fasta_with_hap: 带单倍群后缀的 FASTA 文件路径
    :param metadata_df: 包含 ['SampleID','Region','Haplogroup'] 的 DataFrame
    :param nex_path: 输出的 NEXUS 文件路径
    """
    # 重新解析带单倍群 FASTA
    seqs = read_fasta(fasta_with_hap)
    # 计算 ntax 与 nchar
    ntax = len(seqs)
    nchar = len(next(iter(seqs.values()))) if ntax > 0 else 0

    # 构造 NewID 与 Region 映射
    metadata_df['NewID'] = metadata_df['SampleID'].astype(str) + '_' + metadata_df['Haplogroup'].astype(str)
    region_map = dict(zip(metadata_df['NewID'], metadata_df['Region']))
    regions = sorted({r for r in region_map.values()})
    # Traits 区块参数
    ntraits = len(regions)
    trait_labels = ' '.join(regions)

    with nex_path.open('w', encoding='utf-8') as fh:
        # 写入 Data 区块
        fh.write("#NEXUS\n\n")
        fh.write("Begin Data;\n")
        fh.write(f"\tDimensions ntax={ntax} NCHAR={nchar};\n")
        fh.write("\tFormat datatype=DNA missing=N GAP=-;\n")
        fh.write("\tMATRIX\n")
        for newid, seq in seqs.items():
            fh.write(f"\t{newid}\n{seq}\n")
        fh.write(";\nEND;\n\n")
        # 写入 Traits 区块
        fh.write("Begin Traits;\n")
        fh.write(f"Dimensions NTraits={ntraits};\n")
        fh.write("Format labels=yes missing=? separator=Comma;\n")
        fh.write(f"TraitLabels {trait_labels};\n")
        fh.write("Matrix\n")
        for newid in seqs.keys():
            # 样本对各 Region 的二值标记
            membership = ['1' if region_map.get(newid)==r else '0' for r in regions]
            fh.write(f"{newid} {','.join(membership)}\n")
        fh.write(";\nend;\n")


def main():
    # ——— 请将以下路径替换为您本地实际路径 ——— #
    desktop = Path(r"C:\Users\a\Desktop")
    csv_path   = desktop / "Illumina芯片分析表.csv"      #todo 原始分组表，包含 SampleID, Region, Haplogroup
    fasta_in   = desktop / "extracted_sequences.fasta"  #todo 原 FASTA 文件
    fasta_out  = desktop / "ID_haplogroup.fasta"        # 带单倍群后缀的 FASTA
    nexus_out  = desktop / "final.nex"                  # 生成的 NEXUS 文件

    # 1. 读取分组表，构建 SampleID->Haplogroup 映射
    df = pd.read_csv(csv_path, header=None,
                     names=['SampleID','Region','Haplogroup'],
                     encoding='utf-8')
    haplo_map = dict(zip(df['SampleID'].astype(str), df['Haplogroup'].astype(str)))

    # 2. 读取原始 FASTA
    seqs = read_fasta(fasta_in)
    # 3. 输出带单倍群的 FASTA
    write_id_haplogroup_fasta(seqs, haplo_map, fasta_out)
    # 4. 生成包含 Data 与 Traits 的 NEXUS 文件
    write_nexus_with_traits(fasta_out, df, nexus_out)

    print(f"NEXUS 文件已生成：{nexus_out}")


if __name__ == "__main__":
    main()
