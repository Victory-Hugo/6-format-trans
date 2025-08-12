#!/usr/bin/env python3
# /doc: 根据样本单倍群信息生成 mtDNA 频率分析所需的 ARP 文件

from pathlib import Path
from typing import Dict, Tuple
import pandas as pd


def load_data(file_path: Path) -> pd.DataFrame:
    """
    加载制表符分隔的输入文件，文件应包含 'name' 和 'haplogroup' 列。

    :param file_path: 原始数据文件路径
    :return: DataFrame
    """
    return pd.read_csv(file_path, sep='\t', encoding='utf-8')


def build_haplotype_dict(haplogroups: pd.Index) -> Dict[str, int]:
    """
    为每个唯一的单倍群分配一个整数编号，从 1 开始。

    :param haplogroups: 单倍群名称的唯一列表（Index）
    :return: 单倍群到编号的映射
    """
    return {hg: idx + 1 for idx, hg in enumerate(haplogroups)}


def compute_counts_and_freq(data: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    计算每个组（name）中各单倍群的出现次数及其频率。

    :param data: 原始 DataFrame，含 'name','haplogroup'
    :return: (counts, freqs) 二元组：
        - counts: 行是组名，列是单倍群，值为计数
        - freqs: 同维度，值为频率（计数/组内总数）
    """
    counts = data.groupby(['name', 'haplogroup']).size().unstack(fill_value=0)
    freqs = counts.div(counts.sum(axis=1), axis=0)
    return counts, freqs


def generate_arp_content(
    counts: pd.DataFrame,
    freqs: pd.DataFrame,
    hapl_dict: Dict[str, int],
    title: str
) -> str:
    """
    根据计数表、频率表和单倍群编号映射，生成完整的 ARP 文件内容。

    :param counts: 计数表
    :param freqs: 频率表
    :param hapl_dict: 单倍群到编号的映射
    :param title: Profile 中的 Title 字符串
    :return: ARP 文件的完整文本
    """
    nb_samples = counts.shape[0]
    lines = []

    # —— Profile 段 —— #
    lines.append("[Profile]")
    lines.append(f'    Title="{title}"')
    lines.append(f"    NbSamples={nb_samples}")
    lines.append("    DataType=FREQUENCY")
    lines.append("    GenotypicData=0")
    lines.append("    LocusSeparator=' '")
    lines.append('    MissingData="?"')
    lines.append("    Frequency=REL")
    lines.append("")

    # —— Data 段 —— #
    lines.append("[Data]")
    lines.append("")

    # —— HaplotypeDefinition 段 —— #
    lines.append("[HaplotypeDefinition]")
    lines.append("    HaplList={")
    for hap, idx in hapl_dict.items():
        lines.append(f"        {idx}\t{hap}")
    lines.append("    }")
    lines.append("")

    # —— Samples 段 —— #
    lines.append("[Samples]")
    for group in counts.index:
        sample_size = int(counts.loc[group].sum())
        lines.append(f'    SampleName="{group}"')
        lines.append(f"    SampleSize={sample_size}")
        lines.append("    SampleData={")
        for hap, idx in hapl_dict.items():
            freq = freqs.at[group, hap] if hap in freqs.columns else 0.0
            lines.append(f"        {idx}\t{freq:.6f}")
        lines.append("    }")
        lines.append("")

    return "\n".join(lines)


def write_arp(file_path: Path, content: str) -> None:
    """
    将生成的 ARP 文本写入文件。

    :param file_path: 输出 ARP 文件路径
    :param content: ARP 文件内容
    """
    file_path.write_text(content, encoding='utf-8')


def main():
    # —— 请将以下路径替换为您的实际路径 —— #
    input_file = Path(r"ID_hap.txt") #todo 输入文件，制表符分隔，包含 'name' 和 'haplogroup' 列
    output_file = Path(r"mtDNA.arp") #todo 输出 ARP 文件路径
    title_text = "The population fixation index(Fst) of mtDNA"

    # 1. 加载数据
    df = load_data(input_file)

    # 2. 构建单倍群字典
    hapl_dict = build_haplotype_dict(df['haplogroup'].unique())

    # 3. 计算计数与频率
    counts, freqs = compute_counts_and_freq(df)

    # 4. 生成 ARP 文本
    arp_text = generate_arp_content(counts, freqs, hapl_dict, title_text)

    # 5. 写入文件
    write_arp(output_file, arp_text)

    print(f"ARP 文件已生成：{output_file}")


if __name__ == "__main__":
    main()
