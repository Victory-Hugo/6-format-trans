#!/usr/bin/env python3
# /doc: 将 FASTA 序列和样本分组信息转换为 .arp 格式，用于遗传多样性分析

from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Iterator


def read_file_lines(file_path: Path, encoding: str = 'utf-8') -> List[str]:
    """
    读取文本文件，自动尝试 utf-8 和 latin-1 编码。

    :param file_path: 要读取的文件路径
    :param encoding: 首选编码，默认为 'utf-8'
    :return: 文件的所有行列表（包含换行符）
    """
    try:
        return file_path.read_text(encoding=encoding).splitlines(keepends=True)
    except UnicodeDecodeError:
        return file_path.read_text(encoding='latin-1').splitlines(keepends=True)


def parse_group(lines: List[str]) -> Dict[str, List[str]]:
    """
    解析分组文件，每行格式：样本名<TAB>组名

    :param lines: 分组文件的行列表
    :return: 以组名为键，样本名列表为值的字典
    """
    groups = defaultdict(list)
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) != 2:
            continue  # 跳过格式不正确的行
        sample, group = parts
        groups[group].append(sample)
    return groups


def parse_fasta(lines: List[str]) -> Dict[str, str]:
    """
    解析 FASTA 格式，将每个样本的序列合并为一行字符串。

    :param lines: FASTA 文件的行列表
    :return: 以样本名为键，序列字符串为值的字典
    """
    sequences: Dict[str, List[str]] = {}
    current_sample = None

    for line in lines:
        if line.startswith('>'):
            # 提取样本名（去掉 '>'，空格后只取第一个字段）
            current_sample = line[1:].strip().split()[0]
            sequences[current_sample] = []
        elif current_sample:
            # 累加序列行（去除末尾换行）
            sequences[current_sample].append(line.strip())

    # 合并为单行序列
    return {s: ''.join(seq_lines) for s, seq_lines in sequences.items()}


def generate_arp_content(groups: Dict[str, List[str]],
                         sequences: Dict[str, str],
                         default_len: int = 100) -> Iterator[str]:
    """
    根据分组和序列信息，生成 arp 文件的文本内容迭代器。

    :param groups: 组名到样本列表的映射
    :param sequences: 样本名到序列字符串的映射
    :param default_len: 若样本缺失序列，使用多少个 'N' 填充
    :yield: 每一行 arp 文件内容
    """
    yield "[Profile]"
    yield '   Title = "Genetic Diversity Analysis"'
    yield f"   NbSamples = {len(groups)}"
    yield "   DataType = DNA"
    yield "   GenotypicData = 0"
    yield "   LocusSeparator = NONE"
    yield '   MissingData = "N"'
    yield "   CompDistMatrix = 1"
    yield ""
    yield "[Data]"
    yield ""

    for group_name, samples in groups.items():
        yield "[[Samples]]"
        yield f'   SampleName = "{group_name}"'
        yield f"   SampleSize = {len(samples)}"
        yield "   SampleData= {"
        for sample in samples:
            seq = sequences.get(sample, "N" * default_len)
            yield f"       {sample} 1 {seq}"
        yield "   }"
        yield ""


def write_arp_file(output_path: Path,
                   content_iter: Iterator[str]) -> None:
    """
    将 arp 内容写入文件。

    :param output_path: arp 文件输出路径
    :param content_iter: arp 文件内容的行迭代器
    """
    with output_path.open('w', encoding='utf-8') as f:
        for line in content_iter:
            f.write(line + '\n')


def main():
    # 配置文件路径
    fasta_path = Path(r"/mnt/c/Users/Administrator/Desktop/Affey_chip_mtDNA.fasta") # todo fasta文件注意对齐
    group_path = Path(r"/mnt/c/Users/Administrator/Desktop/1.txt") # todo txt文件第一列为ID，第二列为组名
    arp_path = Path(r"/mnt/c/Users/Administrator/Desktop/This_study.arp")

    # 读取文件
    fasta_lines = read_file_lines(fasta_path)
    group_lines = read_file_lines(group_path)

    # 解析数据
    group_dict = parse_group(group_lines)
    seq_dict = parse_fasta(fasta_lines)

    # 生成并写入 .arp 文件
    arp_iter = generate_arp_content(group_dict, seq_dict)
    write_arp_file(arp_path, arp_iter)


if __name__ == "__main__":
    main()
