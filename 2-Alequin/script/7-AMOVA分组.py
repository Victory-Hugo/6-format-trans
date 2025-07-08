#!/usr/bin/env python3
# /doc: 根据分组文件生成 .arp 格式的 Structure 定义文件，包含组数和各组成员 ID

from pathlib import Path
from typing import Dict, List
import pandas as pd


def load_groups(file_path: Path) -> Dict[str, List[str]]:
    """
    读取制表符分隔的分组文件，返回 Category -> [ID,...] 的映射。

    :param file_path: 分组文件路径，每行格式 "ID<TAB>Category"
    :return: 字典，键为分类名，值为样本 ID 列表
    """
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['ID', 'Category'], encoding='utf-8')
    return df.groupby('Category')['ID'].apply(list).to_dict()


def build_arp_lines(groups: Dict[str, List[str]], structure_name: str) -> List[str]:
    """
    构建 ARP 文件内容行列表。

    :param groups: Category -> [ID 列表]
    :param structure_name: StructureName 字段值
    :return: ARP 文件内容的行列表
    """
    lines: List[str] = []
    # Structure 段头
    lines.append('[[Structure]]')
    lines.append('')
    lines.append(f'StructureName="{structure_name}"')
    lines.append(f'NbGroups={len(groups)}')
    lines.append('')
    # 各组定义
    for ids in groups.values():
        lines.append('Group={')
        for sample_id in ids:
            lines.append(f'\t"{sample_id}"')
        lines.append('}')
    return lines


def write_arp_file(lines: List[str], output_path: Path) -> None:
    """
    将 ARP 内容写入指定文件。

    :param lines: ARP 文件内容行列表
    :param output_path: 输出文件路径
    """
    content = '\n'.join(lines) + '\n'
    output_path.write_text(content, encoding='utf-8')


def main():
    # —— 请根据实际情况修改以下路径 —— #
    desktop = Path("C:/Users/victo/Desktop")
    input_path = desktop / "新建 Text Document.txt" #todo 第一列是群体名，第二列是Group名
    output_path = desktop / "分组.arp" #todo 生成了arp之后直接将这个arp文件放到放到原始的arp末尾即可
    structure_name = "New Edited Structure"

    # 加载分组并生成 ARP 内容
    groups = load_groups(input_path)
    arp_lines = build_arp_lines(groups, structure_name)
    write_arp_file(arp_lines, output_path)

    print(f"ARP 文件已生成：{output_path}")


if __name__ == "__main__":
    main()
