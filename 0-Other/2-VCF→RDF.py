#!/usr/bin/env python3
# /doc: 从 VCF 文件生成 mtDNA FASTA 并组合生成中间 TXT，再重命名序列输出最终 RDF 和映射关系 CSV，处理流程无需人工干预

from pathlib import Path
from typing import List, Dict, Tuple
import pandas as pd
import csv


def load_id_mapping(csv_path: Path) -> Dict[str, str]:
    """
    加载 CSV 映射表，将 Object_ID 映射到包含省份和单倍群的附加信息字符串。

    :param csv_path: 包含 Object_ID, Re_Population_Province, Haplogrouper 列的 CSV 文件路径
    :return: { Object_ID: ";Province;Haplogrouper" }
    """
    df = pd.read_csv(csv_path, dtype=str)
    return {
        row['Object_ID']: f";{row['Re_Population_Province']};{row['Haplogrouper']}"
        for _, row in df.iterrows()
    }


def parse_vcf(vcf_path: Path, id_to_info: Dict[str, str]) -> Tuple[List[str], Dict[str, List[str]]]:
    """
    解析 VCF 文件，按过滤后的样本生成碱基序列，并收集所有变异位点位置。

    :param vcf_path: 输入 VCF 文件路径
    :param id_to_info: 只保留在此映射中的样本名，并在 FASTA 头中附加额外信息
    :return: 
      - positions: 按出现顺序的位点列表（CHROM 列第2列）
      - sample_sequences: { sample_name: [">name;info;", base1, base2, ...] }
    """
    sample_sequences: Dict[str, List[str]] = {}
    positions: List[str] = []
    sample_names: List[str] = []

    with vcf_path.open('r', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                cols = line.rstrip('\n').split('\t')
                all_names = cols[9:]
                # 仅保留映射表中存在的样本
                sample_names = [n for n in all_names if n in id_to_info]
                # 初始化每个样本的序列列表（包含描述行）
                for n in sample_names:
                    sample_sequences[n] = [f">{n}{id_to_info[n]};\n"]
                continue

            # 数据行：提取位置、参考和替代碱基
            parts = line.rstrip('\n').split('\t')
            pos = parts[1]
            positions.append(pos)
            ref = parts[3]
            alt_list = [] if parts[4] in ('.', 'X') else parts[4].split(',')

            # 遍历每个样本，按基因型添加碱基
            for name, field in zip(sample_names, parts[9:]):
                gt = field.split(':', 1)[0]
                if '.' in gt:
                    base = 'N'
                elif gt in ('0/0', '0'):
                    base = ref
                else:
                    alleles = gt.replace('|', '/').split('/')
                    # 杂合或不一致 -> N
                    if len(set(alleles)) > 1:
                        base = 'N'
                    else:
                        try:
                            idx = int(alleles[0]) - 1
                            base = alt_list[idx]
                        except (ValueError, IndexError):
                            base = 'N'
                sample_sequences[name].append(base)

    return positions, sample_sequences


def write_fasta(sample_sequences: Dict[str, List[str]], fasta_path: Path) -> None:
    """
    将 sample_sequences 写入 FASTA 文件。

    :param sample_sequences: { name: [header, base1, base2, ...] }
    :param fasta_path: 输出 FASTA 路径
    """
    with fasta_path.open('w', encoding='utf-8') as out:
        for seq_list in sample_sequences.values():
            out.write(''.join(seq_list) + "\n")


def write_txt(positions: List[str], sample_sequences: Dict[str, List[str]], txt_path: Path) -> None:
    """
    生成中间 TXT 文件，包含位置列表、常数行 “10”，以及 FASTA 序列。

    :param positions: 变异位点位置列表
    :param sample_sequences: FASTA 序列 dict
    :param txt_path: 中间 TXT 输出路径
    """
    count = len(positions)
    with txt_path.open('w', encoding='utf-8') as out:
        out.write("  ;1.0\n")
        out.write(';'.join(positions) + ';\n')
        out.write(';'.join(['10'] * count) + ';\n\n')
        # 追加 FASTA 序列
        for seq_list in sample_sequences.values():
            out.write(''.join(seq_list) + "\n")


def build_mapping_and_rdf(txt_path: Path, rdf_path: Path, mapping_csv: Path) -> None:
    """
    从中间 TXT 读取描述行，为每条序列生成新名称，写入映射 CSV 并输出最终 RDF 文件。

    :param txt_path: 包含描述行和序列的中间 TXT 文件
    :param rdf_path: 输出最终 RDF 文件路径
    :param mapping_csv: 输出映射关系 CSV 路径
    """
    mapping: List[Tuple[str, str]] = []
    new_lines: List[str] = []
    seq_counter = 1

    # 读取并重命名
    with txt_path.open('r', encoding='latin-1') as inp:
        for line in inp:
            if line.startswith('>'):
                orig = line[1:].split(';', 1)[0]
                # 原始 ID 第一字符 + 下划线后第一字符 + 序号 组成新名称
                parts = orig.split('_')
                second = parts[1][0] if len(parts) > 1 and parts[1] else ''
                newname = f"{orig[0]}{second}{seq_counter}"
                seq_counter += 1
                mapping.append((orig, newname))

                # 在第一个分号后插入 ";1;" 并替换为新名称
                header_rest = line.split(';', 1)[1]
                line = f">{newname};1;{header_rest}"
            new_lines.append(line)

    # 写映射关系 CSV
    with mapping_csv.open('w', newline='', encoding='utf-8') as mcf:
        writer = csv.writer(mcf)
        writer.writerow(['Original Name', 'New Name'])
        writer.writerows(mapping)

    # 写最终 RDF
    rdf_path.write_text(''.join(new_lines), encoding='utf-8')


def main():
    # ——— 请根据实际情况修改以下路径 ——— #
    desk = Path.home() / "Desktop"
    vcf_file = desk / "Illumina_mtDNA.vcf"
    info_csv = desk / "Illumina芯片分析表.csv"
    fasta_tmp = desk / "temp_mtDNA.fasta"
    txt_tmp   = desk / "temp_variants.txt"
    rdf_final = desk / "Final.rdf"
    map_csv   = desk / "映射关系.csv"

    # 1. 加载 ID 到附加信息的映射
    id_map = load_id_mapping(info_csv)

    # 2. 解析 VCF，获取位点列表和样本序列
    positions, sequences = parse_vcf(vcf_file, id_map)
    print("过滤规则：含 “.” 或 杂合 → N；0/0 → 参考碱基；纯合替代 → 替代碱基")

    # 3. 写临时 FASTA 并生成中间 TXT
    write_fasta(sequences, fasta_tmp)
    write_txt(positions, sequences, txt_tmp)

    # 4. 删除临时 FASTA
    fasta_tmp.unlink(missing_ok=True)

    # 5. 生成映射关系 CSV 与最终 RDF
    build_mapping_and_rdf(txt_tmp, rdf_final, map_csv)

    # 6. 删除临时 TXT
    txt_tmp.unlink(missing_ok=True)

    print(f"所有文件已处理完成：\n- 最终 RDF：{rdf_final}\n- 映射 CSV：{map_csv}")


if __name__ == '__main__':
    main()
