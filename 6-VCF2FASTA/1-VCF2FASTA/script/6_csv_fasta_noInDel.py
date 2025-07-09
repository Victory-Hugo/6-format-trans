import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

def process_variants_and_generate_fasta(csv_file: str, fasta_file: str, output_dir: str):
    """
    根据变异数据（仅支持多碱基替换）和参考基因组生成修改后的FASTA文件。

    参数:
        csv_file (str): 包含变异信息的CSV文件路径。
                        要求至少包含: P1(1-based坐标), SUB_REF, SUB_ALT。
        fasta_file (str): 参考基因组的FASTA文件路径。
        output_dir (str): 输出目录，用于保存生成的FASTA文件。
    """
    # 1. 读取CSV文件
    data = pd.read_csv(csv_file)

    # 2. 将 SUB_REF 中的 '.' 替换为 '-'，并将 SUB_ALT 中的 '.' 替换为 'N'
    data['SUB_REF'] = data['SUB_REF'].replace('.', '-')
    data['SUB_ALT'] = data['SUB_ALT'].replace('.', 'N')  # 假设 '.' 表示未知碱基，用 'N' 替代

    # 3. 按照变异位点从小到大排序（以1-based坐标升序排列）
    data.sort_values(by='P1', inplace=True)

    # 4. 读取参考基因组FASTA文件
    reference_record = SeqIO.read(fasta_file, "fasta")
    reference_seq = str(reference_record.seq)  # 转为字符串便于拼接

    # 5. 开始分段拼接
    output_seq_fragments = []  # 用于存放每段拼接结果
    current_pos = 0            # 指示当前在参考序列中处理到的位置(0-based)

    for _, row in data.iterrows():
        # a) 提取坐标、参考、替换序列信息
        pos_1_based = int(row['P1'])
        pos_0_based = pos_1_based - 1  # 转为0-based
        sub_ref = str(row['SUB_REF'])
        sub_alt = str(row['SUB_ALT'])

        # b) 如果 SUB_REF 是 '.'，替换为 '-'
        if sub_ref == '.':
            sub_ref = '-'

        # c) 如果 SUB_REF 和 SUB_ALT 的长度不一致，则忽略该变异
        if len(sub_ref) != len(sub_alt):
            print(f"Ignored variant at position {pos_1_based}: length mismatch (SUB_REF: {sub_ref}, SUB_ALT: {sub_alt})")
            continue

        # d) 将[ current_pos, pos_0_based )这段（未受此变异影响的部分）拼接进输出
        if pos_0_based > current_pos:
            output_seq_fragments.append(reference_seq[current_pos:pos_0_based])

        # e) 进行替换
        output_seq_fragments.append(sub_alt)

        # f) 更新 current_pos
        current_pos = pos_0_based + len(sub_ref)

    # 6. 若 current_pos 还没到参考序列末尾，则把剩余的部分拼接到末尾
    if current_pos < len(reference_seq):
        output_seq_fragments.append(reference_seq[current_pos:])

    # 7. 生成最终的新序列
    modified_sequence = ''.join(output_seq_fragments)

    # 8. 检查并替换非 A/G/C/T/N 的碱基为 N
    allowed_bases = {'A', 'G', 'C', 'T', 'N'}
    modified_sequence_checked = ''.join(
        base if base.upper() in allowed_bases else 'N'
        for base in modified_sequence
    )

    # 9. 创建新的FASTA记录
    csv_basename = os.path.splitext(os.path.basename(csv_file))[0]
    modified_record = SeqRecord(
        Seq(modified_sequence_checked),
        id=csv_basename,
        description=""
    )

    # 10. 保存修改后的FASTA文件
    os.makedirs(output_dir, exist_ok=True)
    output_fasta_file = os.path.join(output_dir, f"{csv_basename}.fasta")
    SeqIO.write(modified_record, output_fasta_file, "fasta")
    print(f"Modified FASTA file has been saved as: {output_fasta_file}")

    return output_fasta_file

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_csv_file> <reference_fasta_file> <output_fasta_directory>")
        sys.exit(1)

    input_csv = sys.argv[1]
    reference_fasta = sys.argv[2]
    output_dir = sys.argv[3]

    process_variants_and_generate_fasta(input_csv, reference_fasta, output_dir)
