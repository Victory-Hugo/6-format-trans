import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="提取核心CDS对齐序列。")
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='输入对齐文件的路径（例如，merged_CDS_1.aln）。'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='输出核心对齐文件的路径（例如，core_CDS_1.aln）。'
    )
    parser.add_argument(
        '-f', '--format',
        default='fasta',
        help='对齐文件的格式（默认: fasta）。'
    )
    parser.add_argument(
        '-t', '--threshold',
        type=float,
        default=0.99,
        help='核心位置的阈值比例（默认: 0.99，即99%%）。'
    )
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    try:
        # 读取对齐文件
        alignment = AlignIO.read(args.input, args.format)
    except Exception as e:
        sys.exit(f"无法读取输入文件: {e}")
    
    num_sequences = len(alignment)
    threshold = num_sequences * args.threshold  # 动态阈值
    
    core_positions = []
    alignment_length = alignment.get_alignment_length()
    
    # 假设每3列为一个密码子
    for i in range(0, alignment_length, 3):
        codons = [str(record.seq[i:i+3]).upper() for record in alignment]
        # 统计缺失和终止密码子
        missing = codons.count('NNN')
        stop_codons = sum(codon in ['TAG', 'TAA', 'TGA'] for codon in codons)
        
        if missing <= num_sequences / 100 and stop_codons == 0:
            core_positions.extend([i, i+1, i+2])
    
    if not core_positions:
        sys.exit("未找到符合条件的核心位置。")
    
    # 提取核心CDS
    core_alignment = MultipleSeqAlignment([
        SeqRecord(
            Seq("".join(str(record.seq[pos]) for pos in core_positions)),
            id=record.id,
            description=record.description
        )
        for record in alignment
    ])
    
    try:
        # 写入核心CDS对齐文件
        AlignIO.write(core_alignment, args.output, args.format)
        print(f"核心CDS对齐文件已成功写入: {args.output}")
    except Exception as e:
        sys.exit(f"无法写入输出文件: {e}")

if __name__ == "__main__":
    main()
