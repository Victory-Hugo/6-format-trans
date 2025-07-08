#!/usr/bin/env python3
# /doc: 将 FASTA 文件中所有描述行（以 '>' 开头）的前缀替换为 '>gi|'，并写入新文件

from pathlib import Path


def replace_fasta_prefix(input_path: Path, output_path: Path) -> None:
    """
    读取 input_path 指定的 FASTA 文件，将每个描述行 '>' 替换为 '>gi|'，
    并将结果写入 output_path 指定的新 FASTA 文件。

    :param input_path: 原 FASTA 文件路径
    :param output_path: 输出 FASTA 文件路径
    """
    with input_path.open('r', encoding='utf-8') as src, \
         output_path.open('w', encoding='utf-8') as dst:
        for line in src:
            # 如果是描述行（以 '>' 开头），则插入 'gi|' 前缀
            if line.startswith('>'):
                dst.write(f">gi|{line[1:]}")
            else:
                dst.write(line)


def main():
    # 获取当前用户的桌面路径
    desktop = Path.home() / "Desktop"
    # 输入、输出文件路径（请根据需要修改文件名）
    input_fasta  = desktop / "20K_A.fasta"
    output_fasta = desktop / "mtPhyl.fasta"

    # 执行前缀替换
    replace_fasta_prefix(input_fasta, output_fasta)

    print(f"替换完成，新的 FASTA 已保存：{output_fasta}")


if __name__ == "__main__":
    main()
