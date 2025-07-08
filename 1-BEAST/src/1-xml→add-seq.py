#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re, argparse

def parse_fasta(path):
    """生成器：边读边产出 (seq_id, seq_str)，不占用额外内存。"""
    with open(path, 'r') as f:
        header = None
        buf = []
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line[0] == '>':
                if header:
                    yield header, ''.join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if header:
            yield header, ''.join(buf)

def main():
    p = argparse.ArgumentParser(
        description="高性能向 BEAST XML 批量添加 FASTA 序列")
    p.add_argument('-x','--xml-in', required=True)
    p.add_argument('-f','--fasta-in', required=True)
    p.add_argument('-o','--xml-out', required=True)
    args = p.parse_args()

    # 1. 预编译正则
    spec_re = re.compile(r'\bsequence\b[^>]*\bspec="([^"]+)"')
    id_re   = re.compile(r'\bid="([^"]+)"')
    tax_re  = re.compile(r'\btaxon="([^"]+)"')

    existing_ids = set()
    existing_taxons = set()
    spec = None

    # 2. 一次扫描——既找 spec，又收集已有 ID/taxon
    with open(args.xml_in, 'r', encoding='utf-8') as fin:
        for line in fin:
            if '<sequence' in line:
                if spec is None:
                    m = spec_re.search(line)
                    if m:
                        spec = m.group(1)
                m1 = id_re.search(line)
                if m1:
                    existing_ids.add(m1.group(1))
                m2 = tax_re.search(line)
                if m2:
                    existing_taxons.add(m2.group(1))
            if spec and existing_ids and existing_taxons:
                # 找到一个完整的 sequence 块后，就可以停掉对 spec 的继续查找
                pass

    if not spec:
        print("无法提取 spec，退出。"); return

    # 3. 快速收集：哪些 FASTA ID 重复、哪些是真新序列
    duplicates = []
    new_xml_snippets = []
    for seq_id, seq_str in parse_fasta(args.fasta_in):
        if seq_id in existing_ids or seq_id in existing_taxons:
            duplicates.append(seq_id)
        else:
            new_xml_snippets.append(
                f'    <sequence id="seq_{seq_id}" spec="{spec}" '
                f'taxon="{seq_id}" totalcount="4" value="{seq_str}" />\n'
            )

    if duplicates:
        print("发现重复的 ID，请删除这些序列：")
        for d in duplicates:
            print("", d)
        return
    if not new_xml_snippets:
        print("无需添加新序列，所有已存在。")
        return

    # 4. 流式写入：读原文件、写新文件，遇到 </data> 先插入
    with open(args.xml_in, 'r', encoding='utf-8') as fin, \
         open(args.xml_out,'w',encoding='utf-8') as fout:
        for line in fin:
            if line.strip() == '</data>':
                fout.writelines(new_xml_snippets)
            fout.write(line)

    print(f"完成！新增 {len(new_xml_snippets)} 条 <sequence>，输出至 {args.xml_out}")

if __name__ == '__main__':
    main()
