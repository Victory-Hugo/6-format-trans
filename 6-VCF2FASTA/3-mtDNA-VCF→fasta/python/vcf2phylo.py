#!/usr/bin/env python3
"""从 mtDNA VCF 生成可直接用于 IQ-TREE 3 与 BEAST X 的比对与分区文件。

设计要点
--------
坐标系
    输出比对由 rCRS (NC_012920.1) 全长 16569 列删去若干 mask 列得到。删列后
    比对列号与 rCRS 位置不再相等，二者的对应关系写入 ``col2pos.tsv``，所有
    分区文件均由该映射生成，不存在第二处硬编码坐标。

只保留 SNV
    VCF 中的 indel 记录一律丢弃。这不仅是方法学上的取舍（GTR 不建模 indel，
    gap 被当作缺失数据，几乎不贡献似然），在本数据上还是数据质量上的必需：
    2159 条 indel 记录里有 652 条的 REF 与 rCRS 不符，而 13053 条纯 SNV 记录
    的 REF 与 rCRS 完全一致。

缺失与歧义码
    ``./.``、IUPAC 歧义码 (R/M/Y 等) 以及非标准的 ``D`` 一律写作 N。缺失不
    回填 rCRS 参考碱基——那会人为制造恒定位点，扭曲速率异质性的估计。

基因型解析
    依赖 bcftools ``%TGT`` 对 SNV 行严格输出 3 字节（``C/C`` / ``./.``）这一
    性质，用 numpy 按 4 字节 stride 切片，避免逐样本 split。
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np

from mt_annotation import (
    DEFAULT_MASK_POSITIONS,
    MT_LENGTH,
    PLACEHOLDER_POSITION,
    build_partition_map,
)

LOG = logging.getLogger("vcf2phylo")

ACGT = np.array([ord(c) for c in "ACGT"], dtype=np.uint8)
N_BYTE = ord("N")

# NEXUS 中必须以单引号包裹的字符（含空格）。
NEXUS_UNSAFE = set(" \t\n()[]{}/\\,;:=*'\"`<>-+")


class Stats:
    """转换过程中的计数器，最终写入 summary.txt。"""

    def __init__(self) -> None:
        self.records_total = 0
        self.records_indel = 0
        self.records_snv = 0
        self.records_snv_masked = 0
        self.records_snv_kept = 0
        self.ref_mismatch: list[tuple[int, str, str]] = []
        self.gt_missing = 0
        self.gt_heterozygous = 0
        self.gt_ambiguous = 0
        self.gt_total = 0
        # 在实际输出的比对上实测（而非从全量 VCF 的位点数推断）
        self.constant_sites = 0
        self.variable_sites = 0


def read_reference(path: Path) -> np.ndarray:
    """读取 rCRS，返回长度 16569 的 uint8 数组。"""
    seq = "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    ).upper()
    if len(seq) != MT_LENGTH:
        raise ValueError(f"参考序列长度为 {len(seq)}，期望 {MT_LENGTH}: {path}")
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()


def is_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and all(len(a) == 1 for a in alt.split(","))


def bcftools(args: list[str], *, stream: bool = False):
    cmd = ["bcftools", *args]
    LOG.debug("运行: %s", " ".join(cmd))
    if stream:
        return subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 22)
    done = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return done.stdout


def list_samples(vcf: Path, subset: Path | None) -> list[str]:
    names = bcftools(["query", "-l", str(vcf)]).split()
    if subset is None:
        return names
    wanted = [n for n in subset.read_text().split() if n]
    known = set(names)
    missing = [n for n in wanted if n not in known]
    if missing:
        raise SystemExit(
            f"--samples 中有 {len(missing)} 个名字不在 VCF 里，例如: {missing[:5]}"
        )
    return wanted


def scan_variant_sites(
    vcf: Path, masked: frozenset[int], reference: np.ndarray, stats: Stats
) -> list[int]:
    """轻量地过一遍 VCF（不读基因型），返回保留下来的 SNV 位置（升序）。"""
    text = bcftools(["query", "-f", "%POS\t%REF\t%ALT\n", str(vcf)])
    kept: list[int] = []
    for line in text.splitlines():
        pos_s, ref, alt = line.split("\t")
        pos = int(pos_s)
        stats.records_total += 1
        if not is_snv(ref, alt):
            stats.records_indel += 1
            continue
        stats.records_snv += 1
        observed = chr(reference[pos - 1])
        if observed != ref:
            stats.ref_mismatch.append((pos, ref, observed))
        if pos in masked:
            stats.records_snv_masked += 1
            continue
        kept.append(pos)
    stats.records_snv_kept = len(kept)
    if stats.ref_mismatch:
        raise SystemExit(
            f"有 {len(stats.ref_mismatch)} 条 SNV 记录的 REF 与参考序列不符，"
            f"例如 {stats.ref_mismatch[:5]}。这说明 VCF 与参考不是同一坐标系，已中止。"
        )
    return sorted(kept)


def load_genotypes(
    vcf: Path,
    subset: Path | None,
    kept_positions: list[int],
    n_samples: int,
    stats: Stats,
) -> np.ndarray:
    """返回 (n_samples, n_kept_sites) 的 uint8 碱基矩阵。"""
    index_of = {pos: i for i, pos in enumerate(kept_positions)}
    matrix = np.empty((n_samples, len(kept_positions)), dtype=np.uint8)

    args = ["query", "-f", "%POS\t%REF\t%ALT[\t%TGT]\n"]
    if subset is not None:
        args += ["-S", str(subset)]
    args.append(str(vcf))

    proc = bcftools(args, stream=True)
    assert proc.stdout is not None
    seen = 0
    for raw in proc.stdout:
        i1 = raw.index(b"\t")
        i2 = raw.index(b"\t", i1 + 1)
        i3 = raw.index(b"\t", i2 + 1)
        pos = int(raw[:i1])
        col = index_of.get(pos)
        if col is None:
            continue  # indel 行，或被 mask 的位置
        if not is_snv(raw[i1 + 1 : i2].decode(), raw[i2 + 1 : i3].decode()):
            continue  # 同一 POS 上的 indel 行

        gt = np.frombuffer(raw[i3 + 1 :].rstrip(b"\n"), dtype=np.uint8)
        if gt.size != 4 * n_samples - 1:
            raise SystemExit(
                f"pos {pos}: 基因型字段长度 {gt.size}，与 {n_samples} 个样本的定长"
                f"假设 (4N-1) 不符。可能存在非纯合或多字符等位。"
            )
        first, second = gt[0::4], gt[2::4]

        heterozygous = first != second
        acgt = np.isin(first, ACGT)
        base = np.where(acgt & ~heterozygous, first, N_BYTE).astype(np.uint8)

        stats.gt_total += n_samples
        stats.gt_missing += int(np.count_nonzero(first == ord(".")))
        stats.gt_heterozygous += int(np.count_nonzero(heterozygous))
        stats.gt_ambiguous += int(np.count_nonzero(~acgt & (first != ord("."))))

        matrix[:, col] = base
        seen += 1

    proc.stdout.close()
    if proc.wait() != 0:
        raise SystemExit("bcftools query 非零退出")
    if seen != len(kept_positions):
        raise SystemExit(f"只读到 {seen} 个位点，期望 {len(kept_positions)}")
    return matrix


def count_constant_sites(matrix: np.ndarray, n_cols: int, stats: Stats) -> None:
    """实测输出比对中的恒定位点数。

    ``matrix`` 只含 VCF 里的变异位点，但在样本子集下其中不少位点可能只剩一种
    等位。忽略 N 之后仍只有一种碱基的列算作恒定。未出现在 matrix 中的列全部
    等于 rCRS，必然恒定。
    """
    n_var = matrix.shape[1]
    variable = 0
    for start in range(0, n_var, 2048):  # 分块以限制临时布尔数组的峰值内存
        block = matrix[:, start : start + 2048]
        distinct = np.zeros(block.shape[1], dtype=np.uint8)
        for base in ACGT:
            distinct += (block == base).any(axis=0)
        variable += int(np.count_nonzero(distinct > 1))
    stats.variable_sites = variable
    stats.constant_sites = n_cols - variable


def compress_columns(cols: list[int]) -> str:
    """把升序列号压成 NEXUS charset 语法：``a-b``、``a-b\\3``、或单点。"""
    parts: list[str] = []
    i, n = 0, len(cols)
    while i < n:
        run_end = i
        step = None
        if i + 1 < n:
            step = cols[i + 1] - cols[i]
            j = i + 1
            while j + 1 < n and cols[j + 1] - cols[j] == step:
                j += 1
            run_end = j
        length = run_end - i + 1
        if step in (1, 3) and length >= 3:
            first, last = cols[i], cols[run_end]
            parts.append(f"{first}-{last}" if step == 1 else f"{first}-{last}\\3")
            i = run_end + 1
        else:
            parts.append(str(cols[i]))
            i += 1
    return " ".join(parts)


def nexus_quote(name: str) -> str:
    if any(c in NEXUS_UNSAFE for c in name):
        return "'" + name.replace("'", "''") + "'"
    return name


def write_col2pos(path: Path, kept_pos_of_col: list[int]) -> None:
    with path.open("w") as fh:
        fh.write("col\trcrs_pos\n")
        for col, pos in enumerate(kept_pos_of_col, start=1):
            fh.write(f"{col}\t{pos}\n")


def build_sets_block(partitions: dict[str, list[int]], col_of_pos: dict[int, int]) -> str:
    lines = ["begin sets;"]
    for name, positions in partitions.items():
        cols = sorted(col_of_pos[p] for p in positions)
        lines.append(f"  charset {name} = {compress_columns(cols)};")
    lines.append("end;")
    return "\n".join(lines)


def write_fasta(
    path: Path, names: list[str], template: np.ndarray, var_cols: np.ndarray, matrix: np.ndarray
) -> None:
    buf = np.empty_like(template)
    with path.open("wb") as fh:
        for i, name in enumerate(names):
            buf[:] = template
            buf[var_cols] = matrix[i]
            fh.write(b">" + name.encode() + b"\n" + buf.tobytes() + b"\n")


def write_nexus(
    path: Path,
    names: list[str],
    template: np.ndarray,
    var_cols: np.ndarray,
    matrix: np.ndarray,
    sets_block: str | None,
) -> None:
    quoted = [nexus_quote(n) for n in names]
    width = max(len(q) for q in quoted)
    buf = np.empty_like(template)
    with path.open("wb") as fh:
        fh.write(b"#NEXUS\n\nbegin data;\n")
        fh.write(f"  dimensions ntax={len(names)} nchar={template.size};\n".encode())
        fh.write(b"  format datatype=nucleotide missing=N gap=-;\n  matrix\n")
        for i, label in enumerate(quoted):
            buf[:] = template
            buf[var_cols] = matrix[i]
            fh.write(label.ljust(width).encode() + b" " + buf.tobytes() + b"\n")
        fh.write(b"  ;\nend;\n")
        if sets_block:
            fh.write(b"\n" + sets_block.encode() + b"\n")


def write_summary(
    path: Path,
    stats: Stats,
    masked: frozenset[int],
    n_cols: int,
    n_samples: int,
    partitions: dict[str, list[int]],
) -> None:
    pct = 100 * stats.gt_missing / stats.gt_total if stats.gt_total else 0.0
    lines = [
        "vcf2phylo 转换报告",
        "=" * 60,
        "",
        f"样本数            {n_samples}",
        f"比对列数          {n_cols}   (= {MT_LENGTH} - {len(masked)} 个 mask 位置)",
        "",
        "VCF 记录",
        f"  总记录          {stats.records_total}",
        f"  indel (丢弃)    {stats.records_indel}",
        f"  SNV             {stats.records_snv}",
        f"    被 mask       {stats.records_snv_masked}",
        f"    保留          {stats.records_snv_kept}",
        "",
        "输出比对（在本次样本集合上实测）",
        f"  可变位点        {stats.variable_sites}"
        f"   ({100 * stats.variable_sites / n_cols:.1f}%)",
        f"  恒定位点        {stats.constant_sites}"
        f"   ({100 * stats.constant_sites / n_cols:.1f}%)",
        "",
        "基因型 (仅统计保留的 SNV 位点)",
        f"  总计            {stats.gt_total}",
        f"  缺失 ./.  -> N  {stats.gt_missing}  ({pct:.3f}%)",
        f"  歧义码    -> N  {stats.gt_ambiguous}",
        f"  杂合      -> N  {stats.gt_heterozygous}",
        "",
        f"mask 的 rCRS 位置 ({len(masked)} 个)",
        f"  {sorted(masked)}",
        "",
        "分区 (列数)",
    ]
    for name, positions in partitions.items():
        lines.append(f"  {name:<18} {len(positions)}")
    lines.append(f"  {'合计':<18} {sum(len(p) for p in partitions.values())}")
    path.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="从 mtDNA VCF 生成 IQ-TREE 3 / BEAST X 的比对与分区文件",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-v", "--vcf", type=Path, required=True, help="输入 VCF (bgzip + index)")
    parser.add_argument("-r", "--ref", type=Path, required=True, help="rCRS FASTA (16569 bp)")
    parser.add_argument("-o", "--outdir", type=Path, required=True, help="输出目录")
    parser.add_argument(
        "-S", "--samples", type=Path, default=None, help="样本子集列表（每行一个名字）"
    )
    parser.add_argument(
        "--extra-mask",
        type=int,
        nargs="*",
        default=(),
        help="在默认 mask 之外额外剔除的 rCRS 位置",
    )
    parser.add_argument(
        "--no-catchall",
        action="store_true",
        help="不生成 noncoding_other 分区（未分配的列将被 IQ-TREE 静默丢弃）",
    )
    parser.add_argument(
        "--max-nexus-taxa",
        type=int,
        default=20000,
        help="超过该样本数则跳过 NEXUS 输出（BEAUti 载不动）",
    )
    parser.add_argument("--no-nexus", action="store_true", help="只产出 IQ-TREE 所需文件")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        datefmt="%H:%M:%S",
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    stats = Stats()

    masked = frozenset({*DEFAULT_MASK_POSITIONS, PLACEHOLDER_POSITION, *args.extra_mask})
    LOG.info("mask %d 个 rCRS 位置: %s", len(masked), sorted(masked))

    reference = read_reference(args.ref)
    samples = list_samples(args.vcf, args.samples)
    LOG.info("样本 %d 条", len(samples))

    want_nexus = not args.no_nexus and len(samples) <= args.max_nexus_taxa

    LOG.info("扫描变异位点 ...")
    kept_positions = scan_variant_sites(args.vcf, masked, reference, stats)
    LOG.info(
        "保留 %d 个 SNV 位点；丢弃 indel %d 条、被 mask 的 SNV %d 条",
        len(kept_positions),
        stats.records_indel,
        stats.records_snv_masked,
    )

    keep_mask = np.ones(MT_LENGTH, dtype=bool)
    keep_mask[[p - 1 for p in masked]] = False
    kept_pos_of_col = [p for p in range(1, MT_LENGTH + 1) if keep_mask[p - 1]]
    col_of_pos = {pos: col for col, pos in enumerate(kept_pos_of_col, start=1)}
    template = reference[keep_mask]
    n_cols = template.size
    LOG.info("比对列数 %d", n_cols)

    LOG.info("读取基因型 ...")
    matrix = load_genotypes(args.vcf, args.samples, kept_positions, len(samples), stats)
    LOG.info(
        "缺失 %d (%.3f%%)、歧义码 %d、杂合 %d，全部写作 N",
        stats.gt_missing,
        100 * stats.gt_missing / stats.gt_total,
        stats.gt_ambiguous,
        stats.gt_heterozygous,
    )

    count_constant_sites(matrix, n_cols, stats)
    LOG.info(
        "实测比对: 可变位点 %d (%.1f%%)、恒定位点 %d (%.1f%%)",
        stats.variable_sites,
        100 * stats.variable_sites / n_cols,
        stats.constant_sites,
        100 * stats.constant_sites / n_cols,
    )

    var_cols = np.array([col_of_pos[p] - 1 for p in kept_positions], dtype=np.int64)

    partitions = build_partition_map(masked, catchall=not args.no_catchall)
    sets_block = build_sets_block(partitions, col_of_pos)

    write_col2pos(args.outdir / "col2pos.tsv", kept_pos_of_col)
    LOG.info("写出 col2pos.tsv")

    fasta = args.outdir / "alignment.fasta"
    write_fasta(fasta, samples, template, var_cols, matrix)
    LOG.info("写出 %s", fasta.name)

    partition_nex = args.outdir / "partition.nex"
    partition_nex.write_text("#nexus\n\n" + sets_block + "\n")
    LOG.info("写出 partition.nex (IQ-TREE -p)")

    if want_nexus:
        write_nexus(
            args.outdir / "beast_partitioned.nex", samples, template, var_cols, matrix, sets_block
        )
        write_nexus(
            args.outdir / "beast_unpartitioned.nex", samples, template, var_cols, matrix, None
        )
        LOG.info("写出 beast_partitioned.nex / beast_unpartitioned.nex")
    elif args.no_nexus:
        LOG.info("--no-nexus，跳过 NEXUS 输出")
    else:
        LOG.warning(
            "样本数 %d 超过 --max-nexus-taxa=%d，跳过 NEXUS 输出。BEAST 请先用 -S 生成子集。",
            len(samples),
            args.max_nexus_taxa,
        )

    write_summary(args.outdir / "summary.txt", stats, masked, n_cols, len(samples), partitions)
    LOG.info("写出 summary.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
