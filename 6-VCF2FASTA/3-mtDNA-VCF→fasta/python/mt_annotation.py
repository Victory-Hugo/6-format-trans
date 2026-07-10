"""人类线粒体基因组 (rCRS / NC_012920.1) 的特征注释与分区构建。

坐标一律为 1-based、闭区间，与 GenBank 标注一致。

分区归属的两条约定（由项目决策确定，见 build_codon_map 的实现）:

1. ND6 位于轻链 (complement)，其密码子第 1 位从基因 3' 端 (14673) 起向回计数。
2. 重叠基因的重叠碱基归属下游基因的读码框：
   ATP8/ATP6 重叠 8527-8572 归 ATP6；ATP6/COX3 重叠 9207 归 COX3；
   ND4L/ND4 重叠 10760-10766 归 ND4。
"""

from __future__ import annotations

MT_LENGTH = 16569

# rCRS 中 3107 是历史遗留的占位符 N，不对应任何真实碱基。
PLACEHOLDER_POSITION = 3107

# (name, start, end, on_light_strand)
CDS_FEATURES: tuple[tuple[str, int, int, bool], ...] = (
    ("ND1", 3307, 4262, False),
    ("ND2", 4470, 5511, False),
    ("COX1", 5904, 7445, False),
    ("COX2", 7586, 8269, False),
    ("ATP8", 8366, 8572, False),
    ("ATP6", 8527, 9207, False),
    ("COX3", 9207, 9990, False),
    ("ND3", 10059, 10404, False),
    ("ND4L", 10470, 10766, False),
    ("ND4", 10760, 12137, False),
    ("ND5", 12337, 14148, False),
    ("ND6", 14149, 14673, True),
    ("CYTB", 14747, 15887, False),
)

RRNA_FEATURES: tuple[tuple[str, int, int], ...] = (
    ("RNR1_12S", 648, 1601),
    ("RNR2_16S", 1671, 3229),
)

TRNA_FEATURES: tuple[tuple[str, int, int], ...] = (
    ("tRNA-Phe", 577, 647),
    ("tRNA-Val", 1602, 1670),
    ("tRNA-Leu_UUR", 3230, 3304),
    ("tRNA-Ile", 4263, 4331),
    ("tRNA-Gln", 4329, 4400),
    ("tRNA-Met", 4402, 4469),
    ("tRNA-Trp", 5512, 5579),
    ("tRNA-Ala", 5587, 5655),
    ("tRNA-Asn", 5657, 5729),
    ("tRNA-Cys", 5761, 5826),
    ("tRNA-Tyr", 5826, 5891),
    ("tRNA-Ser_UCN", 7446, 7514),
    ("tRNA-Asp", 7518, 7585),
    ("tRNA-Lys", 8295, 8364),
    ("tRNA-Gly", 9991, 10058),
    ("tRNA-Arg", 10405, 10469),
    ("tRNA-His", 12138, 12206),
    ("tRNA-Ser_AGY", 12207, 12265),
    ("tRNA-Leu_CUN", 12266, 12336),
    ("tRNA-Glu", 14674, 14742),
    ("tRNA-Thr", 15888, 15953),
    ("tRNA-Pro", 15956, 16023),
)

# 控制区跨越环状基因组的原点，在线性 rCRS 坐标上被劈成两段。
CONTROL_REGION_SEGMENTS: tuple[tuple[int, int], ...] = ((16024, 16569), (1, 576))

# 项目决策：这些 rCRS 位置整体从比对中剔除（不论其上的记录是 indel 还是 SNV）。
# 依据是它们落在 poly-C / AC-repeat 区，比对本身不可靠。
DEFAULT_MASK_POSITIONS: tuple[int, ...] = (
    309,  # 309.1C 插入所在位置
    315,  # 315.1C 插入所在位置
    515, 516, 517, 518, 519, 520, 521, 522,  # 515-522 AC indel 区
    16182,  # A16182c
    16183,  # A16183c
    16193,  # 16193.1C 插入所在位置
    16519,  # C16519T / T16519C
)


def build_codon_map() -> dict[int, int]:
    """返回 rCRS 位置 -> 密码子位置 (1/2/3) 的映射。

    按基因 start 升序赋值，后写覆盖先写，从而让重叠碱基自动归属下游基因。
    """
    codon_of: dict[int, int] = {}
    for _name, start, end, light in sorted(CDS_FEATURES, key=lambda f: f[1]):
        if light:
            # 轻链基因：读码框从 3' 端 (end) 向 5' 端 (start) 推进。
            for offset, pos in enumerate(range(end, start - 1, -1)):
                codon_of[pos] = offset % 3 + 1
        else:
            for offset, pos in enumerate(range(start, end + 1)):
                codon_of[pos] = offset % 3 + 1
    return codon_of


def _expand(segments) -> set[int]:
    return {p for start, end in segments for p in range(start, end + 1)}


def build_partition_map(
    masked: frozenset[int], *, catchall: bool = True
) -> dict[str, list[int]]:
    """返回 分区名 -> 该分区包含的 rCRS 位置列表（升序、已剔除 masked）。

    优先级由低到高：control_region < rRNA < tRNA < CDS 密码子位。
    CDS 最后写入，因此与结构 RNA 若有重叠一律按 CDS 归类。

    catchall 为 True 时，把不属于任何已知特征的位点（基因间隔区等）收进
    ``noncoding_other``，避免它们在 IQ-TREE ``-p`` 下被静默丢弃。
    """
    owner: dict[int, str] = {}

    for pos in _expand(CONTROL_REGION_SEGMENTS):
        owner[pos] = "control_region"
    for _name, start, end in RRNA_FEATURES:
        for pos in range(start, end + 1):
            owner[pos] = "rRNA"
    for _name, start, end in TRNA_FEATURES:
        for pos in range(start, end + 1):
            owner[pos] = "tRNA"
    for pos, codon in build_codon_map().items():
        owner[pos] = f"cds_pos{codon}"

    if catchall:
        for pos in range(1, MT_LENGTH + 1):
            owner.setdefault(pos, "noncoding_other")

    names = ["control_region", "rRNA", "tRNA", "cds_pos1", "cds_pos2", "cds_pos3"]
    if catchall:
        names.append("noncoding_other")

    partitions: dict[str, list[int]] = {name: [] for name in names}
    for pos in range(1, MT_LENGTH + 1):
        if pos in masked:
            continue
        name = owner.get(pos)
        if name is not None:
            partitions[name].append(pos)
    return {name: positions for name, positions in partitions.items() if positions}
