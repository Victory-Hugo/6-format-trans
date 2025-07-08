# # BEAST XML 序列补全脚本
# ## 简介
# 本脚本用于将一个 FASTA 文件中缺失的序列，批量添加到 BEAST 输出的 XML 文件的 `<data>` 节点中。  
# - 自动检测并避免与 XML 中已有的 `sequence id` 和 `taxon` 重复。  
# - 若发现重复 ID，会一次性列出所有冲突并退出，不产生新的 XML。  
# - 若所有 FASTA 序列均已存在，直接提示无需添加。
# ## 环境依赖
# - Python 3.x（已在 shebang 中指定解释器路径，也可手动替换为系统中可用的 Python3）  
# - 标准库：`re`、`argparse`、`xml.etree.ElementTree`  
# - 无需第三方模块
# ## 文件结构



/home/luolintao/miniconda3/envs/pyg/bin/python3 \
    /mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/1-BEAST/src/1-xml→add-seq.py \
    --xml-in /mnt/c/Users/Administrator/Desktop/tree_100.xml \
    --fasta-in /mnt/c/Users/Administrator/Desktop/全部_no100.fasta \
    --xml-out /mnt/c/Users/Administrator/Desktop/tree_all.xml
