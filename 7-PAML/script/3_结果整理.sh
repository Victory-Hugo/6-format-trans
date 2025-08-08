BASE_DIR='/mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/7-PAML/script'

/home/luolintao/miniconda3/bin/python3 \
    ${BASE_DIR}/3_结果整理.py \
    /mnt/f/6_起源地混合地/2-系统发育树/1-ML/output/PAML/Result_Africa \
    /mnt/f/6_起源地混合地/2-系统发育树/1-ML/output/PAML/

/home/luolintao/miniconda3/bin/python3 \
    ${BASE_DIR}/4_分子钟校准.py \
    -t /mnt/f/6_起源地混合地/2-系统发育树/1-ML/output/PAML/ID_Length.tree \
    -c /mnt/f/6_起源地混合地/2-系统发育树/1-ML/output/PAML/TIP_Length.csv \
    -m 2.53e-8 \
    -o /mnt/f/6_起源地混合地/2-系统发育树/1-ML/output/PAML/Revised \

