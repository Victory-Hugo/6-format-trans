BASE_DIR='/mnt/f/0_现代DNA处理流程/script/2_ML_MRCA/script'

/home/luolintao/miniconda3/bin/python3 \
    ${BASE_DIR}/3_结果整理.py \
    ${BASE_DIR}/../OUTPUT/Result \
    ${BASE_DIR}/../OUTPUT/

/home/luolintao/miniconda3/bin/python3 \
    ${BASE_DIR}/4_分子钟校准.py \
    -t ${BASE_DIR}/../OUTPUT/ID_Length.tree \
    -c ${BASE_DIR}/../OUTPUT/TIP_Length.csv \
    -m 2.53e-8 \
    -o ${BASE_DIR}/../Revised \

