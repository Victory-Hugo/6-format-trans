cd /mnt/f/0_现代DNA处理流程/script/2_ML_MRCA/script

# /mnt/e/Scientifc_software/iqtree-3.0.0-Linux-intel/bin/iqtree3 \
#   -redo \
#   -m HKY85+G \
#   -s ../data/Example.fasta \
#   -nt 16 \
#   -o GQ119016 \
#   -pre ../data//Example

/home/luolintao/S07_iqtree3/iqtree3 \
  -redo \
  -m "MIX{TVM+FO}" \
  -s ../data/Example.fasta \
  -nt 16 \
  -o HGDP00453_B2b1b1b,HGDP00479_B2b1b1b,HGDP00478_B2a1a1a1,HGDP00462_B2a1a1a1 \
  -pre ../data/Example