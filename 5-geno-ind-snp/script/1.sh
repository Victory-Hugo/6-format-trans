plink --vcf /mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544.NoN.SNP.maf99.vcf.recode.vcf.gz \
      --allow-extra-chr \          # 支持非人类染色体
      --set-missing-var-ids @:# \  # 自动生成SNP ID（格式：chr:pos）
      --geno 0.1 \                 # 剔除缺失率>10%的位点
      --maf 0.05 \                 # 保留次等位基因频率≥5%的位点
      --mind 0.1 \                 # 剔除缺失率>10%的样本
      --make-bed \                 
      --out filtered

plink --bfile filtered \
      --indep-pairwise 50 5 0.2 \  # 窗口50kb，步长5，r²>0.2剪裁
      --out pruned
plink --bfile filtered \
      --extract pruned.prune.in \   # 保留剪裁后位点
      --make-bed \
      --out pruned_data

convertf -p parfile.txt

#todo 默认生成的.ind文件最后一列为占位符（?），需手动替换为真实群体标签：
# 原始示例
Sample1 M ?
# 修正后（替换"?"为群体名）
Sample1 M POP1
#todo 文件一致性检查​​
#? .geno行数 = .snp文件行数（位点数）
#? .geno列数 = .ind文件行数（样本数）
