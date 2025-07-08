# VCF文件转换为geno, ind, snp格式
# 用于ADMIXTOOLS/f3统计分析

# 检查并安装必要的包
if (!require(vcfR, quietly = TRUE)) {
  cat("正在安装vcfR包...\n")
  if (!require(BiocManager, quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("vcfR")
}

if (!require(data.table, quietly = TRUE)) {
  install.packages("data.table")
}

# 加载必要的包
library(vcfR)
library(data.table)

# 输入和输出路径
input_vcf <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544.NoN.SNP.maf99.vcf.recode.vcf.gz"
output_dir <- "/mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/5-f3/script/"

# 检查输入文件是否存在
if (!file.exists(input_vcf)) {
  stop("输入VCF文件不存在: ", input_vcf)
}

# 读取VCF文件
cat("正在读取VCF文件...\n")
vcf <- read.vcfR(input_vcf)

# 获取基本信息
cat("VCF文件包含:\n")
cat("样本数:", ncol(vcf@gt) - 1, "\n")
cat("SNP数:", nrow(vcf@gt), "\n")

# 提取基因型数据
gt_matrix <- extract.gt(vcf, element = "GT")

# 转换基因型编码为数字格式
# 0/0 -> 2 (纯合子参考)
# 0/1 或 1/0 -> 1 (杂合子)  
# 1/1 -> 0 (纯合子替换)
# ./. 或其他 -> 9 (缺失)

cat("正在转换基因型编码...\n")
geno_matrix <- gt_matrix
geno_matrix[geno_matrix == "0/0" | geno_matrix == "0|0"] <- "2"
geno_matrix[geno_matrix == "0/1" | geno_matrix == "1/0" | geno_matrix == "0|1" | geno_matrix == "1|0"] <- "1"
geno_matrix[geno_matrix == "1/1" | geno_matrix == "1|1"] <- "0"
geno_matrix[is.na(geno_matrix) | geno_matrix == "./." | geno_matrix == ".|."] <- "9"

# 确保所有值都是数字
geno_matrix[!geno_matrix %in% c("0", "1", "2", "9")] <- "9"

# 1. 创建geno文件
cat("正在创建geno文件...\n")
geno_file <- file.path(output_dir, "hp_7544.geno")

# geno文件格式：每行是一个SNP，每列是一个个体
# 需要转置矩阵（因为我们要按SNP写行）
geno_t <- t(geno_matrix)
write.table(geno_t, file = geno_file, 
            sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 2. 创建ind文件
cat("正在创建ind文件...\n")
ind_file <- file.path(output_dir, "hp_7544.ind")

# ind文件格式：样本名 性别 群体
# 从VCF文件获取样本名
sample_names <- colnames(gt_matrix)

# 创建ind数据框
ind_data <- data.frame(
  sample = sample_names,
  sex = "U",  # U表示未知性别
  population = sample_names  # 可以根据需要修改群体标签
)

write.table(ind_data, file = ind_file, 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 3. 创建snp文件
cat("正在创建snp文件...\n")
snp_file <- file.path(output_dir, "hp_7544.snp")

# snp文件格式：SNP名 染色体 遗传距离(cM) 物理位置(bp) 参考等位基因 替换等位基因
# 从VCF文件获取SNP信息
snp_info <- data.frame(
  snp_id = paste0("snp_", 1:nrow(vcf@fix)),  # SNP ID
  chr = vcf@fix[, "CHROM"],                   # 染色体
  genetic_dist = 0,                          # 遗传距离（设为0）
  position = as.numeric(vcf@fix[, "POS"]),   # 物理位置
  ref_allele = vcf@fix[, "REF"],             # 参考等位基因
  alt_allele = vcf@fix[, "ALT"]              # 替换等位基因
)

write.table(snp_info, file = snp_file, 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 输出统计信息
cat("\n转换完成!\n")
cat("输出文件:\n")
cat("- geno文件:", geno_file, "\n")
cat("- ind文件:", ind_file, "\n") 
cat("- snp文件:", snp_file, "\n")
cat("\n文件统计:\n")
cat("- 样本数:", nrow(ind_data), "\n")
cat("- SNP数:", nrow(snp_info), "\n")

# 检查文件大小
file_sizes <- sapply(c(geno_file, ind_file, snp_file), function(x) {
  round(file.size(x) / 1024^2, 2)  # MB
})
names(file_sizes) <- c("geno (MB)", "ind (MB)", "snp (MB)")
print(file_sizes)

cat("\n注意事项:\n")
cat("1. geno文件编码: 0=纯合子替换, 1=杂合子, 2=纯合子参考, 9=缺失\n")
cat("2. ind文件中所有样本性别设为'U'(未知)，群体标签为样本名\n")
cat("3. snp文件中遗传距离设为0，可根据需要修改\n")
cat("4. 如需修改群体标签，请编辑ind文件的第三列\n")