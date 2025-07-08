# VCF文件转换为geno, ind, snp格式 - 高效版本
# 用于处理大型VCF文件，使用分块读取避免内存溢出

# 检查并安装必要的包
packages <- c("data.table")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 输入和输出路径
input_vcf <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544.NoN.SNP.maf99.vcf.recode.vcf.gz"
output_dir <- "/mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/5-f3/script/"

# 检查输入文件是否存在
if (!file.exists(input_vcf)) {
  stop("输入VCF文件不存在: ", input_vcf)
}

cat("开始处理VCF文件:", input_vcf, "\n")

# 使用系统命令读取VCF头部信息
cat("读取VCF头部信息...\n")
header_cmd <- paste0("zcat '", input_vcf, "' | grep '^#CHROM' | head -1")
header_line <- system(header_cmd, intern = TRUE)
headers <- strsplit(header_line, "\t")[[1]]

# 提取样本名（从第10列开始）
sample_names <- headers[10:length(headers)]
n_samples <- length(sample_names)

cat("样本数:", n_samples, "\n")

# 获取SNP总数
snp_count_cmd <- paste0("zcat '", input_vcf, "' | grep -v '^#' | wc -l")
n_snps <- as.numeric(system(snp_count_cmd, intern = TRUE))

cat("SNP数:", n_snps, "\n")

# 创建输出文件
geno_file <- file.path(output_dir, "hp_7544.geno")
ind_file <- file.path(output_dir, "hp_7544.ind")
snp_file <- file.path(output_dir, "hp_7544.snp")

# 创建ind文件
cat("创建ind文件...\n")
ind_data <- data.frame(
  sample = sample_names,
  sex = "U",  # U表示未知性别
  population = sample_names  # 可以根据需要修改群体标签
)
write.table(ind_data, file = ind_file, 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 处理VCF文件并创建geno和snp文件
cat("开始分块处理VCF数据...\n")

# 打开输出文件连接
geno_conn <- file(geno_file, "w")
snp_conn <- file(snp_file, "w")

# 处理函数
process_vcf_line <- function(line) {
  fields <- strsplit(line, "\t")[[1]]
  
  # 提取SNP信息
  chrom <- fields[1]
  pos <- as.numeric(fields[2])
  ref_allele <- fields[4]
  alt_allele <- fields[5]
  
  # 提取基因型（从第10列开始）
  genotypes <- fields[10:length(fields)]
  
  # 转换基因型编码
  geno_codes <- character(length(genotypes))
  for (i in 1:length(genotypes)) {
    gt <- strsplit(genotypes[i], ":")[[1]][1]  # 提取GT字段
    if (gt == "0/0" || gt == "0|0") {
      geno_codes[i] <- "2"
    } else if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
      geno_codes[i] <- "1"
    } else if (gt == "1/1" || gt == "1|1") {
      geno_codes[i] <- "0"
    } else {
      geno_codes[i] <- "9"  # 缺失或其他
    }
  }
  
  return(list(
    geno_line = paste(geno_codes, collapse = ""),
    snp_info = c(paste0("snp_", pos), chrom, "0", pos, ref_allele, alt_allele)
  ))
}

# 分块读取VCF文件
chunk_size <- 1000  # 每次处理1000行
processed_snps <- 0

# 使用系统命令逐行处理
vcf_cmd <- paste0("zcat '", input_vcf, "' | grep -v '^#'")
con <- pipe(vcf_cmd, "r")

cat("处理进度: ")
while (TRUE) {
  # 读取一批行
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0) break
  
  # 处理每一行
  for (line in lines) {
    if (nchar(line) > 0) {
      result <- process_vcf_line(line)
      
      # 写入geno文件
      writeLines(result$geno_line, geno_conn)
      
      # 写入snp文件
      writeLines(paste(result$snp_info, collapse = "\t"), snp_conn)
      
      processed_snps <- processed_snps + 1
    }
  }
  
  # 显示进度
  if (processed_snps %% 10000 == 0) {
    cat(processed_snps, "")
  }
}

close(con)
close(geno_conn)
close(snp_conn)

cat("\n处理完成!\n")

# 验证结果
cat("\n验证输出文件:\n")
if (file.exists(geno_file)) {
  geno_lines <- length(readLines(geno_file, warn = FALSE))
  cat("geno文件行数:", geno_lines, "\n")
}

if (file.exists(snp_file)) {
  snp_lines <- length(readLines(snp_file, warn = FALSE))
  cat("snp文件行数:", snp_lines, "\n")
}

if (file.exists(ind_file)) {
  ind_lines <- length(readLines(ind_file, warn = FALSE))
  cat("ind文件行数:", ind_lines, "\n")
}

# 检查文件大小
file_sizes <- sapply(c(geno_file, ind_file, snp_file), function(x) {
  if (file.exists(x)) {
    round(file.size(x) / 1024^2, 2)  # MB
  } else {
    0
  }
})
names(file_sizes) <- c("geno (MB)", "ind (MB)", "snp (MB)")
print(file_sizes)

cat("\n输出文件:\n")
cat("- geno文件:", geno_file, "\n")
cat("- ind文件:", ind_file, "\n") 
cat("- snp文件:", snp_file, "\n")

cat("\n注意事项:\n")
cat("1. geno文件编码: 0=纯合子替换, 1=杂合子, 2=纯合子参考, 9=缺失\n")
cat("2. ind文件中所有样本性别设为'U'(未知)，群体标签为样本名\n")
cat("3. snp文件中遗传距离设为0，可根据需要修改\n")
cat("4. 如需修改群体标签，请编辑ind文件的第三列\n")
