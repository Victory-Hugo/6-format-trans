# 在开始之前请确保你已经安装了下列package
# 载入需要的package
# install.packages("ape")
# install.packages("igraph")
# install.packages("ggplot2")
# install.packages("pheatmap")
# install.packages("reshape2")
# install.packages("ggsci")
# install.packages("gridExtra")

library(ape)
library(igraph)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggsci)
library(gridExtra)
# library(httpgd)
# 设置工作路径，你需要把Fst文件和分组文件放到工作路径
setwd("/mnt/c/Users/Administrator/Desktop")

# 清除内存
rm(list=ls())

# 读取数据
mydata<-read.table("/mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/2-Alequin/conf/Fst.csv",header=TRUE,sep=",",row.names = 1) # 输入Fst文件
group <-read.table("/mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/2-Alequin/conf/group.csv",header=TRUE,sep=",", row.names = 1) # 输入标签文件

# 创建绘图PDF
#pdf("FstMatrix1.pdf", width=11, height=8.5)
dev.off()
pdf("Matrix1.pdf",
      width=11, height=8.5)
# ! 下面的样式任选其一进行绘制！
# 绘制图像
pheatmap(mydata, 
         cluster_cols=TRUE, 
         cluster_rows=TRUE, 
         angle_col = c("45"),
         fontsize = 8,
         fontsize_row =8,
         fontsize_col =6,
         annotation_col = group, 
         annotation_row = group,
         cellwidth =8, 
         cellheight = 8, 
         cutree_cols=4, # 这个数字决定最终的矩形被分成几个部分
         cutree_rows=4, # 这个数字决定最终的矩形被分成几个部分
         main = "FstMatrix",
         color = colorRampPalette(c("#00A087FF","#3C5488FF","#F39B7FFF"))(10000),
         display_numbers = matrix(ifelse(abs(mydata)> 50, "++", ifelse(abs(mydata)>=40,"+"," ")), nrow(mydata)))
# 结束绘图
dev.off()

# 换个样式
pheatmap(mydata, 
         cluster_cols=TRUE, 
         cluster_rows=TRUE, 
         angle_col = c("45"),
         fontsize = 8,
         fontsize_row =8,
         fontsize_col =6,
         annotation_col = group, 
         annotation_row = group,
         cellwidth =8, 
         cellheight = 8, 
         cutree_cols=4,
         cutree_rows=4, 
         main = "FstMatrix",
         color = colorRampPalette(c("#F8F8FF","#91D1C2FF","#3C5488FF"))(10000),
         display_numbers = matrix(ifelse(abs(mydata)> 50, "++", ifelse(abs(mydata)>=40,"+"," ")), nrow(mydata)))

# 换个样式
pheatmap(mydata,
         cluster_cols=TRUE, # 是否进行列聚类
         cluster_rows=TRUE, # 是否进行行聚类
         angle_col = c("45"), # x轴文字角度
         fontsize = 1, # 字体大小
         fontsize_row = 1, # 横轴文字大小
         fontsize_col = 1, # 纵轴字体大小
         annotation_col = group, # 列是否进行标签的标记
         annotation_row = group, # 行是否进行标签的标记
         cellwidth = 1, # 单元格的宽度
         cellheight = 1, # 单元格的高度
         cutree_cols = 10, # 手动调整列聚类的数量 # ! 聚类的数量并不代表群体真实的分层，受限于算法的局限性。
         cutree_rows = 10, # 手动调整行聚类的数量 # ! 并不一定聚类到一起的群体都是相似的群体!
         main = "Fst matrix",
         color = colorRampPalette(c("#20364F","#31646C","#4E9280","#96B89B","#DCDFD2","#ECD9CF","#D49C87","#B86265","#8B345E","#50184E"))(100), # !选择你喜欢的颜色列表！
         #color = colorRampPalette(c("#023047","#126883","#279EBC","#90C9E6","#FC9E7F","#F75B41","#D52120"))(100),
         #color = colorRampPalette(c("#50184E","#B86265","#D49C87","#DCDFD2","#96B89B","#4E9280","#31646C","#003936"))(10000),
         #display_numbers = matrix(ifelse(abs(mydata) < 0.1, "**", ifelse(abs(mydata) < 0.05, "***", ifelse(abs(mydata) > 90, "++", ifelse(abs(mydata) >= 80, "+", " ")))), nrow(mydata)), #! 是否在单元格中展示具体的数字或者符号
         #border_color = "#FFFFFF",  # 设置边框颜色为白色
         lwd = 0.10)  # 设置边框线宽为0.25

dev.off()

# 换个样式
pheatmap(mydata,
         cluster_cols=TRUE,
         cluster_rows=TRUE,
         angle_col = c("45"),
         fontsize = 8,
         fontsize_row =8,
         fontsize_col =6,
         annotation_col = group,
         annotation_row = group,
         cellwidth =8,
         cellheight = 8,
         cutree_cols=4,
         cutree_rows=4,
         main = "FstMatrix",
         color = colorRampPalette(c("#023047","#126883","#279EBC","#90C9E6","#FC9E7F","#F75B41","#D52120","#20364F","#31646C","#4E9280","#96B89B",))(10000),
         display_numbers = matrix(ifelse(abs(mydata)> 50, "++", ifelse(abs(mydata)>=40,"+"," ")), nrow(mydata)))
