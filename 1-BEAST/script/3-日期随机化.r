# 载入TipDatingBeast库，如果没有安装则先安装
if (!require("TipDatingBeast")) {
  install.packages("TipDatingBeast", dependencies = TRUE)
  library(TipDatingBeast)
}

# 指定工作路径
setwd("C:/Users/victo/Desktop")

# 指定你要产生多少份xml文件用于检测，默认20个
rep = 20
# 指定你要烧掉前%多少的树，默认0.1即10%
burnin = 0.1 

# 执行 RandomDates 函数，生成随机化日期的 XML 文件
RandomDates(name="TIP", rep, writeTrees=FALSE)
# 使用 PlotDRT 函数分析和绘制结果
PlotDRT(name="TIP", rep, burnin=0.1)
