

#install.packages("corrplot")


#引用包
library(corrplot)

riskFile="risk.TCGAall.txt"             #风险文件
immuneFile="MCPcounter.result.txt"      #免疫细胞浸润结果
setwd("C:\\biowolf\\NMF\\35.riskImmCor")     #设置工作目录

#读取免疫细胞浸润结果文件
data=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#相关性矩阵
M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)

#绘制相关性图形
pdf(file="immuneCor.pdf", width=7, height=7)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()




