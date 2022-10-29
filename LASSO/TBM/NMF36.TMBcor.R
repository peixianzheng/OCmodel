

#install.packages("corrplot")
#install.packages("circlize")


#引用包
library(corrplot)
library(circlize)

tmbFile="TMB.txt"                       #肿瘤突变文件
riskFile="risk.TCGAall.txt"             #风险文件
immuneFile="MCPcounter.result.txt"      #免疫细胞浸润结果
setwd("C:\\biowolf\\NMF\\36.TMBcor")     #设置工作目录

#读取免疫细胞浸润结果文件
data=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#合并肿瘤突变符合文件
TMB=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(TMB), row.names(data))
rt=cbind(TMB[sameSample,,drop=F], data[sameSample,,drop=F])

#计算相关相关性矩阵
cor1=cor(rt)

#设置图形颜色
col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

#绘制圈图
pdf(file="circos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))       #绘制图例
dev.off()
circos.clear()



