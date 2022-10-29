

#install.packages("pheatmap")

library(pheatmap)
setwd("D:\\biowolf\\methySubtype\\15.cliHeatmap")                        #设置工作目录
rt=read.table("multiIndepSigExp.txt",sep="\t",header=T,row.names=1,check.names=F)      #读取甲基化文件
rt=t(rt[,c(-1,-2)])

#读取临床数据和分型文件
clinical=read.table("clinical.txt",sep="\t",header=T,row.names=1,check.names=F)
cluster=read.table("cluster.txt",sep="\t",header=T,row.names=1,check.names=F)
Type=clinical[row.names(cluster),]
Type=cbind(Type,Cluster=paste0("C",cluster[,"cluster"]))
#按照分型对样品排序
Type=Type[order(Type$Cluster),]
rt=rt[,row.names(Type)]
#绘制热图
pdf(file="heatmap.pdf",height=6.5,width=10)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("blue", "yellow", "red"))(50),
         cluster_cols =F,
         fontsize=6.3,
         fontsize_row=6.3,
         #scale="row",
         show_colnames=F,
         show_rownames=F,
         fontsize_col=3)
dev.off()

