

#install.packages("pheatmap")

setwd("C:\\Users\\lexb4\\Desktop\\m6A\\26.RiskClinicalHeatmap")      #设置工作目录
rt=read.table("riskCliExp.txt",sep="\t",header=T,row.names=1,check.names=F)    #读取文件
rt=t(rt)
outpdf="heatmap.pdf"

library(pheatmap)
Type=read.table("riskCliGroup.sig.txt",sep="\t",header=T,row.names=1,check.names=F)
Type=Type[order(Type$Risk),]
rt=rt[,row.names(Type)]

pdf(outpdf,height=6.3,width=10)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("green", "white", "red"))(50),
         cluster_cols =F,
         fontsize=7.5,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()


