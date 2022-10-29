

#install.packages("pheatmap")

library(pheatmap)
setwd("D:\\biowolf\\methySubtype\\15.cliHeatmap")                        #���ù���Ŀ¼
rt=read.table("multiIndepSigExp.txt",sep="\t",header=T,row.names=1,check.names=F)      #��ȡ�׻����ļ�
rt=t(rt[,c(-1,-2)])

#��ȡ�ٴ����ݺͷ����ļ�
clinical=read.table("clinical.txt",sep="\t",header=T,row.names=1,check.names=F)
cluster=read.table("cluster.txt",sep="\t",header=T,row.names=1,check.names=F)
Type=clinical[row.names(cluster),]
Type=cbind(Type,Cluster=paste0("C",cluster[,"cluster"]))
#���շ��Ͷ���Ʒ����
Type=Type[order(Type$Cluster),]
rt=rt[,row.names(Type)]
#������ͼ
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
