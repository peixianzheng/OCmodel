

#install.packages("survival")
#install.packages("NMF")


#引用包
library(survival)
library(NMF)

setwd("C:\\biowolf\\NMF\\16.NMF")
rt=read.table("TCGA.expTime.txt", header=T, sep="\t", check.names=F, row.names=1)

#单因素COX分析
sigGenes=c()
for(i in colnames(rt)[3:ncol(rt)]){
	cox=coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary=summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<0.01){ sigGenes=c(sigGenes,i) }
}

#NMF分析
data=t(rt[,sigGenes])
res=nmf(data, rank=2:10, method="brunet", nrun=10, seed=123456)
pdf(file="cophenetic.pdf", width=8, height=7, onefile=F)
plot(res)
dev.off()

#输出所有分型的热图
pdf(file="heatmap.all.pdf", width=15, height=15, onefile=F)
consensusmap(res,
             annRow=NA,
             annCol=NA,
             #tracks=c("consensus:"),
             main="Consensus matrix",
             info=FALSE)
dev.off()


#输出分型的结果
clusterNum=2        #分几类，根据判断标准判断
res=nmf(data, rank=clusterNum, method="brunet", nrun=10, seed=123456)
Cluster=predict(res)
Cluster=as.data.frame(Cluster)
Cluster$Cluster=paste0("C", Cluster$Cluster)
clusterOut=rbind(ID=colnames(Cluster), Cluster)
write.table(clusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)

#输出单独分型的热图
pdf(file="heatmap.pdf", width=6, height=6, onefile=F)
consensusmap(res,
             annRow=NA,
             annCol=NA,
             #tracks=c("consensus:"),
             main="Consensus matrix", 
             info=FALSE)
dev.off()




