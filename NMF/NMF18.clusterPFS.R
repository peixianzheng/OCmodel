

#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)

clusterFile="cluster.txt"      #分型结果文件
cliFile="Survival_SupplementalTable_S1_20171025_xena_sp"      #临床数据文件
setwd("C:\\biowolf\\NMF\\18.clusterPFS")        #设置工作目录

#读取分型结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("PFI.time", "PFI")]
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365
cli=as.matrix(cli)
row.names(cli)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(cli))

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#生存差异统计
length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit=survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Cluster",
		           legend.labs=levels(factor(rt[,"Cluster"])),
		           legend = c(0.8, 0.8),
		           font.legend=10,
		           xlab="Time(years)",
		           ylab="Progression free survival",
		           break.time.by = 1,
		           palette = bioCol,
		           #surv.median.line = "hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.3)

#输出图形
pdf(file="PFS.pdf",onefile=FALSE, width=7, height=6)
print(surPlot)
dev.off()



