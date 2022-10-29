

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
riskFile="risk.IMvigor.txt"      #风险文件
cliFile="clinical.txt"           #临床数据文件
setwd("C:\\biowolf\\NMF\\38.IMvigorCliCor")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,"riskScore",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk, cli)

#临床相关性分析，输出图形结果
for(clinical in colnames(rt)[2:ncol(rt)]){
	data=rt[c("riskScore", clinical)]
	colnames(data)=c("riskScore", "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	#设置比较组
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	#绘制箱线图
	boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
		          xlab="",
		          ylab="Risk score",
		          legend.title=clinical,
		          add = "jitter")+ 
	    stat_compare_means(comparisons = my_comparisons)
	    #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	#输出图片
	pdf(file=paste0("cliCor.", clinical, ".pdf"), width=6, height=5)
	print(boxplot)
	dev.off()
}




