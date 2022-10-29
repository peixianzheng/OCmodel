

#install.packages("ggpubr")

setwd("C:\\Users\\biowolf\\Desktop\\ssGSEA\\12.estimateVioplot")      #设置工作目录

library(ggpubr)
Type=read.table("cluster.Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
#Type=Type[order(Type[,2]),]

score=read.table("scores.txt",sep="\t",check.names=F,row.names=1,header=T)
score=score[row.names(Type),]
colnames(Type)=c("cluster","Subtype")
cluster=cbind(Type,score)
cluster=cluster[,-1]

cluster$Subtype=factor(cluster$Subtype, levels=c("Immunity_L","Immunity_M","Immunity_H"))
my_comparisons=list(c("Immunity_L","Immunity_M"),c("Immunity_M","Immunity_H"),c("Immunity_H","Immunity_L"))

pdf(file="vioplot.pdf",width=6,height=5)
ggviolin(cluster, x="Subtype", y="TumorPurity", fill = "Subtype", 
         palette = c("green", "blue", "red"), 
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()
