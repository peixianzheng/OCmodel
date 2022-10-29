#install.packages("ggpubr")

library(ggpubr)

setwd("C:\\Users\\biowolf\\Desktop\\ssGSEA\\14.HLAboxplot")                          #设置工作目录
rt=read.table("HLAexp.txt",sep="\t",header=T,row.names=1,check.names=F)    #读取文件

Type=read.table("cluster.Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[order(Type[,2]),]
rt=t(rt[,row.names(Type)])

#准备箱线图的输入文件
data=data.frame()
for(i in colnames(rt)){
  data=rbind(data,cbind(expression=log2(rt[,i]+1),gene=i,Subtype=as.vector(Type[,2])))
}
write.table(data,file="data.txt",sep="\t",row.names=F,quote=F)

#绘制箱型图
data=read.table("data.txt",sep="\t",header=T,check.names=F)       #读取箱线图输入文件
data$Subtype=factor(data$Subtype, levels=c("Immunity_L","Immunity_M","Immunity_H"))
p=ggboxplot(data, x="gene", y="expression", color = "Subtype", orientation = "horizontal",
     ylab="Gene expression",
     xlab="",
     palette = c("green","blue","red") )
p=p+rotate_x_text(60)
pdf(file="boxplot.pdf",width=6,height=6)                          #输出图片文件
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

