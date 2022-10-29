
#install.packages("ggpubr")

library(ggpubr)

pFilter=0.05
setwd("C:\\Users\\biowolf\\Desktop\\ssGSEA\\19.CIBERSORTcor")                           #设置工作目录
rt=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)    #读取文件
data=rt[rt[,"P-value"]<0.05,]
data=data[,1:(ncol(rt)-3)]

Type=read.table("cluster.Immunity.txt",sep="\t",check.names=F,row.names=1,header=F)
Type=Type[row.names(data),]
colnames(Type)=c("cluster","Subtype")

outTab=data.frame()
data=cbind(data,Type)
for(i in colnames(data[,1:(ncol(data)-2)])){
  rt1=data[,c(i,"Subtype")]
  colnames(rt1)=c("expression","Subtype")
  ksTest<-kruskal.test(expression ~ Subtype, data = rt1)
  pValue=ksTest$p.value
  if(pValue<pFilter){
      outTab=rbind(outTab,cbind(rt1,gene=i))
      print(pValue)
  }
}
write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)

#绘制箱型图
data=read.table("data.txt",sep="\t",header=T,check.names=F)       #读取箱线图输入文件
data$Subtype=factor(data$Subtype, levels=c("Immunity_L","Immunity_M","Immunity_H"))
p=ggboxplot(data, x="gene", y="expression", color = "Subtype",
     ylab="Fraction",
     xlab="",
     palette = c("green","blue","red") )
p=p+rotate_x_text(45)
pdf(file="boxplot.pdf",width=7.5,height=5.5)                          #输出图片文件
p+stat_compare_means(aes(group=Subtype),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

