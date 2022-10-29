

#install.packages("corrplot")

library(corrplot)
setwd("C:\\Users\\lexb4\\Desktop\\m6A\\09.correlation")      #设置工作目录
rt=read.table("m6Aexp.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=t(rt)
res1 <- cor.mtest(rt, conf.level = 0.95)

pdf("correlation.pdf",height=8,width=8)              #保存图片的文件名称
corrplot(corr=cor(rt),
         method = "circle",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         p.mat = res1$p,
         sig.level = 0.001,
         insig = "pch",
         number.cex = 1,
         type = "upper",
         col=colorRampPalette(c("blue", "white", "red"))(50),
         )
dev.off()
