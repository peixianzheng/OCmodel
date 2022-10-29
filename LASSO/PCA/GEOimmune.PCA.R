

#install.packages("ggplot2")

#pca analysis
setwd("C:\\Users\\lexb4\\Desktop\\GEOimmune\\13.PCA")                   #设置工作目录
data=read.table("CIBERSORT.filter.txt",header=T,sep="\t",row.names=1)   #读取表格
data=as.matrix(data)                                                    #矩阵转置
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)                                 #PCA分析
write.table(predict(data.pca),file="newTab.xls",quote=F,sep="\t")       #输出新表

#pca 2d plot
library(ggplot2)
group=c(rep("Normal",6),rep("Tumor",6))  #对照组和实验组的样品数目
pcaPredict=predict(data.pca)
PCA = data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)
PCA.mean=aggregate(PCA[,1:2],list(group=PCA$group),mean)

#自定义函数
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(PCA$group)){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$group==g,],
                  veganCovEllipse(cov.wt(cbind(PCA1,PCA2),
                  wt=rep(1/length(PCA1),length(PCA1)))$cov,
                  center=c(mean(PCA1),mean(PCA2))))),group=g))
}

pdf(file="PCA.pdf",height=6,width=7)
ggplot(data = PCA, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
    geom_path(data=df_ell, aes(x=PCA1, y=PCA2,colour=group), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PCA1,y=PCA.mean$PCA2,label=PCA.mean$group)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

