

#install.packages("corrplot")


#���ð�
library(corrplot)

riskFile="risk.TCGAall.txt"             #�����ļ�
immuneFile="MCPcounter.result.txt"      #����ϸ��������
setwd("C:\\biowolf\\NMF\\35.riskImmCor")     #���ù���Ŀ¼

#��ȡ����ϸ���������ļ�
data=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

#��ȡ�����ļ�
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#���ݺϲ�
sameSample=intersect(row.names(data), row.names(risk))
data=cbind(risk[sameSample,"riskScore",drop=F], data[sameSample,,drop=F])

#����Ծ���
M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)

#���������ͼ��
pdf(file="immuneCor.pdf", width=7, height=7)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()



