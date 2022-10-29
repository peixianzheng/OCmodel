

#install.packages("survivalROC")

library(survivalROC)
setwd("C:\\Users\\lexb4\\Desktop\\m6A\\23.ROC")      #设置工作目录
rt=read.table("lassoRisk.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取lasso回归风险文件
pdf(file="ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
      predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
  xlab="False positive rate", ylab="True positive rate",
  main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

