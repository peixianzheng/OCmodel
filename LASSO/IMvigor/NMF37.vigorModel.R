

#install.packages("survival")
#install.packages("caret")
#install.packages("glmnet")
#install.packages("survminer")
#install.packages("timeROC")


#引用包
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

setwd("C:\\biowolf\\NMF\\37.IMvigor210")      #设置工作目录
rt=read.table("lasso.SigExp.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件

#模型构建
multiCox <- coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
score=predict(multiCox, type="risk", newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`", "", coxGene)
	
#读取IMvigor210表达数据
vigor=read.table("exp.txt", header=T, sep="\t", check.names=F, row.names=1)
vigor=t(vigor[coxGene,])

#读取生存数据文件
cli=read.table("time.txt", header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")

#合并IMvigor210数据
sameSample=intersect(row.names(vigor), row.names(cli))
vigorTime=cbind(cli[sameSample,,drop=F], vigor[sameSample,,drop=F])
#vigorTime[,3:ncol(vigorTime)]=vigorTime[,3:ncol(vigorTime)]*median(as.matrix(rt[,coxGene]))/median(as.matrix(vigorTime[,3:ncol(vigorTime)]))

#输出IMvigor210的风险值
vigorScore=predict(multiCox, type="risk", newdata=vigorTime)
Risk=as.vector(ifelse(vigorScore>median(vigorScore), "high", "low"))
vigorRiskOut=cbind(vigorTime, riskScore=as.vector(vigorScore), Risk)
vigorRiskOut=cbind(id=rownames(vigorRiskOut), vigorRiskOut)
write.table(vigorRiskOut,file="risk.IMvigor.txt",sep="\t",quote=F,row.names=F)

#绘制生存曲线函数
bioSurvival=function(inputFile=null, outFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#比较高低风险组生存差异，得到显著性p值
	diff=survdiff(Surv(futime, fustat) ~Risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
		
	#绘制生存曲线
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 2,
		           palette=c("red", "blue"),
		           risk.table=F,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)
	pdf(file=outFile,onefile = FALSE,width = 6,height =5)
	print(surPlot)
	dev.off()
}

#定义绘制ROC曲线函数
bioROC=function(inputFile=null, outFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#ROC曲线
	ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
	               marker=rt$riskScore,cause=1,
	               weighting='aalen',
	               times=c(1,3,5),ROC=TRUE)
	pdf(file=outFile, width=5, height=5)
	plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("green",'blue','red'),lwd=2,bty = 'n')
	dev.off()
}

#调用函数，绘制生存曲线和ROC曲线
bioSurvival(inputFile="risk.IMvigor.txt", outFile="sur.IMvigor.pdf")
bioROC(inputFile="risk.IMvigor.txt", outFile="ROC.IMvigor.pdf")



