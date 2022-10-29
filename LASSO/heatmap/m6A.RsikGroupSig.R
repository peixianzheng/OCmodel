

setwd("C:\\Users\\lexb4\\Desktop\\m6A\\25.RsikGroupSig")
field="Risk"
flag1="low"
flag2="high"

rt=read.table("riskCliGroup.txt",sep="\t",header=T,check.names=F)
trainFlag=rt[rt[,field]==flag1,]
trainFlag=cbind(trainFlag,flag="Group1")
testFlag=rt[rt[,field]==flag2,]
testFlag=cbind(testFlag,flag="Group2")
newTable=rbind(trainFlag,testFlag)

newLabels=c("id")
for(i in 2:(ncol(rt)-1) ){
  nameStat=colnames(newTable)[i]
  tableStat=table(newTable[,c(nameStat,"flag")])
  pStat=chisq.test(tableStat)
  pvalue=pStat$p.value
  if(pvalue<0.001){
	  newLabels=c(newLabels,paste0(colnames(newTable)[i],"***"))
	}else if(pvalue<0.01){
	  newLabels=c(newLabels,paste0(colnames(newTable)[i],"**"))
	}else if(pvalue<0.05){
	  newLabels=c(newLabels,paste0(colnames(newTable)[i],"*"))
	}else{
	  newLabels=c(newLabels,colnames(newTable)[i])
	}
	print(paste(colnames(newTable)[i],pvalue,sep=" "))
}
newLabels=c(newLabels,colnames(newTable)[ncol(rt)])
colnames(rt)=newLabels
write.table(rt,file="riskCliGroup.sig.txt",sep="\t",row.names=F,quote=F)
