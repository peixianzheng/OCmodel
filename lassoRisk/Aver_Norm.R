setwd("~/desktop/1")
tt=read.table( "symbol.txt",header = T,row.names =1)

  

ttMean=apply(tt,1,mean)

ttAVer=sweep(tt, 1, ttMean, '-')  ##每列减去其平均值

gene=rownames(ttAVer)
ttRes=cbind(gene,ttAVer)
write.table( ttRes, file="symbol_Aver.txt",sep="\t",row.names=F ,quote=F)


