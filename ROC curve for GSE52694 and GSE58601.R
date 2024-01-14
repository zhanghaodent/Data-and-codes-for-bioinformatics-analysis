
rm(list=ls())

library(pROC)

setwd("D:\\抗体介导排斥课题\\生信分析\\8.新加临床相关绘图")  

rt=read.table("GSE52697_merge.txt",header=T,sep="\t",row.names=1)   

View(rt)
rt=rt[1:13,]


gfit.train <- roc(Group~HALLMARK_INTERFERON_GAMMA_RESPONSE, data = rt,smooth=F)
pdf(file="GSE52697_ROC_2month.pdf",width=5.45,height=4.75)
plot(gfit.train, print.auc=TRUE,main = "", col= "#AD002AFF",identity.col="black",
     identity.lty=1,identity.lwd=1)
text(0.39, 0.43, "95% CI: 0.867-0.979", col="#AD002AFF")
dev.off()
