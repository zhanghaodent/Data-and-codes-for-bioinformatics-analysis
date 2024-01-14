
rm(list=ls())
library(survival)
library(timeROC)
library(ggsci)

setwd("")      #???ù???Ŀ¼
rt=read.table("Interferon_hallmark_merge_survival.txt",header=T,sep="\t",check.names=F,row.names=1)    #??ȡcox?ع??????ļ?

#lancet中的配色,"#00468BFF"为深蓝色，"#42B540FF"为绿色，"#925E9FFF"为紫色，"#AD002AFF"为红色
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,marker=rt$riskScore,cause=1,
               weighting='aalen',times=c(1,2,3,4),ROC=TRUE)
pdf(file="GSE21374_survivalROC.pdf",width=5,height=5)
plot(ROC_rt,time=1,col="#00468BFF",title=FALSE,lwd=2)
plot(ROC_rt,time=2,col="#42B540FF",title=FALSE,lwd=2,add=T)
plot(ROC_rt,time=3,col="#925E9FFF",title=FALSE,lwd=2,add=T)
plot(ROC_rt,time=4,col="#AD002AFF",title=FALSE,lwd=2,add=T)
legend("bottomright", 
       c(paste0("1-Year AUC = ",sprintf("%.03f",ROC_rt$AUC[1])),
       paste0('2-Year AUC = ',sprintf("%.03f",ROC_rt$AUC[2])),
       paste0("3-Year AUC = ",sprintf("%.03f",ROC_rt$AUC[3])),
       paste0("4-Year AUC = ",sprintf("%.03f",ROC_rt$AUC[4]))),
       col = c("#00468BFF","#42B540FF","#925E9FFF","#AD002AFF"),lwd=2,bty="n")
dev.off()




rm(list=ls())
library(survival)
library(timeROC)
library(ggsci)

setwd("D:\\缺血再灌注课题\\正式绘图_ER_stress\\22.KM曲线_ROC")      #???ù???Ŀ¼
rt=read.table("GSE21374_threegene_input.txt",header=T,sep="\t",check.names=F,row.names=1)    #??ȡcox?ع??????ļ?

#lancet中的配色,"#00468BFF"为深蓝色，"#42B540FF"为绿色，"#925E9FFF"为紫色，"#AD002AFF"为红色
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,marker=rt$ATF3,cause=1,
               weighting='aalen',times=c(1,2,3,4),ROC=TRUE)
pdf(file="ATF3_survivalROC_new.pdf",width=4.85,height=5)
plot(ROC_rt,time=1,col="#00468BFF",title=FALSE,lwd=2)
plot(ROC_rt,time=3,col="#925E9FFF",title=FALSE,lwd=2,add=T)
plot(ROC_rt,time=4,col="#AD002AFF",title=FALSE,lwd=2,add=T)
legend("bottomright", 
       c(paste0("1-Year AUC = ",sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('2-Year AUC = ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0("3-Year AUC = ",sprintf("%.03f",ROC_rt$AUC[3]))),
       col = c("#00468BFF","#925E9FFF","#AD002AFF"),lwd=2,bty="n")
dev.off()







setwd("D:\\肿瘤代谢与免疫-BLCA\\14.合并文件")      #???ù???Ŀ¼
rt=read.table("merge_risk.txt",header=T,sep="\t",check.names=F,row.names=1)    #??ȡcox?ع??????ļ?

ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,marker=rt$riskScore,cause=1,
               weighting='aalen',times=c(1,2,3,5,10),ROC=TRUE)
pdf(file="merge_ROC.pdf",width=6,height=6)
plot(ROC_rt,time=1,col="green",title=FALSE,lwd=2)
plot(ROC_rt,time=3,col="red",title=FALSE,lwd=2,add=T)
plot(ROC_rt,time=5,col="yellow",title=FALSE,lwd=2,add=T)
legend("bottomright", 
       c(paste0("AUC at 1 years: ",sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3])),
         paste0("AUC at 5 years: ",sprintf("%.03f",ROC_rt$AUC[4]))),
       col = c("green","red","yellow"),lwd=2,bty="n")
dev.off()



