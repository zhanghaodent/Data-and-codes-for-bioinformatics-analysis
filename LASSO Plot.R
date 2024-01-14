
rm(list = ls())

library(caret)
library(glmnet)
library(reshape)
library(ggsci)
setwd("")  

set.seed(123)
rt=read.table("",header=T,sep="\t",row.names=1,check.names = F)   
rt <- t(rt)

gene=read.table("selected_gene.txt",header=T,check.names = F)
train=rt[,c("fustat",as.vector(gene[,1]))]

x=as.matrix(train[,c(2:ncol(train))])
y = train[,1]
cv_fit <- cv.glmnet(x=x, y=y, alpha = 1,family="binomial")
plot(cv_fit)
abline(a=)
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se,family="binomial")

coef=coef(model_lasso_min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef_new_1se.txt",sep="\t",quote=F,row.names=F)

#求生存数据的每个样本的风险得分
rt=read.table("GSE126805_merge_fustat.txt",header=T,sep="\t",row.names=1)
rt <- t(rt)

#下面进行一次就可以不运行了，为了调节协变量
lassoGene1 <- lassoGene[-1]
actCoef1 <- actCoef[-1]


myFun=function(x){crossprod(as.numeric(x),actCoef1)+actCoef[1]}
testFinalGeneExp=rt[,lassoGene1] 
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("fustat",lassoGene1)
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="GSE126805_risk.txt",sep="\t",quote=F,row.names=F)

#求非生存数据的每个样本的风险得分
rt=read.table("GSE21374_merge_survival_t.txt",header=T,sep="\t",row.names=1)
rt <- t(rt)
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="9.9GSE13507_risk.txt",sep="\t",quote=F,row.names=F)


#????train??????ֵ
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="TCGA_risk.txt",sep="\t",quote=F,row.names=F)




#进行lasso第二种图的绘制
model_lasso_min_n <- glmnet(x=x, y=y, alpha = 1, nlambda = 100,family="binomial")

pdf("3.lasso_lambda_new.pdf")
plot(model_lasso_min_n,xvar = "lambda",label=T)+
  abline(v=log(cv_fit$lambda.1se),col = "gray", lwd = 2, lty = 2)



x <- coef(model_lasso_min_n)  
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- model_lasso_min_n$lambda[tmp$variable+1] # extract the lambda values
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm  


ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cv_fit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Lambda (log scale)") + 
  #xlab("L1 norm")+
  ylab('Coefficients')+
  theme_bw(base_rect_size = 2)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(10),pal_aaas()(10),pal_gsea()(30)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+
  #annotate('text',x = -3.3,y=0.26,label='Optimal Lambda = 0.012',color='black')+
  guides(col=guide_legend(ncol = 1))





