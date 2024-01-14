###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######??????ѧ??: http://www.biowolf.cn/
######???????䣺2749657388@qq.com
######????΢??: 18520221056


install.packages('forestplot')

setwd("")
library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="green",line="darkblue", summary="royalblue")             #????ɭ??ͼ??ɫ
rt=read.table("pre_merge1.txt",header=T,sep="\t",check.names=F,row.names=1)

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 outTab=rbind(outTab,
              cbind(id=i,
              HR=coxSummary$conf.int[,"exp(coef)"],
              HR.95L=coxSummary$conf.int[,"lower .95"],
              HR.95H=coxSummary$conf.int[,"upper .95"],
              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
              )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)

#????ɭ??ͼ
rt=read.table("uniCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #????ͼƬ????
pdf(file="forest.pdf",
       width = 6,             #ͼƬ?Ŀ???
       height = 4,            #ͼƬ?ĸ߶?
       )
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4,
           xlab="Hazard ratio"
           )
dev.off()


###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######??????ѧ??: http://www.biowolf.cn/
######???????䣺2749657388@qq.com
######????΢??: 18520221056