######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("survival")
#install.packages("survminer")
rm(list=ls())

library(survival)
library("survminer")
setwd("")              #???ù???Ŀ¼

		rt=read.table("GSE21374_threegene_input.txt",header=T,sep="\t")                   #??ȡ?????ļ?
		#?Ƚϸߵͷ????????????죬?õ???????pֵ
		diff=survdiff(Surv(futime, fustat) ~ATF3_group,data = rt)
		pValue=1-pchisq(diff$chisq,df=1)
		pValue=signif(pValue,4)
		pValue=format(pValue, scientific = TRUE)
		fit <- survfit(Surv(futime, fustat) ~ ATF3_group, data = rt)
		
		#????????????
		surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=paste0("p=",pValue),
		           pval.size=4,
		           risk.table=F,
		           legend.labs=c("Low risk","High risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 1,
		           risk.table.title="",
		           palette=c("#0072B5FF","#AD002AFF"),
		           risk.table.height=.25,censor=F,legend="right")
surPlot

ggsave("ATF3_KM.pdf",width = 7,height = 5)




