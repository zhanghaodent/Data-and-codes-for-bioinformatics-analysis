

rm(list = ls())

library(ggpubr)
library(ggplot2)
library(reshape)

setwd("")      #???ù???Ŀ¼
rt=read.table("GSE98622_target_gene.txt",header=T,sep="\t",check.names=F,row.names=1)     #??ȡ?????ļ?

rt1 <- melt(rt)

ggline(rt1, x="group", y="value", add = "mean_se", color = "variable", 
       palette = "lancet")+ 
  stat_compare_means(aes(group=variable), label = "p.signif", hide.ns = TRUE,
                     label.y = c(3.7, 4.6, 5.8))+
  theme(legend.title=element_blank())+   ##把所有图例的标题去掉
  theme(legend.background = element_rect(color = "steelblue", 
                                         linetype = "solid", size = 0.1, ## 添加外框
                                         fill = "white"), ## 填充颜色
        legend.key = element_rect(fill = "lightblue"), ## 修改示例填充色
        legend.key.size = unit(.5, "cm"),## 修改文字高度
        legend.key.width = unit(1,"cm"))+## 修改文字宽度
  theme(legend.title=element_blank())   ##把所有图例的标题去掉




