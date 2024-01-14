
rm(list = ls())

library(caret)
library(glmnet)
library(reshape2)
library(ggsci)
setwd("")  

rt=read.table("GSE43974_merge_fustat.txt",header=T,sep="\t",row.names=1,check.names = F)   
rt <- t(rt)

gene=read.table("three_gene.txt",header=T,check.names = F)
rt1=rt[,c("group",as.vector(gene[,1]))]
head(rt1)
rt1 <- as.data.frame(rt1)
data_long1=melt(rt1,id.vars=c("group")) #需保留的不参与聚合的变量列名
data_long1=data_long1 %>%mutate(group2=rep("GSE43974",nrow(data_long1)))
head(data_long1)
    
rt=read.table("GSE90861_merge_fustat.txt",header=T,sep="\t",row.names=1,check.names = F)   
rt <- t(rt)
gene=read.table("three_gene.txt",header=T,check.names = F)
rt2=rt[,c("group",as.vector(gene[,1]))]

rt2 <- as.data.frame(rt2)
data_long2=melt(rt2,id.vars=c("group")) #需保留的不参与聚合的变量列名
data_long2=data_long2 %>%mutate(group2=rep("GSE90861",nrow(data_long2)))
head(data_long2)
                               
rt=read.table("GSE126805_merge_fustat.txt",header=T,sep="\t",row.names=1,check.names = F)   
rt <- t(rt)
gene=read.table("three_gene.txt",header=T,check.names = F)
rt3=rt[,c("group",as.vector(gene[,1]))]

rt3 <- as.data.frame(rt3)
data_long3=melt(rt3,id.vars=c("group")) #需保留的不参与聚合的变量列名
data_long3=data_long3 %>%mutate(group2=rep("GSE126805",nrow(data_long3)))
head(data_long3)

data_long4 <- rbind(data_long1,data_long2,data_long3)
data_long4$value <- as.numeric(data_long4$value)#注意value的格式是什么，需要改成数字格式
#调节变量顺序
data_long4$group2 <- factor(data_long4$group2,levels = c("GSE43974", "GSE90861", "GSE126805"))

library(ggpubr)
p=ggboxplot(data_long4, x="variable", y="value", color = "group",width=0.8, 
            ylab="Gene expression",
            xlab="",
            legend.title="variable",
            palette = c("#00468BFF","#AD002AFF"),
            add = "none")+facet_grid(.~group2)


p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0,0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", " ")),
                        label = "p.signif")

p2 <- p1 + ggprism::theme_prism(palette = "black_and_white",
                                base_size = 12, 
                                #base_line_size = 1,#base_family = "Roboto Condensed",
                                axis_text_angle = 0,
                                border = F,base_rect_size = 1.2)



