
#绘制分面箱线图

rm(list = ls())

library(caret)
library(glmnet)
library(reshape2)
library(ggsci)
setwd("")  

rt=read.table("barplot_input.txt",header=T,sep="\t",row.names=1,check.names = F)   
rt <- rt[1:13,]

rt1 <- as.data.frame(rt)
data_long1=melt(rt1,id.vars=c("group")) #需保留的不参与聚合的变量列名
head(data_long1)
data_long1$variable <-factor(data_long1$variable,levels =c("PPP1R15A","JUN","ATF3"))

library(ggpubr)
p=ggboxplot(data_long1, x="variable", y="value", color = "group",width=0.8, 
            ylab="Gene expression",
            xlab="",
            legend.title="variable",
            palette = c("#00468BFF","#AD002AFF"),
            add = "jitter")


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
ggsave("2month_threegene_new.pdf",height=5,width=7)








