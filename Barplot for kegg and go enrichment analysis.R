
rm(list=ls())

library(ggplot2)
library(ggsci)
setwd("")
data=read.table("",sep="\t",header=T,check.names=F)

data$ID <- factor(data$ID,levels = data$ID)

#lancet中"#AD002AFF"是红色，#00468BFF"是蓝色
#最终绘图代码
p <- ggplot(data)+
  geom_col(aes(ID,NES,fill=qvalues))+
  scale_fill_gradient2(low = "#00468B", high = "#AD002AFF",mid="#00468B")+
  #添加一条竖线
  theme_classic()+
  ylim(0,10)+
  coord_flip()+
  #调整主题
  theme(
    #去除图例
    legend.position="none",
    #标题居中
    plot.title=element_text(hjust=0.5),
    axis.line.y=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
  )+
  ylab("Normalized Enrichemnt Score")

p <- p + ggprism::theme_prism(palette = "black_and_white",
                              base_size = 8, 
                              #base_line_size = 1,#base_family = "Roboto Condensed",
                              axis_text_angle = 0,
                              border = TRUE,base_rect_size = 2)


#或者根据 Category 绘制分面图
p + facet_grid(Category~.,scale = 'free_y', space = 'free_y')

ggsave("5K_GSEA_5.pdf",height=5,width=7)


