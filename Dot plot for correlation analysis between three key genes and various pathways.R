
rm(list=ls())
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
immFile="cor_merge_nohallmark.txt"       #?????ļ?
setwd("D:\\缺血再灌注课题\\正式绘图_ER_stress\\15.关键hallmark_三基因相关性热图")     #???ù???Ŀ¼

data=read.table(immFile, header=T, sep="\t", check.names=F)

data$Catogery <- factor(data$Catogery,levels = c("GGSE43974","GSE90861","GSE126805"))
ggplot(data, aes(x = Immune, y = Gene, size = text_new, color=cor)) + 
  geom_point()+scale_color_gradient(high="red",low="blue")+
  ggprism::theme_prism(palette = "black_and_white",
                       base_size = 10, 
                       #base_line_size = 1,#base_family = "Roboto Condensed",
                       axis_text_angle = 45,
                       border = TRUE,base_rect_size = 2)+facet_grid(Catogery~., scale = 'free_x', space = 'free_x')











