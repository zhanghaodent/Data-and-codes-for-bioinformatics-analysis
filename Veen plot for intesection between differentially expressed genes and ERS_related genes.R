
# 安装并加载包
# install.packages("ggVennDiagram")
library(ggVennDiagram)
# 示例数据准备
genes <- paste("gene",1:1000,sep="")
set.seed(2022)
x <- list(A=sample(genes,300),
          B=sample(genes,525),
          C=sample(genes,440))

setwd("")
rt=read.table("veen_input.txt",sep="\t",header=T,check.names=F)
GSE43974 <- rt$GSE43974[1:219]
GSE90861 <- rt$GSE90861
GSE126805 <- rt$GSE126805[1:327]
ERS <- rt$ERS_gene[1:258]
x <- list(GSE43974=GSE43974,
          ERS=ERS)


library(ggplot2)
ggVennDiagram(x,category.names = c("GSE43974","ERS"),
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              edge_lty = "solid", 
              edge_size = 1) +
  scale_fill_gradient(low="white",high = "#AD002AFF",name = "gene count")





