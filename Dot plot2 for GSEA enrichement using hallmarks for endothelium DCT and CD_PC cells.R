

rm(list=ls())
#???ð?
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggsci)

setwd("D:\\抗体介导排斥课题\\生信分析\\1.肾活检\\6.GSEA_NA_ABMR")      #???ù???Ŀ¼
dd <- read.table("小彭_new.txt",header=T,sep="\t",check.names=F)


geneList <- dd$logFC
## 2.命名
names(geneList) = dd$id
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
#write.table(geneList,"geneList.txt",sep = "\t")
#####运行GSEA分析
## 读入hallmarks gene set
hallmarks <- read.gmt("h.all.v2022.1.Hs.symbols.gmt")
# 需要网络
kk <- GSEA(
  geneList,
  exponent = 1,
  minGSSize = 5,
  maxGSSize = 1000,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = hallmarks,
  verbose = TRUE,
  seed = "0404",
  by = "fgsea"
)

kkTab=as.data.frame(kk)
write.table(kkTab, file="GSEA_hallmark.xls", sep="\t", quote=F, row.names=F)

#对GSEA结果进行绘图，pvalue可不进行展示
gseaplot2(kk,geneSetID = c("GOBP_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",  
                           "GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS"),pvalue_table = F,
          rel_heights = c(0.8,0.2,0.5),base_size=11,
          color = pal_npg("nrc")(10))





