####merge the scRNA-data####
rm(list=ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
library(clusterProfiler)

setwd("")
#saveRDS(sce.mergeTEN,file="mergeTEN_2_after_anno.rds")

sce.mergeTEN<-readRDS(file="mergeTEN_2_after_anno.rds")

#Idents(sce.mergeTEN)<-sce.mergeTEN$RNA_snn_res.1.5
#VlnPlot(sce.mergeTEN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)
sce.mergeTEN@meta.data[1:5,]
features = c("Ppp1r15a","Jun","Atf3")
Idents(sce.mergeTEN)<-sce.mergeTEN$Seurat_harmony
library(ggplot2)
library(dittoSeq)
library(ggpubr)

p1 <- VlnPlot(sce.mergeTEN,features="Ppp1r15a",idents = c("Endothelium","DCT","DC-PC"),
        split.by="orig.ident",pt.size = 0)+stat_compare_means(label = "p.signif")+
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE)

p2 <- VlnPlot(sce.mergeTEN,features="Jun",idents = c("Endothelium","DCT","DC-PC"),
              split.by="orig.ident",pt.size = 0)+stat_compare_means(label = "p.signif")+
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE)

p3 <- VlnPlot(sce.mergeTEN,features="Atf3",idents = c("Endothelium","DCT","DC-PC"),
              split.by="orig.ident",pt.size = 0)+stat_compare_means(label = "p.signif")+
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE)


DotPlot(sce.mergeTEN,features = markers,cols = c("grey","darkblue"))+RotatedAxis()

##########绘制目标基因的密度图
library("Nebulosa")
plot_density(sce.mergeTEN, "Ppp1r15a")########没有split.by参数

