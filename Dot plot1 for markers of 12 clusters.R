####merge the scRNA-data####
rm(list=ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
library(clusterProfiler)


plan("multiprocess", workers = 3) ###set the compute core
options(future.globals.maxSize = 10000 * 1024^2)
getwd()

setwd("D:\\ischemic_reperfusion\\GSE161201\\2.merge")
control<-readRDS("control.rds")
IRI2<-readRDS("IRI6h.rds")
TR3 <- readRDS("IRId1.rds")

control@meta.data[1:5,1:5]
#####取小的对象###如果电脑内存较小的话只用1000个细胞，这里我们不用
SRR780<-SRR780[,1:1000]
SRR781<-SRR781[,1:1000]
SRR782<-SRR782[,1:1000]
SRR783<-SRR783[,1:1000]
##########


sce.mergeTEN<- merge(control,y=c(IRI2,TR3),project = "scTEN")

sce.mergeTEN@meta.data[1:5,]

rm(SRR780,SRR781,SRR782,SRR783)
gc()

#saveRDS(sce.mergeTEN,file="mergeTEN_2_after_anno.rds")

sce.mergeTEN<-readRDS(file="mergeTEN_2_after_anno.rds")

####integration with harmony####
library(devtools)
install_github("immunogenomics/harmony")
library(harmony)
gc()

sce.mergeTEN@meta.data[1:5,]


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sce.mergeTEN <- CellCycleScoring(sce.mergeTEN, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sce.mergeTEN@meta.data[1:5,]

##############查看各个样本之间差异,批次矫正前
sce.mergeTEN<-NormalizeData(sce.mergeTEN,verbose = T) 
sce.mergeTEN<-FindVariableFeatures(sce.mergeTEN,selection.method = "vst", nfeatures = 2000)
sce.mergeTEN<-ScaleData(sce.mergeTEN,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = FALSE)
sce.mergeTEN<-RunPCA(sce.mergeTEN,verbose = T,npcs = 50)
ElbowPlot(sce.mergeTEN,ndims = 50)
p1 <- DimPlot(object = sce.mergeTEN, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sce.mergeTEN, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
########开始批次矫正
sce.mergeTEN<-RunHarmony(sce.mergeTEN,group.by.vars = c("orig.ident"), plot_convergence = TRUE)

harmony_embeddings <- Embeddings(sce.mergeTEN, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = sce.mergeTEN, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = sce.mergeTEN, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))

###########批次矫正前后测对比
p1+p3

sce.mergeTEN <- sce.mergeTEN %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  RunTSNE(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)



sce.mergeTEN<-FindClusters(sce.mergeTEN,resolution = 1.0)

table(Idents(sce.mergeTEN))
Idents(sce.mergeTEN)<-sce.mergeTEN$RNA_snn_res.1
sce.mergeTEN@meta.data[1:5,]
p1<-DimPlot(sce.mergeTEN,reduction = "tsne",label = T)
p1
DimPlot(sce.mergeTEN,reduction = "umap",label = T)
DimPlot(sce.mergeTEN,reduction = "tsne",label = T,split.by = "tissue_type")
DimPlot(sce.mergeTEN,reduction = "umap",label = T,split.by = "tissue_type")
table(sce.mergeTEN@meta.data$cell_cluster_origin)



#Idents(sce.mergeTEN)<-sce.mergeTEN$RNA_snn_res.1.5
#VlnPlot(sce.mergeTEN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)
p1+p2

features = c("Cd79a","Cd79b","Atp6v1g3","Slc12a3","Kdr","Emcn","Pdgfra",
             "Pdgfrb","Slc12a1","Fcgr1","Nphs1","Nphs2","Slc13a3","Slc17a3","Slc22a30",
             "Cd247","S100a9","S100a8")

features = c("Nrp1","Kdr","Nphs1","Nphs2","Slc27a2","Lrp2","Slc12a3",
             "Pvalb","Aqp2","Hsd11b2","Atp6v1g3","Atp6v0d2","Insrr",
             "Rhbg","Mki67","Cdca3","Plac8","S100a4","C1qb","C1qa",
             "S100a8","S100a9","Cd79a","Cd79b","Ltb","Cxcr6","Gzma",
             "Nkg7","Stmn1","Cldn1","Spp2","Sptssb","Aqp1","Aqp4","C1qc","Cnn1",
             "Dcn","Thy1","Umod","Slc12a1","Emcn","Cd3d","Cd3g","Cd247","Slc34a1","Malat1",
             "Fbln5","Akap12","Rgs5","Tpm2","Acsm2","Slc4a4","Ccl5","Ctsw","Gimap3")

markers = c("Slc27a2","Lrp2","Acsm2",#PT
            "Slc12a3","Umod",#DCT
            "Aqp2","Hsd11b2",#CD-PC
            "Atp6v1g3","Atp6v0d2",#CD-IC
            "Nphs1","Nphs2",#Podocyte
            "Rgs5","Tpm2",#Mesangial cell
            "Akap12",#Medullary cell
            "S100a8","S100a9",#Neutrophil
            "C1qb","C1qa",#Macrophage
            "Ccl5","Ctsw",#NK
            "Plac8","S100a4",#Fibroblast
            "Kdr","Emcn","Nrp1","Ptprb")#Endothelium

markers <- c("Lcn2","Havcr1","Igfbp7","Timp2","Il18","Chi3l1","Spp1","S100A9")
Idents(sce.mergeTEN)<-sce.mergeTEN$Seurat_harmony

DotPlot(sce.mergeTEN,features = markers,cols = c("grey","darkblue"))+RotatedAxis()



############Dotplot_anno是封装好的函数，进行marker热图的绘制
library(Seurat)
library(ggplot2)
#BiocManager::install("dittoSeq")
library(dittoSeq)
library(viridis)
source('D:\\ischemic_reperfusion\\GSE161201\\2.merge\\Dotplot_anno.R')
#单细胞的效果
DefaultAssay(sce.mergeTEN) <- "RNA"
Dotplot_anno(sce.mergeTEN, features = markers, celltype_color = dittoColors(),
             group = c(rep('PT',3), rep('DCT',2),rep('CD-PC',2),
                       rep('CD-IC',2), rep('Pod',2),rep("Mes",2), rep('Med',1),rep('Neu',2),rep('Mac',2),
                       rep('NK',2),rep('Fib',2),rep('Endo',4)),
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             order = T)


