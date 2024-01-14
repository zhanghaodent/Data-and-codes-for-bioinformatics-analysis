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
DimPlot(sce.mergeTEN,reduction = "tsne",label = T)

