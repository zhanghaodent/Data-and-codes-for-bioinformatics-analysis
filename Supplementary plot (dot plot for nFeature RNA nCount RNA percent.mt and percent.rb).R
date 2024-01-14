
rm(list=ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
library(clusterProfiler)

setwd("D:\\ischemic_reperfusion\\GSE161201\\1.single_sample")

matrix<-read.table("GSM4891840_Homeo_umi_expression_matrix.tsv",header=T,row.names = 1)
library(stringr)
genename<-str_replace_all(colnames(matrix),pattern = "_",replacement = "-")
colnames(matrix)<-genename

#colnames(matrix)<-genename
#matrix[1:5,1:5]
#matrix<-t(matrix)
#matrix[1:5,1:5]
IRIsham_object<- CreateSeuratObject(counts = matrix,project = "control", min.cells = 0, min.features = 200)

IRIsham_object[["percent.mt"]] <- PercentageFeatureSet(IRIsham_object, pattern = "^mt-")
hist(IRIsham_object[["percent.mt"]]$percent.mt)
IRIsham_object[["percent.rb"]] <- PercentageFeatureSet(object = IRIsham_object, pattern = "^Rp[sl]") 


matrix<-read.table("GSM4891841_IRI6h_umi_expression_matrix.tsv",header=T,row.names = 1)
library(stringr)
genename<-str_replace_all(colnames(matrix),pattern = "_",replacement = "-")
colnames(matrix)<-genename

#colnames(matrix)<-genename
#matrix[1:5,1:5]
#matrix<-t(matrix)
#matrix[1:5,1:5]
IRIsham6h_object<- CreateSeuratObject(counts = matrix,project = "IRI6h", min.cells = 0, min.features = 200)

IRIsham6h_object[["percent.mt"]] <- PercentageFeatureSet(IRIsham6h_object, pattern = "^mt-")
hist(IRIsham6h_object[["percent.mt"]]$percent.mt)
IRIsham6h_object[["percent.rb"]] <- PercentageFeatureSet(object = IRIsham6h_object, pattern = "^Rp[sl]") 


sce.mergeTEN<- merge(IRIsham_object,y=c(IRIsham6h_object,IRIsham1b1_object),project = "scTEN")


VlnPlot(sce.mergeTEN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)









