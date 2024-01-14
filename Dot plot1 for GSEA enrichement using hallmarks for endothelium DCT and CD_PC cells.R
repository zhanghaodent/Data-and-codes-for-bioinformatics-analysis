######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("Seurat")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle")


#读取数据
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(monocle)

setwd("C:\\biowolf\\sCell\\13.scRNAdiff")             #设置工作目录

#读取对照组的数据,并对数据进行整理
rt=read.table("control.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
conData=avereps(data)
colnames(conData)=paste0("C.", colnames(conData))

#读取实验组的数据,并对数据进行整理
rt=read.table("treat.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
treatData=avereps(data)
colnames(treatData)=paste0("T.", colnames(treatData))

#数据合并
sameGene=intersect(row.names(conData), row.names(treatData))
data=cbind(conData[sameGene,], treatData[sameGene,])

#将矩阵转换为Seurat对象，并对数据进行过滤
pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_")
#使用PercentageFeatureSet函数计算线粒体基因的百分比
pbmc[["percent.mt"]]=PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmcCon=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #对数据进行过滤

#对数据进行标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)

#组间差异分析
logFCfilter=1
adjPvalFilter=0.05
groups=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
names(groups)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
pbmc.markers=FindMarkers(pbmc, ident.1 = "T", ident.2 = "C", group.by = 'group')
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
write.table(sig.markers,file="diffGene.txt",sep="\t",row.names=F,quote=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
