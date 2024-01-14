######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("Seurat")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle")


#��ȡ����
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(monocle)

setwd("C:\\biowolf\\sCell\\13.scRNAdiff")             #���ù���Ŀ¼

#��ȡ�����������,�������ݽ�������
rt=read.table("control.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
conData=avereps(data)
colnames(conData)=paste0("C.", colnames(conData))

#��ȡʵ���������,�������ݽ�������
rt=read.table("treat.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
treatData=avereps(data)
colnames(treatData)=paste0("T.", colnames(treatData))

#���ݺϲ�
sameGene=intersect(row.names(conData), row.names(treatData))
data=cbind(conData[sameGene,], treatData[sameGene,])

#������ת��ΪSeurat���󣬲������ݽ��й���
pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_")
#ʹ��PercentageFeatureSet�����������������İٷֱ�
pbmc[["percent.mt"]]=PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmcCon=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #�����ݽ��й���

#�����ݽ��б�׼��
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#��ȡ��Щ��ϸ�������ϵ���ϴ�Ļ���
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)

#���������
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
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056