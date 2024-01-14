rm(list=ls())
#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 6) ###set the compute core
options(future.globals.maxSize = 60000 * 1024^2)
getwd()
library(stringr)

#memory.limit()
#BiocManager::install("igraph",force = T)

setwd("D:\\ischemic_reperfusion\\GSE161201\\3.cellchat_individual")
####prepare the data for cellchat analysis####
sce.mergeTEN<-readRDS("mergeTEN_2_after_anno.rds")

table(Idents(sce.mergeTEN))

#########去除相应细胞
sce.mergeTEN<-sce.mergeTEN[,Idents(sce.mergeTEN)!="Doublets"]

sce.mergeTEN@meta.data[1:5,]

table(sce.mergeTEN$orig.ident)
levels(Idents(sce.mergeTEN))

normaldata<-sce.mergeTEN[,sce.mergeTEN$orig.ident=="control"]
tumordata<-sce.mergeTEN[,sce.mergeTEN$orig.ident=="IRI6h"]
tumordata_1d<-sce.mergeTEN[,sce.mergeTEN$orig.ident=="IRId1"]


##减少计算，只分析前面2000个细胞
normaldata<-normaldata[,1:2000]
table(Idents(normaldata))

tumordata<-tumordata[,1:2000]
table(Idents(tumordata))
###perform the cellchat analysis##

#devtools::install_github("sqjin/CellChat")
library(CellChat)

#########需要标准化以后的数据
normal.input <- GetAssayData(normaldata, assay = "RNA", slot = "data") # normalized data matrix

labels <- factor(normaldata$Seurat_harmony,levels=levels(Idents(sce.mergeTEN)))
labels

meta <- data.frame(group = labels, row.names = rownames(normaldata@meta.data)) # create a dataframe of the cell labels
meta$group<-str_replace_all(meta$group,pattern = "/",replacement = "_")
meta$group <- factor(meta$group,levels=c("Epithelial","T_NK","B","Plasma_B","Myeloid","MAST", "Fibroblast","Endothelial","Pericytes","Neuron"))

cellchat_normal <- createCellChat(object = normal.input, meta = meta, group.by = "group")
saveRDS(cellchat_normal,file="cellchat_normal.rds")


tumor.input <- GetAssayData(tumordata, assay = "RNA", slot = "data") # normalized data matrix
labels <- factor(tumordata$Seurat_harmony,levels=levels(Idents(sce.mergeTEN)))
meta <- data.frame(group = labels, row.names = rownames(tumordata@meta.data)) # create a dataframe of the cell labels
meta$group<-str_replace_all(meta$group,pattern = "/",replacement = "_")
meta$group <- factor(meta$group,levels=c("Epithelial","T_NK","B","Plasma_B","Myeloid","MAST", "Fibroblast","Endothelial","Pericytes","Neuron"))
cellchat_tumor <- createCellChat(object = tumor.input, meta = meta, group.by = "group")
saveRDS(cellchat_tumor,file="cellchat_tumor.rds")

tumor.input <- GetAssayData(tumordata_1d, assay = "RNA", slot = "data") # normalized data matrix
labels <- factor(tumordata_1d$Seurat_harmony,levels=levels(Idents(sce.mergeTEN)))
meta <- data.frame(group = labels, row.names = rownames(tumordata_1d@meta.data)) # create a dataframe of the cell labels
meta$group<-str_replace_all(meta$group,pattern = "/",replacement = "_")
meta$group <- factor(meta$group,levels=c("Epithelial","T_NK","B","Plasma_B","Myeloid","MAST", "Fibroblast","Endothelial","Pericytes","Neuron"))
cellchat_tumor <- createCellChat(object = tumor.input, meta = meta, group.by = "group")
saveRDS(cellchat_tumor,file="cellchat_tumor_1d.rds")





####prepare the deconvolution analysis gene signature####


rm(list=ls())
setwd("~/Documents_PC/scRNA-seq/Data")
normal_cellchat<-readRDS(file="cellchat_normal.rds")
tumor_cellchat<-readRDS(file="cellchat_tumor.rds")
tumor_cellchat_1d<-readRDS(file="cellchat_tumor_1d.rds")


normal_cellchat@DB <- CellChatDB.mouse
gc()##释放内存
normal_cellchat <- subsetData(normal_cellchat) # subset the expression data of signaling genes for saving computation cost
#future::plan("multiprocess", workers = 2) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
normal_cellchat <- identifyOverExpressedGenes(normal_cellchat)
normal_cellchat <- identifyOverExpressedInteractions(normal_cellchat)
normal_cellchat <- projectData(normal_cellchat, PPI.mouse)
normal_cellchat <- computeCommunProb(normal_cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
normal_cellchat <- filterCommunication(normal_cellchat, min.cells = 1)
normal_cellchat <- computeCommunProbPathway(normal_cellchat)
normal_cellchat <- netAnalysis_computeCentrality(normal_cellchat, slot.name = "netP")


saveRDS(normal_cellchat,file="normalcellchat_1.rds")
normal_cellchat<-readRDS("./Cellchat/normalcellchat_1.rds")


normal_cellchat <- aggregateNet(normal_cellchat)
groupSize_control <- as.numeric(table(normal_cellchat@idents))
groupSize_control
table(normal_cellchat@idents)

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(normal_cellchat@net$count,arrow.size = 0.01, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(normal_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- normal_cellchat@net$count

#或者看weight
mat <- normal_cellchat@net$weight
dev.off()

####下面针对其中一个细胞亚型分析其对其他细胞的interaction的Ligand-receptor数量
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[1, ] <- mat[1, ]
netVisual_circle(mat2, vertex.weight = groupSize,arrow.size = 0.2, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[1])


###查看具体的信号通路
df.net <- subsetCommunication(normal_cellchat)
#df.net_tumor <- subsetCommunication(tumor_cellchat)
View(df.net)

levels(df.net$source)

###########关注某些信号通路
normal_cellchat@netP$pathways#######关注其中有哪些信号通路
pathways.show <- c("GAS") 

#table(df.net$pathway_name)
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#dev.new()
seq(1,4)

vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(normal_cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")

netVisual_aggregate(normal_cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "chord")

#?netVisual_aggregate
#par(mfrow=c(1,1))

#netVisual_aggregate(normal_cellchat, signaling = pathways.show,vertex.receiver = vertex.receiver, layout = "circle")
netVisual_aggregate(normal_cellchat, signaling = pathways.show, layout = "circle")




netAnalysis_contribution(normal_cellchat, signaling = pathways.show)


####get specific  cell 
levels(normal_cellchat@idents)

##某个细胞亚群对其他信号通路有哪些通路
netVisual_bubble(normal_cellchat, sources.use = c(1), targets.use = c(2:10), remove.isolate = FALSE)###epithelial
netVisual_bubble(normal_cellchat, sources.use = c(2), targets.use = c(8,9), remove.isolate = FALSE)##TNK
netVisual_bubble(normal_cellchat, sources.use = c(5), targets.use = c(1:4,6:10), remove.isolate = FALSE)##myeloid cell
#saveRDS(normal_cellchat,file="./Cellchat/normalcellchat_analysis.rds")





###########进行实验组cellchat
###We have run these in the workstation###

tumor_cellchat<-readRDS(file="cellchat_tumor.rds")
tumor_cellchat@DB <- CellChatDB.mouse

#gc()
tumor_cellchat <- subsetData(tumor_cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel

tumor_cellchat <- identifyOverExpressedGenes(tumor_cellchat)
tumor_cellchat <- identifyOverExpressedInteractions(tumor_cellchat)
tumor_cellchat <- projectData(tumor_cellchat, PPI.mouse)
tumor_cellchat <- computeCommunProb(tumor_cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
tumor_cellchat <- filterCommunication(tumor_cellchat, min.cells = 1)
tumor_cellchat <- computeCommunProbPathway(tumor_cellchat)
tumor_cellchat <- netAnalysis_computeCentrality(tumor_cellchat, slot.name = "netP") 


saveRDS(tumor_cellchat,file="tumorcellchat_6h.rds")
tumor_cellchat<-readRDS("./Cellchat/tumorcellchat_1.rds")
tumor_cellchat <- aggregateNet(tumor_cellchat)
groupSize_IRI6h <- as.numeric(table(tumor_cellchat@idents))
groupSize_IRI6h
table(tumor_cellchat@idents)

###########另一个样本
tumor_cellchat_1d@DB <- CellChatDB.mouse

#gc()
tumor_cellchat_1d <- subsetData(tumor_cellchat_1d) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel

tumor_cellchat_1d <- identifyOverExpressedGenes(tumor_cellchat_1d)
tumor_cellchat_1d <- identifyOverExpressedInteractions(tumor_cellchat_1d)
tumor_cellchat_1d <- projectData(tumor_cellchat_1d, PPI.mouse)
tumor_cellchat_1d <- computeCommunProb(tumor_cellchat_1d, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
tumor_cellchat_1d <- filterCommunication(tumor_cellchat_1d, min.cells = 1)
tumor_cellchat_1d <- computeCommunProbPathway(tumor_cellchat_1d)
tumor_cellchat_1d <- netAnalysis_computeCentrality(tumor_cellchat_1d, slot.name = "netP") 

saveRDS(tumor_cellchat_1d,file="tumorcellchat_1d.rds")
tumor_cellchat_1d <- aggregateNet(tumor_cellchat_1d)
groupSize_1d <- as.numeric(table(tumor_cellchat_1d@idents))
groupSize_1d
table(tumor_cellchat@idents)


groupSize_control

par(mfrow = c(2,2), xpd=TRUE)




par(mfrow = c(2,3), xpd=TRUE)
netVisual_circle(normal_cellchat@net$count,arrow.size = 0.01, vertex.weight = groupSize_control, weight.scale = T, label.edge= F, title.name = "Number of interactions control")
netVisual_circle(tumor_cellchat@net$count,arrow.size = 0.01, vertex.weight = groupSize_IRI6h, weight.scale = T, label.edge= F, title.name = "Number of interactions IRI6h")
netVisual_circle(tumor_cellchat_1d@net$count,arrow.size = 0.01, vertex.weight = groupSize_1d, weight.scale = T, label.edge= F, title.name = "Number of interactions IRI1d")
netVisual_circle(normal_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize_control, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength control")
netVisual_circle(tumor_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize_IRI6h, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength IRI6h")
netVisual_circle(tumor_cellchat_1d@net$weight, arrow.size = 0.01,vertex.weight = groupSize_1d, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength IRI1d")



par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(tumor_cellchat_1d@net$count,arrow.size = 0.01, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(tumor_cellchat_1d@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



mat <- tumor_cellchat@net$count

##mat <- tumor_cellchat@net$weight
#dev.off()

####下面针对其中一个细胞亚型分析其对其他细胞的interaction的Ligand-receptor数量
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[1, ] <- mat[1, ]
netVisual_circle(mat2, vertex.weight = groupSize,arrow.size = 0.2, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[1])


###查看具体的信号通路
df.net <- subsetCommunication(tumor_cellchat)
#df.net_tumor <- subsetCommunication(tumor_cellchat)

normal_cellchat@netP$pathways
tumor_cellchat@netP$pathways
tumor_cellchat_1d@netP$pathways

levels(df.net$source)
pathways.show <- c("PDGF") 

#table(df.net$pathway_name)
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#dev.new()
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(tumor_cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")
?netVisual_aggregate

par(mfrow=c(1,1)) 
netVisual_aggregate(tumor_cellchat, signaling = pathways.show, layout = "circle")
netAnalysis_contribution(tumor_cellchat, signaling = pathways.show)


####get specific  cell 

netVisual_bubble(normal_cellchat, sources.use = c(1), targets.use = c(2:10), remove.isolate = FALSE)###epithelial
netVisual_bubble(tumor_cellchat, sources.use = c(1), targets.use = c(2:10), remove.isolate = FALSE)###epithelial
netVisual_bubble(normal_cellchat, sources.use = c(2), targets.use = c(1,3:10), remove.isolate = FALSE)##TNK
netVisual_bubble(tumor_cellchat, sources.use = c(2), targets.use = c(1,3:10), remove.isolate = FALSE)##TNK
netVisual_bubble(normal_cellchat, sources.use = c(5), targets.use = c(1:4,6:10), remove.isolate = FALSE)##myeloid cell
netVisual_bubble(tumor_cellchat, sources.use = c(5), targets.use = c(1:4,6:10), remove.isolate = FALSE)##myeloid cell


#netAnalysis_signalingRole_network(normal_cellchat, signaling = c("IGF"), width = 8, height = 2.5, font.size = 10)#,"VEGF"



# the slot 'netP' means the inferred intercellular communication network of signaling pathways

##########一种热图，某个通路下sender、receiver
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(tumor_cellchat, signaling = c("CXCL"), width = 8, height = 2.5, font.size = 10)#,"VEGF"
netAnalysis_signalingRole_network(tumor_cellchat, signaling = c("IFN-II"), width = 8, height = 2.5, font.size = 10)

ht1 <- netAnalysis_signalingRole_heatmap(tumor_cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(tumor_cellchat, pattern = "incoming")
ht1 + ht2


#normal_cellchat <- netAnalysis_computeCentrality(normal_cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(normal_cellchat, signaling = c("IGF"), width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(normal_cellchat, signaling = c("IFN-II"), width = 8, height = 2.5, font.size = 10)

ht1 <- netAnalysis_signalingRole_heatmap(normal_cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(normal_cellchat, pattern = "incoming")
ht1 + ht2



#saveRDS(normal_cellchat,file="./Cellchat/normalcellchat_analysis.rds")
#saveRDS(tumor_cellchat,file="./Cellchat/tumorcellchat_analysis.rds")

###we next compare the cell signaling pathways between the normal and tumor tissues##
data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)
object.list <- list(normal = normal_cellchat, tumor = tumor_cellchat,tumor_1d=tumor_cellchat_1d)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)###counts 差异
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")###weight的差异

gg1 <- netVisual_heatmap(cellchat)

#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2



library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.7.1.1010
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
########上面是一个组织类型中的incoming和outgoing比较，下面是两种组织的比较
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"),height=unit(12, "cm"),width=unit(16, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"),height=unit(12, "cm"),width=unit(16, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"),height=unit(12, "cm"))


########针对某种细胞进行分析
cellchat@idents
##Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(4,7),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2,3),  comparison = c(1, 2), angle.x = 45)
###指定信号通路
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(4,7), signaling = c("MIF"), comparison = c(1, 2), angle.x = 90)




#########探究某个信号通路下所有基因进行不同细胞群表达
colnames(mat2)
##compare specific gene clusters### vlnplot
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("normal", "tumor")) # set factor level

plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "VEGF", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "IGF", split.by = "datasets", colors.ggplot = T)

save.image(file="./Cell_chat.RData")



