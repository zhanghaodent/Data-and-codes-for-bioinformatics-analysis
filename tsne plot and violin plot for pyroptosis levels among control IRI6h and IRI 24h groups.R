

###########irGSEA分析自定义数据组间集富集分析
rm(list=ls())
setwd("")
sce.mergeTEN<-readRDS(file="mergeTEN_2_after_anno.rds")
library(UCell)
library(irGSEA)
library(ggpubr)
gene <- read.table("细胞焦亡.txt",header=T,sep="\t",row.names = 1)
gene_1 <- list()
gene_1$TcellN <- gene$gene
sce.final <- irGSEA.score(object = sce.mergeTEN, assay = "RNA", slot = "data", 
                          seeds = 123, ncores = 1,msigdb=F, custom = T, geneset = gene_1, 
                          method = c("AUCell", "UCell", "singscore",
                                     "ssgsea"), kcdf = 'Gaussian')
sce.final@meta.data[1:5,]
##########富集评分的密度tsne或umap图
scatterplot <- irGSEA.density.scatterplot(object = sce.final,
                                          method = "singscore",
                                          show.geneset = c("TcellN"),
                                          reduction = "tsne")
scatterplot 

ridgeplot <- irGSEA.ridgeplot(object = sce.final,
                              method = "AUCell",
                              show.geneset = c("TcellN"))
ridgeplot

#########计算每群细胞富集评分高低的热图
densityheatmap <- irGSEA.densityheatmap(object = sce.final,
                                        method = "UCell",
                                        show.geneset = c("TcellN"))

#########计算每群细胞富集评分的高低
halfvlnplot <- irGSEA.halfvlnplot(object = sce.final,
                                  method = "UCell",
                                  show.geneset = c("TcellN"))+stat_compare_means()
halfvlnplot

######取出富集分数构建数据框
df <- as.data.frame(sce.final@meta.data)
colnames(df)

########只选择其中特定细胞进行绘图分析
df_new <- subset(df,Seurat_harmony=="DCT"|Seurat_harmony=="DC-PC")

############绘制分组各个细胞群的富集分数箱线图
p=ggboxplot(df_new, x="Seurat_harmony", y="nCount_singscore", color = "orig.ident", 
            ylab="Gene expression",
            xlab="",
            legend.title="",
            palette = c("#00468BFF","#E69F00","#AD002AFF"),
            width=0.6, add = "mean_sd")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=orig.ident),
                        method="kruskal.test",
                        symnum.args=list(cutpoints = c(0,0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", " ")),
                        label = "p.signif")
p2 <- p1 + ggprism::theme_prism(palette = "black_and_white",
                                base_size = 12, 
                                #base_line_size = 1,#base_family = "Roboto Condensed",
                                axis_text_angle = 0,
                                border = TRUE,base_rect_size = 1.2)
p2

#############绘制各个细胞群的富集分数小提琴图
d <- ggviolin(df_new, x="Seurat_harmony", y="nCount_singscore", 
              width =0.8, #设置小提琴宽度
              fill="orig.ident",#填充
              palette =c("#00468BFF","#E69F00","#AD002AFF"),#设置颜色
              add = 'mean_sd',#添加均值和标准差
              xlab = F, #不显示x轴的标签
              legend = "right")#图例显示在右侧)
d1 <- d+stat_compare_means(aes(group=orig.ident),
                           method="kruskal.test",
                           symnum.args=list(cutpoints = c(0,0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", " ")),
                           label = "p.signif")
d2 <- d1 + ggprism::theme_prism(palette = "black_and_white",
                                base_size = 12, 
                                #base_line_size = 1,#base_family = "Roboto Condensed",
                                axis_text_angle = 0,
                                border = TRUE,base_rect_size = 1.2)
d2

