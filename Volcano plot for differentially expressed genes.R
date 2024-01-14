
rm(list=ls())
devtools::install_github("BioSenior/ggvolcano")
library(ggVolcano)

## basic example code
setwd("")     #???ù???Ŀ¼
deg_data=read.table("GSE43974_Diff_火山图_input.txt", header=T, sep="\t", check.names=F)


# use the function -- add_regulate to add a regulate column 
# to the DEG result data. 
data <- add_regulate(deg_data, log2FC_name = "logFC",
                     fdr_name = "adj.P.Val",log2FC = 0.5, fdr = 0.05)


# plot
ggvolcano(data, x = "log2FoldChange", y = "padj",log2FC_cut = 0.5,FDR_cut = 0.05,
          colors = c("#00468BFF","#999999", "#AD002AFF"),fills = c("#00468BFF","#999999", "#AD002AFF"),
          label = "id",label_number = 0,output = FALSE,custom_label = NULL)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave("GSE43974_Diff_火山图.pdf",height=5,width=7)


