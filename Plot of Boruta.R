#这里引用earth包的内置数据集“ozone1”作示例，节选自洛杉矶的臭氧观测数据，并包含部分其它环境因子。

#洛杉矶的臭氧数据集，详情加载 earth 包后 ?ozone1
library(earth)
data(ozone1)
head(ozone1)
#该数据集共由10列变量的330行观测值组成
rm(list=ls())
library(Boruta)
#详情 ?Boruta
#pValue 指定置信水平，mcAdj=TRUE 意为将应用 Bonferroni 方法校正 p 值
#此外还可以提供 mtry 和 ntree 参数的值，这些值将传递给随机森林函数 randomForest()
setwd("")
rt=read.table("GSE53605_merge_fustat.txt",header=T,sep="\t",row.names=1)   
gene=read.table("target_gene.txt",header=T)
rt=rt[,c("fustat",as.vector(gene[,1]))]

set.seed(123)
boruta_output <- Boruta(fustat ~ ., data=rt, pValue = 0.05, mcAdj = TRUE, doTrace = 2)
summary(boruta_output)

importance <- boruta_output$ImpHistory  #该结果项可以理解为 Boruta 获得的变量重要性得分
importance
#write.csv(importance, 'Variable_Importance.csv', row.names = FALSE)
pdf("boruta_output.pdf")
plot(boruta_output, las = 2, xlab = '', main = 'Variable Importance')
dev.off()


#在这个示例中，默认执行了14次迭代运算，对于每个变量，在每次迭代中各获得一个重要性得分。
#最后通过比较这些变量的重要性得分的高低，即可评估哪些环境环境变量对于预测或者解释臭氧浓度是更重要的。
#例如根据排名可知，前5个重要的变量依次为:
#temp（温度，℉）> ibt（逆温基底温度，℉）> ibh（逆温基底高度，英尺）> dpg（压力梯度，mm Hg）> humidity（湿度百分比）

#查看Boruta 获得的变量重要性得分的统计（均值、中位数、最大值、最小值等）
#此外如若构建模型，Boruta 还可帮助评估哪些变量被建议在模型中保留（Confirmed）还是剔除（Rejected）
boruta_df <- attStats(boruta_output)
boruta_df
write.table(boruta_df,file="geneCoef_gene.xls",sep="\t",quote=F,row.names=T)





