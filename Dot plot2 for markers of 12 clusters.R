#Dotplot_anno本函数的主要作用是做单细胞基因表达气泡图或者普通数据气泡图，可添加注释
#参数解释
#object：单细胞seurat对象或者普通data frame。是seurat对象的时候，参数single_type=T
#single_type=T：数据类型判断，如果是单细胞数据这里必须是T，否则选择False
#features，单细胞数据的时候为基因向量，否则为NULL
#group：分组注释向量，这个视基因组别而定
#color：表达量密度连续色彩，类似于普通热图颜色注释。
#order：单细胞数据情况下（需要设置T，不需要F）,普通数据情况下使用自行选择，是否按照表达量从大到小排列。类似于热图每个细胞类型高表达marker从大到小排列展示
#celltype_color：分组注释颜色
#level：针对普通数据，在order=T情况下使用，自定义y轴顺序
#legend.postion:legend的位置
#heatmap:针对普通表达数据，是否使用热图形式展示
#bulk_feature:非单细胞数据中基因或者GO terms的列名
#bulk_samples：非单细胞数据中样本列名
Dotplot_anno <- function(object,
                         single_type=T,
                         features=NULL,
                         bulk_feature=NULL,
                         bulk_samples=NULL,
                         group,
                         color,
                         order,
                         celltype_color,
                         level,
                         legend.postion=NULL,
                         heatmap=T){
  if (single_type==T){
    
    data.usage <- DotPlot(object,features = markers)$data
    data.anno <- data.frame(
      features.plot = unique(data.usage$features.plot),
      label = group)
    
    df.plot <- plyr::join(data.usage,data.anno)
    
    if(order==T){
      df.plot$id <- factor(df.plot$id,levels = sort(levels(df.plot$id),decreasing = T))
      
      p <- ggplot(df.plot,aes(x=features.plot,
                              y =  as.numeric(id),
                              size = pct.exp, 
                              color = avg.exp.scaled))+
        geom_point() + 
        scale_size(range = c(0,10)) + 
        scale_color_gradientn(colours = color,
                              guide = guide_colorbar(ticks.colour = "black",
                                                     frame.colour = "black"),
                              name = "Average\nexpression") +
        ylab("") + xlab("") +
        scale_y_continuous(breaks = 1:length(levels(df.plot$id)),
                           labels = levels(df.plot$id),
                           sec.axis = dup_axis())+ 
        facet_grid(~label, scales="free_x",space = "free")+
        theme_classic() +
        theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
              
              axis.ticks.y = element_blank(),
              axis.text.y.right = element_blank(),
              axis.ticks.x = element_blank(),
              
              
              panel.spacing=unit(0, "mm"), 
              strip.text.x = element_text(size=12,color = "black",
                                          vjust = 0.5,margin = margin(b = 3,t=3)),
              strip.background = element_rect(colour="black",size = 1),
              plot.margin=unit(c(1, 1, 1, 1),'cm'),
              legend.position = legend.postion)
      
      
    }else{
      
      p <- ggplot(df.plot,aes(x=features.plot,
                              y =  as.numeric(id),
                              size = pct.exp, 
                              color = avg.exp.scaled))+
        geom_point() + 
        scale_size(range = c(0,100)) + 
        scale_color_gradientn(colours = color,
                              guide = guide_colorbar(ticks.colour = "black",
                                                     frame.colour = "black"),
                              name = "Average\nexpression") +
        ylab("") + xlab("") +
        scale_y_continuous(breaks = 1:length(levels(df.plot$id)),
                           labels = levels(df.plot$id),
                           sec.axis = dup_axis())+ 
        facet_grid(~label, scales="free_x",space = "free")+
        theme_classic() +
        theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
              
              axis.ticks.y = element_blank(),
              axis.text.y.right = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_line(colour = 'black',linewidth = 0.5), 
              
              panel.spacing=unit(0, "mm"), 
              strip.text.x = element_text(size=12, color = "black",
                                          vjust = 0.5,margin = margin(b = 3,t=3)),
              strip.background = element_rect(colour="black",size = 1),
              plot.margin=unit(c(1, 1, 1, 1),'cm'),
              legend.position = legend.postion)
      
    }
  }else{
    
    data.usage <- object
    data.anno <- data.frame(
      features.plot = unique(data.usage[, colnames(data.usage)%in% bulk_feature]),
      label = group)
    
    df.plot <- plyr::join(data.usage,data.anno)
    df.plot[, colnames(df.plot)%in% bulk_samples]<- as.factor(df.plot[, colnames(df.plot)%in% bulk_samples])
    
    if(order==T){
      df.plot[, colnames(df.plot)%in% bulk_samples] <- factor(df.plot[, colnames(df.plot)%in% bulk_samples],levels = level)
    }else{
      df.plot = df.plot
    }
    
    
    if('ratio' %in% colnames(df.plot)){
      p <- ggplot(df.plot,aes(x=features.plot,
                              y =  id,
                              size=ratio,
                              color = exp))+
        geom_point() + 
        scale_size(range = c(0,5))+
        scale_color_gradientn(colours = color,
                              guide = guide_colorbar(ticks.colour = "black",
                                                     frame.colour = "black"),
                              name = "-Log10(P)") +
        ylab("") + xlab("") +
        facet_grid(~label, scales="free_x",space = "free")+
        theme_classic() +
        theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
              
              axis.ticks.y = element_blank(),
              axis.text.y.right = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_line(colour = 'black',linewidth = 0.5), 
              
              panel.spacing=unit(0, "mm"), 
              strip.text.x = element_text(size=12, color = "black",
                                          vjust = 0.5,margin = margin(b = 3,t=3)),
              strip.background = element_rect(colour="black",size = 1),
              plot.margin=unit(c(1, 1, 1, 1),'cm'),
              legend.margin = margin(-0.2,-0.2,0,0,'cm'),
              legend.key.height = unit(0.4,'cm'),
              legend.key.width = unit(0.4,'cm'),
              legend.position = legend.postion)+
        guides(size=guide_legend(title="Ratio"))
    }else{
      
      if (heatmap==T){
      
      p <- ggplot(df.plot,aes(x=features.plot,
                              y =  id))+
        geom_tile(aes(fill = exp), colour = "black", size = 0.5) + 
        scale_fill_gradientn(colours = color,
                             guide = guide_colorbar(ticks.colour = "black",
                                                    frame.colour = "black"),
                             name = "Relative\nexpression") +
        ylab("") + xlab("") +
        facet_grid(~label, scales="free_x",space = "free")+
        theme_classic()+
        theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
              
              axis.ticks.y = element_blank(),
              axis.text.y.right = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_blank(), 
              
              panel.spacing=unit(0, "mm"), 
              strip.text.x = element_text(size=12, color = "black",
                                          vjust = 0.5,margin = margin(b = 2,t=2)),
              strip.background = element_rect(colour="black",size = 0.5),
              plot.margin=unit(c(1, 1, 1, 1),'cm'),
              legend.position = legend.postion)+
        scale_x_discrete(expand = c(0,0))+
        scale_y_discrete(expand = c(0,0))
    }else{
      
      p <- ggplot(df.plot,aes(x=features.plot,
                              y =  id,
                              size=5,
                              color = exp))+
        geom_point() + 
        guides(size="none")+
        scale_color_gradientn(colours = color,
                              guide = guide_colorbar(ticks.colour = "black",
                                                     frame.colour = "black"),
                              name = "Relative\nexpression") +
        ylab("") + xlab("") +
        facet_grid(~label, scales="free_x",space = "free")+
        theme_classic() +
        theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
              
              axis.ticks.y = element_blank(),
              axis.text.y.right = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_line(colour = 'black',linewidth = 0.5), 
              
              panel.spacing=unit(0, "mm"), 
              strip.text.x = element_text(size=12, color = "black",
                                          vjust = 0.5,margin = margin(b = 3,t=3)),
              strip.background = element_rect(colour="black",size = 1),
              plot.margin=unit(c(1, 1, 1, 1),'cm'),
              legend.position = legend.postion)
    }
      
    
    }

    
  }
  
  
  g <- ggplot_gtable(ggplot_build(p))
  strips <- which(grepl('strip-', g$layout$name))
  fills <- celltype_color
  k = 1
  
  for (i in strips) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot(g)
  
}

