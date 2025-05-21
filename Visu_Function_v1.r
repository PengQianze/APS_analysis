###########函数文件
###########function field

STACK_ANA <- function (taxnonmy,metadata,analysisgroup, label = "ana_", level = level,
                       colorsF = c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                                   "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                                   "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39") ,
                       topN = 10 
){
  library(reshape2)
  library(ggplot2)
  library(tibble)
  library( rstatix)
  idx = rownames(metadata) %in% colnames(taxonomy)
  metadata = metadata[idx, , drop = F]
  taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  colnames(sampFile)[1] = "Group"
  
  mean_sort = as.data.frame(taxonomy[(order(-rowSums(taxonomy))), ])
  idx = grepl("unassigned|unclassified|unknown", rownames(mean_sort), ignore.case = T)
  mean_sort = rbind(mean_sort[!idx, ], mean_sort[idx, ])
  other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
  mean_sort = mean_sort[1:(topN - 1), ]
  mean_sort = rbind(mean_sort, other)
  rownames(mean_sort)[topN] = c("Other")
  colnames(sampFile)[1] = "Group"
  amplic2 = t(mean_sort)
  
  mat_t2 = merge(sampFile, amplic2, by = "row.names")
  mat_t2 = mat_t2[, c(-1)]
  mat_mean = aggregate(mat_t2[, -1], by = mat_t2[1], FUN = mean)
  mat_mean_final = do.call(rbind, mat_mean)[-1, ]
  colnames(mat_mean_final) = mat_mean$Group
  mean_sort2 = as.data.frame(mat_mean_final)
  mean_sort2$Taxonomy = rownames(mean_sort2)
  data_all = as.data.frame(melt(mean_sort2, id.vars = c("Taxonomy")))
  data_all$Taxonomy = factor(data_all$Taxonomy, levels = rownames(mean_sort2))
  data_all$value = as.numeric(data_all$value)
  
  #mat_mean$Group <- reorder(mat_mean$Group, -mat_mean$Pseudomonadales)
 #data_all <- data_all %>%mutate(variable=factor(variable,levels = levels(mat_mean$Group)))
  
  data_all <- data_all %>%mutate(variable=factor(variable,levels =level))
  
  p7 = ggplot(data_all, aes(x = variable, y = value, fill = Taxonomy)) + 
#   geom_bar(stat = "identity", position = "fill", width = 0.8,colour ="white") + 
    geom_bar(stat = "identity", position = "stack", width = 0.7,colour ="white") + 
#    scale_y_continuous(labels = scales::percent) +
    xlab("") + ylab("Relative abundance (%)") +
    scale_fill_manual(values = colorsF)+
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
          strip.text = element_text(size=12),
          axis.line = element_line(color = "black",linewidth = 0.25),
          axis.ticks =  element_line(color = "black",linewidth = 0.25),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(color="black",size=10),
          axis.text.x = element_text(color="black",size=10,angle = 0),
          panel.spacing.x = unit(0,"cm"),
          panel.border = element_blank(),
          panel.spacing = unit(0,"lines")
    )
  
#  pcxlong<-length(level)*2
#  pcylong<-topN*0.6
  ggsave(paste0(label,".",analysisgroup,"stack.pdf"), p7, width=8, height=6, units="in")
  ggsave(paste0(label,".",analysisgroup,"stack.png"), p7, width=8, height=6, units="in")
  return(p7)
}


PCOA_ANA <- function (taxonomy,metadata,analysisgroup = "TimeGrowth",
                      colors= c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                   "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                              "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39"),
                      shapes = c(21,22,23,24,25),label="ana_"
){
  
  library(vegan)
  library(ggplot2)
  
  idx = rownames(metadata) %in% colnames(taxonomy)
  metadata = metadata[idx, , drop = F]
  taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  colnames(sampFile)[1] = "Group"
  
  data <- vegdist(t(taxonomy),method = "bray")
  write.table(as.matrix(data), 'bray.txt', sep = '\t', col.names = NA, quote = FALSE)#输出距离矩阵
  otu_matrix<-read.table("bray.txt", header=T, row.names=1, sep="\t", comment.char="")
  colnames(sampFile)[1] = "group"
  
  #非限制性排序及组间显著性分析
  #sampFile = as.data.frame(metadata[, "group"], row.names = row.names(metadata))
  pcoa = cmdscale(otu_matrix, k = 3, eig = T)
  points = as.data.frame(pcoa$points)
  eig = pcoa$eig
  points = cbind(points, sampFile)
  colnames(points) = c("x", "y", "z", "group")
  
  #显著性
  adonis_table = adonis2(otu_matrix ~ group, data = sampFile, permutations = 10000)
  adonis_pvalue = adonis_table$`Pr(>F)`[1]
  adonis_R2=adonis_table$`R2`[1]
  adonis_R2pvalue = paste("group","R2:", adonis_R2, "pvalue:", adonis_pvalue, sep = "\t")
  write.table(adonis_R2pvalue, file = paste0(label,".",analysisgroup,"-adonis2test.txt"), append = TRUE, 
              sep = "\t", quote = F, row.names = F, col.names = F)
  
  pcoa.fit <- envfit(pcoa ~ group, data=sampFile, perm=999)
  envfitR2 = paste("group R2:", pcoa.fit$factors$r, "Pvalue:", pcoa.fit$factors$pvals,sep = "\t")
  write.table(envfitR2,paste0(label,".",analysisgroup,"-envfit.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  
  #可视化
  p2 = ggplot(points, aes(x = x, y = y, color = group,shape = group,fill = group)) + 
    labs(x = paste("PCo 1 (", format(100 * eig[1]/sum(eig), digits = 4), "%)", sep = ""), 
         y = paste("PCo 2 (", format(100 * eig[2]/sum(eig), digits = 4), "%)", sep = ""), color = "group")+
    scale_fill_manual(values = colors)+
    scale_colour_manual(values = colors)+
    scale_shape_manual(values = shapes)
  p2 = p2+ geom_point(alpha = 0.8, size = 4,color = "black") + 
    theme(
      strip.text = element_text(size=12),
      axis.line = element_line(color = "black",linewidth = 0.25),
      axis.ticks= element_line(colour = "black",size = 0.25),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.y = element_text(color="black",size=10),
      axis.text.x = element_text(color="black",size=10,angle = 0),
      panel.spacing.x = unit(0,"cm"),
      panel.border = element_blank(),
      panel.spacing = unit(0,"lines")
    )+
    stat_ellipse(level = 0.68)+ 
    ggtitle(paste( " P = ", pcoa.fit$factors$pvals,  " R2 = ", pcoa.fit$factors$r, sep = ""))
  p2

  
  ggsave(paste0(label,".",analysisgroup,"-PCOA12.pdf"), p2, width=4, height=4, units="in")
  ggsave(paste0(label,".",analysisgroup,"-PCOA12.png"), p2, width=4, height=4, units="in")
  
  return(p2)
}


CPCOA_ANA <- function (taxonomy,metadata,analysisgroup = "TimeGrowth",
                       colors= c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                                "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                                "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39"),
                        shapes = c(21,22,23,24,25),label="ana_"
){
  
  library(vegan)
  library(ggplot2)
  
  idx = rownames(metadata) %in% colnames(taxonomy)
  metadata = metadata[idx, , drop = F]
  taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"
  
  #CPCOA
  otu_cap <- capscale(t(taxonomy) ~ group, data = sampFile, add = F, sqrt.dist = T, distance = "bray")
  variance = (otu_cap$CCA$tot.chi/otu_cap$tot.chi)
  eig = otu_cap$CCA$eig
  points = as.data.frame(otu_cap$CCA$wa)
  points = cbind(sampFile, points[rownames(points), ])
  
  #显著性分析及R2
  otu_cap_envfit <- envfit(otu_cap,metadata,permu=1000)
  cpcoaenvfitR2 = paste("CAPgroup R2:", otu_cap_envfit$factors$r, "Pvalue:", otu_cap_envfit$factors$pvals,
                        "variance:",variance,sep = "\t")
  write.table(cpcoaenvfitR2,paste0(label,".",analysisgroup,"-CCAenvfit.txt"),append = TRUE, sep = "\t", quote = F, row.names = F, col.names = F)
  perm_otu_cap = anova.cca(otu_cap, permutations = 1000, parallel = 4)
  p.val = perm_otu_cap[1, 4]
  
  
  #可视化
  p3 = ggplot(points, aes(x = CAP1, y = CAP2, color = group,shape= group,fill=group)) + 
    geom_point(alpha = 0.8, size = 4,color = "black") + 
    labs(x = paste("CPCo 1 (", format(100 * eig[1]/sum(eig), digits = 4), "%)", sep = ""), 
         y = paste("CPCo 2 (", format(100 * eig[2]/sum(eig), digits = 4), "%)", sep = ""), color = "group") + 
    ggtitle(paste(format(100 * variance, digits = 3), " % of variance; P = ", format(p.val, digits = 2), sep = "")) +
    theme_classic() + 
    theme(
      strip.text = element_text(size=12),
      axis.line = element_line(color = "black",linewidth = 0.25),
      axis.ticks= element_line(colour = "black",size = 0.25),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.y = element_text(color="black",size=10),
      axis.text.x = element_text(color="black",size=10,angle = 0),
      panel.spacing.x = unit(0,"cm"),
      panel.border = element_blank(),
      panel.spacing = unit(0,"lines")
    )+
    scale_fill_manual(values = colors)+
    scale_colour_manual(values = colors )+
    scale_shape_manual(values = shapes)
#  p3 = p3 + stat_ellipse(level = 0.68)
  
  pcxlong<-length(level)*1.8
  pcylong<-length(level)*1.2
  ggsave(paste0(label,".",analysisgroup,"-CPCOA12.pdf"), p3, width=6, height=4, units="in")
  ggsave(paste0(label,".",analysisgroup,"-CPCOA12.png"), p3, width=6, height=4, units="in")
  
  return(p3)
}

box_ANA <- function (alphadiversity,metadata,analysisgroup = "TimeGrowth", col_name = "HCA",
                     colors = c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                                "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96"),level,shapes,label="ana_"
){
  library(rstatix)
  library(ggplot2)
  library( ggpubr )
  idx = rownames(metadata) %in% rownames(alphadiversity)
  metadata1 = metadata[idx, , drop = F]
  alpha = alphadiversity[rownames(metadata1),]
  sampFile2 = as.data.frame(metadata1[, c(analysisgroup)], row.names = row.names(metadata1))
  colnames(sampFile2)[1] = "Group"
  alpha_dif = cbind(sampFile2, alpha)
  df <- alpha_dif%>%mutate(Group=factor(Group,levels =level))
  
  dfn <- df[c(col_name, "Group")]
  names(dfn) <- c("Content", "Group")
  stat.test<- dfn %>% t_test(Content ~ Group)%>%
  add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")


  
  theme<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
               strip.text = element_text(size=12),
               axis.line = element_line(color = "black",linewidth = 0.25),
               axis.ticks =  element_line(color = "black",linewidth = 0.25),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(color="black",size=10),
               axis.text.x = element_text(color="black",size=10,angle = 0),
               panel.spacing.x = unit(0,"cm"),
               panel.border = element_blank(),
               panel.spacing = unit(0,"lines")
  )
    
  p1 <-ggplot(dfn,aes(Group,Content)) +
    geom_boxplot(aes(fill=Group,alpha = 0.8),alpha=1, width = 0.4,outlier.shape = NA, outlier.size = 0)+
#    geom_point(aes(fill=Group,alpha=0.7),colour ="black",alpha=0.9,size =1.5, shape = 21, stroke = 0.5,position = position_jitter(0.1))+
    
#    geom_boxplot(aes(alpha = 1,colour=Group),alpha=1, width = 0.7,outlier.shape = NA, outlier.size = 0,fill= "transparent")+
#    geom_point(aes(fill=Group,alpha=0.7,colour =Group),alpha=0.8,size =0.8, shape = 21, stroke = 0.5,position = position_jitter(0.25))+
    scale_fill_manual(values = colors)+
#    scale_colour_manual(values = colors)+
    stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
    scale_size_continuous(range=c(1,3)) +  
    ylab(col_name)+
    ggtitle(paste(col_name ))+theme
  p1
  
  ggsave(paste0(label,".",analysisgroup,col_name,".pdf"), p1, width=6, height=6, units="in")
  ggsave(paste0(label,".",analysisgroup,col_name,".png"), p1, width=6, height=6, units="in")

  return(p1)
} 


volcano_plot <- function (data1, cut_off_FDR, cut_off_log2FC, COL = c("FC","Pvalue","name"),
                          color,label="ana_",topn=5
                          
) {
  library(ggplot2)
  library(ggrepel)
  
  data = as.data.frame(data1[, COL], row.names = row.names(data1))
  colnames(data)<- c("v1", "v2" ,"FDR","gene_name")
  data$FC <- data$v1/data$v2
  data$log2FC <- log2(data$FC)
  data$log10FDR <- log10(data$FDR)
  # 设定差异基因阈值
  
  data$Sig = ifelse(data$FDR < cut_off_FDR & abs(data$log2FC) >= cut_off_log2FC,
                    ifelse(data$log2FC > cut_off_log2FC, 'Up', 'Down'), 'no')
  
  producedata = as.data.frame(data[, !(names(data) %in% c("FDR","gene_name"))], row.names = row.names(data))
  newdata <- cbind(data1,producedata)
  
  write.csv(as.data.frame(newdata), file =paste0(label,"core_otu.csv") , row.names = FALSE)
  # 绘制火山图
  p1 <- ggplot(data, aes(x = log2FC, y = -log10FDR, colour = Sig)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = color) +
    #  xlim(c(-30, 30)) +
    geom_vline(xintercept = c(-cut_off_log2FC, cut_off_log2FC), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(cut_off_FDR), lty = 4, col = "black", lwd = 0.8) +
    labs(x = "log2FC", y = "-log10FDR") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank())
  p1

  up <- data[data$Sig  == 'Up', ]
  top_5_up <- order(up$log10FDR, decreasing = FALSE)[1:topn]
  top_5_rows_up <- up[top_5_up, ]
  
  down <- data[data$Sig  == 'Down', ]
  top_5_down <- order(down$log10FDR, decreasing = FALSE)[1:topn]
  top_5_rows_down <- down[top_5_down, ]
  
  
  p2 <- p1 + geom_text_repel(
    data = top_5_rows_up,
    aes(label = gene_name), size = 1,colour = "black",
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
  p2 <-p2 + geom_text_repel(
    data = top_5_rows_down,
    aes(label = gene_name), size = 1,colour = "black",
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
  
  p2
  
  
  ggsave(paste0(label,".","volcano.pdf"), p1, width=5, height=4, units="in")
  ggsave(paste0(label,".","volcano.png"), p1, width=5, height=4, units="in")
  
  ggsave(paste0(label,"lable.","volcano.pdf"), p2, width=5, height=4, units="in")
  ggsave(paste0(label,"lable.","volcano.png"), p2, width=5, height=4, units="in")
  
  return(p2)
}



VEEN_ANA <- function (taxonomy,metadata,analysisgroup = "ORG"
){

library( reshape2 )
library( rstatix  )
library( plyr   )
sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
colnames(sampFile)[1] = "Group"

dat1 <- as.data.frame(base::merge(t(taxonomy),sampFile,by='row.names'))
df1 <-dat1[,-1]
dfnorder = as.data.frame(melt(df1, id.vars = c("Group")))
dfnorder[,3][dfnorder[,3] > 0] = 1

dfnorder = dfnorder[dfnorder$value == 1, ]
un <- unique(dfnorder, incomparables = FALSE)
wider <- spread(un, Group, value)

for (i in 1:nrow(wider)) {
  otu <- as.vector(wider[i, 1]) 
  wider[i, ][!is.na(wider[i, ])] <- otu
}

veen <- wider[, !names(wider) %in% "variable"]
write.csv(as.data.frame(veen ),paste0(analysisgroup,"-venn.csv"))

###############################取veen交补集

core <- Reduce(intersect, veen)
veen_out <- list(core = core)

for (col in colnames(veen)) {
    colindex <- which(names(veen) == {{col}})
    df <- veen[, c(colindex, setdiff(1:ncol(veen), colindex))]
    special<-Reduce(setdiff, df)
#    names( special)<-{{col}}
    veen_out <- append(veen_out, list(special), length(veen_out) + 1)
    names(veen_out)[colindex+1] <- {{col}}
}
veen_out <- do.call(cbind, lapply(lapply(veen_out, unlist), `length<-`, max(lengths(veen_out))))
write.csv(as.data.frame(veen_out),paste0(analysisgroup,"-veen_out.csv"))

########################计算共有，特有物种差异

idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(metadata)]

sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
colnames(sampFile)[1] = "Group"

dat1 <- as.data.frame(merge(t(taxonomy),sampFile,by='row.names'))
df1 <-dat1[,-1]

difana = as.data.frame(melt(df1, id.vars = c("Group")))

for (class in colnames(veen_out)){
    analysis_tax <- veen_out[,{{class}}]
    difana1 <- difana[difana$variable %in% analysis_tax, ]
    orderstat.test<- difana1 %>% group_by(variable)%>% t_test(value ~ Group)%>%
    add_significance(p.col = "p")
    write.csv(as.data.frame(orderstat.test),paste0(class,analysisgroup,"-t-test.csv"))
    sumorder<-ddply(dfnorder,c("Group","variable"),summarise,allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
                            n1=sum(!is.na(value)),se1=sd1/sqrt(n1))
    write.csv(as.data.frame(sumorder),paste0(class,analysisgroup,"-sum_.csv"))  
}

###############################venn图
library(ggVennDiagram)
p1 <- ggVennDiagram(veen,label = "count", label_alpha = 0,edge_lty = "dashed", edge_size = 0)

library(ggplot2)
p2 <-p1 + 
scale_fill_distiller(palette = "Set2", direction = -1) +
scale_color_brewer(palette = "Set2")

ggsave(paste0(analysisgroup,"-venn",".pdf"), p2,width=5, height=5, units="in")
    return(core)
}


pie_plot <- function (analysis, COL, topN,
                      colorsw = c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                               "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                               "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39"),
                      label="ana_"
) {
  library(ggplot2)
  library(ggrepel)
  library(scales)
  data = as.data.frame(analysis[ ,COL ], row.names = row.names(analysis))
  colnames(data) <- c("Order","allmean")
  result <- aggregate(allmean ~ Order, data = data , FUN = sum)
  mean_sort <- result[order(result$allmean, decreasing = TRUE), ]
  
  idx = grepl("Unassigned|unassigned|unclassified|unknown", rownames(mean_sort), ignore.case = T)
  mean_sort = rbind(mean_sort[!idx, ], mean_sort[idx, ])
  other = sum(mean_sort$allmean[topN:dim(mean_sort)[1]])
  other = data.frame(Order = "Other", allmean = other)
  mean_sort = mean_sort[1:(topN - 1), ]
  mean_sort = rbind(mean_sort, other)
  rownames(mean_sort)[topN] = c("Other")
  
  #mean_sort$Order = factor(mean_sort$Order, levels = mean_sort$Order )
  mean_sort$Order = factor(mean_sort$Order, levels = rev(mean_sort$Order) )
  # 将数据框转换为长格式，因为ggplot2更适合处理长格式数据
  
  # 绘制饼图
  p1 <-
    ggplot(mean_sort, aes(x = "", y = allmean, fill = Order)) +
    geom_bar(width = 0.5, stat = "identity",colour ="white") +
    scale_fill_manual(values = colorsw) +
    ##geom_text(aes(y = sum(allmean) - c(cumsum(allmean)[-length(allmean)],0), #+ allmean/10,
    geom_text(aes(y = c(0,cumsum(allmean)[-length(allmean)])+ allmean/2,
                  label = percent(allmean/sum(allmean))), size=3)+
    coord_polar(theta = "y", start=0) + # 将条形图转换为极坐标，形成饼图
    theme_void() + # 移除多余的坐标轴和标签
    labs(fill = "Category") # 设置图例标题
  
  ggsave(paste0(label,".","pie.pdf"), p1, width=5, height=4, units="in")
  ggsave(paste0(label,".","pie.png"), p1, width=5, height=4, units="in")
  
  return(p1)
}


bb_plot <- function (data1, COL, Topn,shapes=c(21,22,23,24),
                     colors=c("#00CD66","#7FFF00","#C0FF3E","#CAFF70","#8B864E","#CDCDB4","#CDCD00",
                              "#FFD700","#CD950C","#CD9B9B","#FF8247","#CDBA96","#FF00FF","#EE82EE","#FFA500",
                              "#009ACD","#AB82FF","#00FFFF","#8B668B","#CD4F39"),
                     label="ana_"
) {
  
  library(ggplot2)
  library(tibble)  
  #names(data)[names(data) == analysisgroup] <- "Order1"
  
  data = as.data.frame(data1[ ,COL ], row.names = row.names(data1))
  colnames(data) <- c("Order","MeanDecrease","shape1")
  data  <- rownames_to_column(data , var = "row")
  data <- data[order(data$MeanDecrease, decreasing = TRUE), ]
  top_MeanDecAcc <- head(data, Topn)
  
  p1<-ggplot(top_MeanDecAcc,aes(x = MeanDecrease , y = reorder(row,MeanDecrease)) )+
    #geom_col(aes(x = MeanDecrease , y = reorder(row,MeanDecrease) , fill =  Order) )+ 
    geom_segment(aes(x=0,xend=MeanDecrease,y=row,yend=row,color=Order,),linewidth=0.3,linetype="solid")+
    geom_point(aes(color=Order,fill=Order,shape=shape1),size=5)+
    scale_fill_manual(values = colors)+
    scale_colour_manual(values = colors)+ 
    scale_shape_manual(values = shapes)+
    xlab("MeanDecrease")+ylab("")+theme(
      axis.line.x = element_line(linewidth = 0.25, colour = "black"),
      axis.line.y = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black",angle=0))
  
  ggsave(paste0(label,".","bbplot.pdf"), p1, width=4.5, height=5, units="in")
  ggsave(paste0(label,".","bbplot.png"), p1, width=4.5, height=5, units="in")
  
  return(p1)
}

###########################################差异分析

DF_ANA <- function (taxonomy,metadata,analysisgroup = "TimeGrowth",repea=6,label="ana_"
){
  library(dplyr )
  library(rstatix )
  library(tidyr)
  
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  colnames(sampFile)[1] = "Group"
  taxonomy = taxonomy[rowSums(taxonomy)!=0, rownames(sampFile)]
  
  uniqlegth <- function (x){
    y = length(unique(x))
    return (y)
  }
  
  uniq.num.df  <- data.frame(num = apply(taxonomy, 1, uniqlegth))
  uniq.num.df.2 <- filter(uniq.num.df, num > repea*3)
  taxonomy2 <- taxonomy[rownames(uniq.num.df.2), ]
  dat1 <- as.data.frame(merge(t(taxonomy2),sampFile,by='row.names'))
  
  df1 <-dat1[,-1]
  dfnorder = as.data.frame(reshape2::melt(df1, id.vars = c("Group")))

  sumorder<-plyr::ddply(dfnorder,c("Group","variable"),summarise,
                        allmean=mean(value,na.rm=TRUE),sd1=sd(value,na.rm = TRUE),
                        n1=sum(!is.na(value)),se1=sd1/sqrt(n1))
  
  max_mean_group<- sumorder %>% group_by(variable)  %>% slice_max(allmean, n =1) %>%
   select(variable, enrichment = Group, max_mean = allmean)
  
  sumorder2<-
    full_join(sumorder, max_mean_group, by = 'variable',relationship = "many-to-many")
  write.csv(as.data.frame(sumorder2),paste0(label,".",analysisgroup,"-SUM.csv"))
  
  sumorder_wide<- sumorder %>%  select(variable, Group, allmean)  %>%
   pivot_wider(names_from = Group, values_from = c(allmean))
  write.csv(as.data.frame(sumorder_wide),paste0(label,".",analysisgroup,"-wider-means.csv"))
  
  sumorder_wide<- sumorder %>%  select(variable, Group, se1)  %>%
  pivot_wider(names_from = Group, values_from = c(se1))
 write.csv(as.data.frame(sumorder_wide),paste0(label,".",analysisgroup,"-wider-se.csv"))
  
  
  orderstat.test<- dfnorder %>% group_by(variable)%>% t_test(value ~ Group)%>%adjust_pvalue(method = "fdr") %>%
    add_significance(p.col = "p") # Benjamini & Yekutieli 
  write.csv(as.data.frame(orderstat.test),paste0(label,".",analysisgroup,"-df-t-test.csv"))
  
  return(sumorder)
}


NET_ANA <- function (otu,metadata,analysisgroup = "TimeGrowth" ,fileter = 0.6){
    
    library("phyloseq")
    library("NetCoMi")
  sampFile = as.data.frame(metadata[, c(analysisgroup)], row.names = row.names(metadata))
  colnames(sampFile)[1] = "Group"
  groupall<-unique(sampFile$Group)
  
  for (grou in groupall) {
    print(grou)
    sampFile1<-sampFile[sampFile$Group==grou,,drop = F]
    indxx<-rownames(sampFile1)
    otu1<-otu[,indxx, drop = F]
    filterrownum <- ncol(otu1)*fileter
    zero_counts <- as.data.frame(rowSums(otu1==0))
    colnames(zero_counts)[1] = "taxnum"
  #filtered_rows <- as.data.frame(rowSums(zero_counts>filterrownum ))
    filtered_rows <- subset(zero_counts, taxnum <= filterrownum )
    idx = rownames(otu1) %in% rownames(filtered_rows)
    otu2 = t(as.matrix(otu1[idx, , drop = F]))
    OTU =  otu_table(otu2, taxa_are_rows = TRUE)
    physeq = phyloseq(OTU)  
    net_genus <- netConstruct(physeq, measure="sparcc",zeroMethod = "pseudoZO",normMethod = "clr",sparsMethod="threshold",thresh=0.8)
    edge<-net_genus$edgelist1
    write.csv(as.data.frame(edge), file =paste0(grou,".","edge.csv") , row.names = FALSE)
  }
    return(groupall)
                     }