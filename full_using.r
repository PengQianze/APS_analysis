
##############示例文件 
##############EXAMPLE feild


######加载函数

source("Visu_Function_v1.R")
################################################################################FIG1

###################稀释曲线
library("amplicon")

leaf1=read.csv("leaf_microbiome.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
root1=read.csv("root_microbiome.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

leaf = alpha_sample_rare(leaf1, length=10, rep=10, count_cutoff=1)
root = alpha_sample_rare(root1, length=10, rep=10, count_cutoff=1)
leaf
root
ggsave(paste0("leaf.pdf"), leaf, width=7, height=3, units="in")
ggsave(paste0("root.pdf"), root, width=7, height=3, units="in")

#####################################
library(reshape2)
library(ggplot2)
library(tibble)

metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("Order RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

colorscat <-c("#a4cfa9","#B2FA3E","#ABE6F4","#7FCFF5","#3DCBF7",
"#90FF8C","#C4C4C2","#FFD071","#01FF7C","#A2A2A2",
"#FFFB00","#FFD071","#7E7E7E","#BBBBBB",
"#dca759","#956632","#b38b4e","#b27e36","#b29569","#74c3aa","#FF8247","#cd903d",
"#686046","#41584c","#6c4a30","#83b072","#76cad9","#92a159","#89854d")
level =c("R","L")
shapes<-c(L=21, R =21)
analysisgroup =  "ORG"

#堆叠图 

STACK_ANA(taxnonmy,metadata,analysisgroup, label = "Organ_order_", level = level,
colorsF = colorscat,   topN = 10 )

#alpha （根/叶，根际，叶际PCOA） 

alphadiversity=read.csv("alpha.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

colors <-c(L="#a4cfa9",R ="#B2FA3E")
level =c("R","L")
shapes<-c(L=21, R =21)
analysisgroup = "ORG"

box_ANA (alphadiversity,metadata,analysisgroup, col_name= "shannon", colors ,level,shapes,label="Organ")
box_ANA (alphadiversity,metadata,analysisgroup, col_name = "chao1", colors ,level,shapes,label="Organ")
box_ANA (alphadiversity,metadata,analysisgroup, col_name = "simpson", colors ,level,shapes,label="Organ")
box_ANA (alphadiversity,metadata,analysisgroup, col_name = "richness", colors ,level,shapes,label="Organ")

#多样性差异 
taxonomy=read.csv("all_OTU .RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
colors <-c(L="#a4cfa9",R ="#B2FA3E")
level =c("L","R")
shapes<-c(L=21, R =21)
analysisgroup = "ORG"
PCOA_ANA  (taxonomy,metadata,analysisgroup ,colors,level,shapes,label="Organ_OTU_")

#多样性差异 
analysisgroup = "ORG"
VEEN_ANA (taxonomy,metadata,analysisgroup = analysisgroup)

#############共有饼图差异图
data1=read.csv("core.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
cut_off_FDR = 0.05
cut_off_log2FC = 3
COL <-  c("L_allmean","R_allmean","Pvalue","name")
color <- c("#B2FA3E", "#d2dae2", "#a4cfa9")

volcano_plot  (data1, cut_off_FDR, cut_off_log2FC, COL = COL,  color,label="ana_",topn=5)


data1=read.csv("ana_core_otu.csv", header=T, row.names= NULL , comment.char="", stringsAsFactors=F)
#analysis <- data1[data1$Group  == 'L', ]

topN  = 10

COL <- c("Order","L_allmean")
analysis <- data1[data1$Sig  == 'Up', ]
pie_plot  (analysis, COL, topN,colorsw = colorscat,label="L_UP_core")


COL <- c("Order","R_allmean")
analysis <- data1[data1$Sig  == 'Down', ]
pie_plot  (analysis, COL, topN,colorsw = colorscat,label="R_UP_core")

################################特异OTU

data1=read.csv("R_sp.csv", header=T, row.names= NULL , comment.char="", stringsAsFactors=F)
topN  = 10
COL <- c("Order","allmean")
pie_plot  (analysis=data1, COL, topN,colorsw = colorscat,label="R_sp")

data1=read.csv("L_sp.csv", header=T, row.names= NULL , comment.char="", stringsAsFactors=F)
pie_plot  (analysis=data1, COL, topN,colorsw = colorscat,label="L_sp")

################################################################################FIG2

#alpha （根/叶，根际，叶际PCOA） 

metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

#多样性差异 
alpha=read.csv("alpha.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata2 <- metadata[metadata$ORG == "L", ]
metadata3 <- metadata[metadata$ORG == "R", ]
taxonomy=read.csv("all_OTU .RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

#################cpcoa
analysisgroup = "Group"
colors <-c(Seedling="#a4cfa9",Tillering ="#a4cfa9",Elongation ="#a4cfa9")
level =c("Seedling","Tillering","Elongation")
shapes<-c(Seedling=22,Tillering =23,Elongation =24)
CPCOA_ANA  (taxonomy,metadata2,analysisgroup ,colors,shapes = shapes,label="L")

colors <-c(Seedling="#B2FA3E",Tillering ="#B2FA3E",Elongation ="#B2FA3E")
CPCOA_ANA  (taxonomy,metadata3,analysisgroup ,colors,shapes,label="R")

#####差异分析
#ALPHA
analysisgroup = "Group"
DF_ANA (taxonomy = t(alpha),metadata2,analysisgroup,repea=2,label="L_alpha")
DF_ANA (taxonomy = t(alpha),metadata3,analysisgroup,repea=2,label="R_alpha")

DF_ANA (taxonomy =taxonomy,metadata2,analysisgroup,repea=0,label="L_otu")
DF_ANA (taxonomy =taxonomy,metadata3,analysisgroup,repea=0,label="R_otu")


#RF特征棒状图

COL <- c("Order","MeanDecreaseGini","Group")
Topn  = 20
shapes<-c(Seedling=22,Tillering =23,Elongation =24)

data1=read.csv("L_RF_OTU.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
bb_plot (data1, COL, Topn,shapes=shapes, colors=colorscat,label="otu_L_RF") 

data1=read.csv("R_RF_OTU.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
bb_plot (data1, COL, Topn,shapes=shapes, colors=colorscat,label="otu_R_RF") 

###########################all_network
otu <- read.table("OTU.csv", header = TRUE, sep = ",", row.names = 1)
metadata <- read.table("metadata.csv", header = TRUE, sep = ",", row.names = 1)
analysisgroup = "TimeGrowth"

NET_ANA (otu,metadata,analysisgroup)

############################no_aps_network
otu <- read.table("no_ animal _OTU .csv", header = TRUE, sep = ",", row.names = 1)
metadata <- read.table("metadata.csv", header = TRUE, sep = ",", row.names = 1)
analysisgroup = "TimeGrowth"

NET_ANA (otu,metadata,analysisgroup)

###################################################网络特征气泡图

library(ggplot2)

net<-read.csv("net.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
colors1 <-c(L="#a4cfa9",R ="#B2FA3E")
shapes<-c(Seedling=22,Tillering =23,Elongation =24)
colors2 <- c(no="#000000",all ="#00AEFF")


p <-ggplot(net,aes(clustering_coefficient,node)) +
geom_point(aes(fill=Organ,colour =anmi,size = average_degree,shape = Develop),alpha=0.9)+
scale_fill_manual(values = colors1)+
scale_colour_manual(values = colors2)+
scale_shape_manual(values = shapes)+
scale_size_continuous(range=c(5,10)) +  
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
p

ggsave(paste0("net_struc",".pdf"), p, width=5, height=5.5, units="in")
ggsave(paste0("net_struc",".png"), p, width=5, height=5.5, units="in")

#################hub
#################hub

library(ggrepel)
library(ggplot2)
library( tibble)

net2<-read.csv("hub.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)


p2 <-ggplot(net2,aes(betweenesscentrality,closnesscentrality)) +
geom_point(aes(fill=organ,colour = anmimal,size = degree,shape = development),alpha=0.8)+
scale_fill_manual(values = colors1)+
scale_colour_manual(values = colors1)+
scale_shape_manual(values = shapes)+
scale_size_continuous(range=c(5,10)) +  
theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",linewidth = 0.25),
        axis.ticks = element_line(color = "black",linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle = 0),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       )+
geom_label_repel(data = net2, aes(label=label, fill=order), size=1.8, color="white",
segment.size=0.25,segment.color = "black", nudge_y=0.05, nudge_x=100,direction="x", hjust=0.8)
p2
p2
ggsave(paste0("hubs",".pdf"), p2, width=8, height=8, units="in")
ggsave(paste0("hubs",".png"), p2, width=8, height=8, units="in")

##########################################################FIG3
##########################################################FIG3

library(reshape2)
library(ggplot2)
library(tibble)

metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("FAPROTAX _.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

level =c("R","L")
shapes<-c(L=21, R =21)
analysisgroup =  "ORG"
#堆叠图 
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "RA_FUNC_", level = level,
colorsF = colorscat,   topN = 10 )

#功能多样性差异 （根/叶，根际，叶际PCOA）
metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata2 <- metadata[metadata$ORG == "L", ]
metadata3 <- metadata[metadata$ORG == "R", ]

taxonomy=read.csv("faprotax.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

colors <-c(L="#a4cfa9",R ="#B2FA3E")
level =c("L","R")
shapes<-c(L=21, R =21)
analysisgroup = "ORG"
PCOA_ANA  (taxonomy,metadata,analysisgroup ,colors,level,shapes,label="organ")

analysisgroup = "Group"
colors <-c(Seedling="#a4cfa9",Tillering ="#a4cfa9",Elongation ="#a4cfa9")
level =c("Seedling","Tillering","Elongation")
shapes<-c(Seedling=22,Tillering =23,Elongation =24)

PCOA_ANA  (taxonomy,metadata2,analysisgroup ,colors,level,shapes,label="L")
CPCOA_ANA  (taxonomy,metadata2,analysisgroup ,colors,level,shapes,label="L")

colors <-c(Seedling="#B2FA3E",Tillering ="#B2FA3E",Elongation ="#B2FA3E")
PCOA_ANA  (taxonomy,metadata3,analysisgroup ,colors,level,shapes,label="R")
CPCOA_ANA  (taxonomy,metadata3,analysisgroup ,colors,level,shapes,label="R")

#热图，柱状图
library(reshape2)
library(ggplot2)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(pheatmap)

amplic1<-read.csv("FUNCL.csv",header=T,row.name=1)
col_plot <-amplic1[,1:2]
col_plot <- cbind(Functionname = row.names(col_plot), col_plot)

class= c( C = "#00C3FF", N ="#E2F9AE" , S = "#83DDF8",
animal_associated = "#a4cfa9" , Other = "#B2FA3E")
p3 <- ggplot(col_plot ) +
    geom_col(aes(x = MeanDecreaseGini , y = reorder(Functionname, MeanDecreaseGini, decreasing = FALSE), fill =  class) )+
    scale_fill_manual(values = class) +
    labs(x = 'MeanDecreaseGini', y = '')+
    theme_classic()
p3  

ggsave(paste0("l_MeanDecGini",".pdf"), p3, width=6, height=4, units="in")
ggsave(paste0("l_MeanDecGini",".png"), p3, width=6, height=4, units="in")  

#热图 10:2 p

amplic<-as.matrix(amplic1[,3:5])
annotation_col = as.data.frame(amplic1[,2],row.names = row.names(amplic1))#创建分组列
colnames(annotation_col)[1] = "class"
ann_colors = list(
#Group = c(L = "#5e9162", R = "#f5c870"),
class= c( C = "#00C3FF", N ="#E2F9AE" , S = "#83DDF8",
animal_associated = "#a4cfa9" , Other = "#B2FA3E")) #定义分组颜色

p <- pheatmap(amplic, 
        color = colorRampPalette(c("white","#a4cfa9"))(100),
        scale = "row",
         cluster_cols = FALSE,  #取消列聚类，表示行聚类使用皮尔森相关系数聚类
        cluster_rows = FALSE,  #取消列聚类，表示行聚类使用皮尔森相关系数聚类  
         treeheight_row = 15, # 设置行聚类树高
         clustering_distance_rows = "manhattan", # "euclidean",
         clustering_method = "ward.D", # "single" ,"ward.D"
#         treeheight_col = 15,
       #   cutree_rows = FALSE, #根据样品列聚类情况将热图的行方向隔开为3份
         cellwidth = 30,cellheight = 20, # 设置热图方块宽度和高度
         display_numbers = F, # 热图上显示数值
       #   fontsize_number = 8, #热图上数值的字体大小
         main="Function ", # 设置图形标题
         show_colnames = T, # 设置行列标签的显示
         show_rownames = T,
         border= NA, # 设置边框为白色
         legend = T, # FALSE去除图例; T显示图例
         fontsize_row = 10, # 分别设置行列标签字体大小
         fontsize_col = 10,
#         angle_col = 45, # 设置标签显示角度
         annotation_row  = annotation_col, #显示样品列的分组信息及图例
         annotation_colors = ann_colors, #使用annotation_colors参数设定样品列分组的颜色
#       filename ="FUNCL.pdf" # 自动保存到设置路径下
         )

##############################################################

amplic1<-read.csv("FUNCR.csv",header=T,row.name=1)
col_plot <-amplic1[,1:2]
col_plot <- cbind(Functionname = row.names(col_plot), col_plot)
amplic<-as.matrix(amplic1[,3:5])


class= c( C = "#00C3FF", N ="#E2F9AE" , S = "#83DDF8",H = "#83b072",
animal_associated = "#a4cfa9" , Other = "#B2FA3E")


p3 <- ggplot(col_plot ) +
    geom_col(aes(x = MeanDecreaseGini , y = reorder(Functionname, MeanDecreaseGini, decreasing = FALSE), fill =  class) )+
    scale_fill_manual(values = class) +
    labs(x = 'MeanDecreaseGini', y = '')+
    theme_classic()
p3  

ggsave(paste0("R_MeanDecGini",".pdf"), p3, width=6, height=4, units="in")
ggsave(paste0("R_MeanDecGini",".png"), p3, width=6, height=4, units="in")  

#热图 10:2 p

annotation_col = as.data.frame(amplic1[,2],row.names = row.names(amplic1))#创建分组列
colnames(annotation_col)[1] = "class"
ann_colors = list(
#Group = c(L = "#5e9162", R = "#f5c870"),
class= c( C = "#00C3FF", N ="#E2F9AE" , S = "#83DDF8",H = "#83b072",
animal_associated = "#a4cfa9" , Other = "#B2FA3E")) #定义分组颜色

p <- pheatmap(amplic, 
        color = colorRampPalette(c("white","#a4cfa9"))(100),
        scale = "row",
         cluster_cols = FALSE,  #取消列聚类，表示行聚类使用皮尔森相关系数聚类
        cluster_rows = FALSE,  #取消列聚类，表示行聚类使用皮尔森相关系数聚类  
         treeheight_row = 15, # 设置行聚类树高
         clustering_distance_rows = "manhattan", # "euclidean",
         clustering_method = "ward.D", # "single" ,"ward.D"
#         treeheight_col = 15,
       #   cutree_rows = FALSE, #根据样品列聚类情况将热图的行方向隔开为3份
         cellwidth = 30,cellheight = 20, # 设置热图方块宽度和高度
         display_numbers = F, # 热图上显示数值
       #   fontsize_number = 8, #热图上数值的字体大小
         main="Function ", # 设置图形标题
         show_colnames = T, # 设置行列标签的显示
         show_rownames = T,
         border= NA, # 设置边框为白色
         legend = T, # FALSE去除图例; T显示图例
         fontsize_row = 10, # 分别设置行列标签字体大小
         fontsize_col = 10,
#         angle_col = 45, # 设置标签显示角度
         annotation_row  = annotation_col, #显示样品列的分组信息及图例
         annotation_colors = ann_colors, #使用annotation_colors参数设定样品列分组的颜色
#       filename ="FUNCL.pdf" # 自动保存到设置路径下
         )         


################################################################################FIG4


#功能pie_lopt——绝对量

#多样性差异 OTU 4:4

metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata2 <- metadata[metadata$ORG == "L", ]
metadata3 <- metadata[metadata$ORG == "R", ]

taxonomy=read.csv("animal _OTU .RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

colors <-c(L="#a4cfa9",R ="#B2FA3E")
level =c("L","R")
shapes<-c(L=21, R =21)
analysisgroup = "ORG"
PCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,shapes,label="organ")

analysisgroup = "Group"
colors <-c(Seedling="#a4cfa9",Tillering ="#a4cfa9",Elongation ="#a4cfa9")
level =c("Seedling","Tillering","Elongation")
shapes<-c(Seedling=22,Tillering =23,Elongation =24)
CPCOA_ANA  (taxonomy,metadata2,analysisgroup ,colors,shapes,label="L")

colors <-c(Seedling="#B2FA3E",Tillering ="#B2FA3E",Elongation ="#B2FA3E")
CPCOA_ANA  (taxonomy,metadata3,analysisgroup ,colors,shapes,label="R")

#物种分布图-

#堆叠图 TOP20 GENUS 8:8

level =c("L","R")
analysisgroup =  "ORG"
metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("animal . Genus .RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="Genus_animal")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "Genus_animal", level = level,colorsF = colorscat,  topN = 21 )

analysisgroup =  "TimeGrowth"
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="Genus_animal")

#物种分别热图

###数据整理 6:8

####作图
library(ggplot2)
dynamic= read.csv("Genus_ani.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

shapes<-c(Seedling=22,Tillering =23,Elongation =24)
colors <-c(L="#a4cfa9",R ="#B2FA3E")

p9<- ggplot(dynamic,aes(y=reorder(variable,top, decreasing = TRUE),x=reorder(STAGES,Stageord)))+
  geom_point(aes(size=allmean,
                 fill=enrichcompartment,
                 shape=enrichstages,
                 color=enrichcompartment
                 ))+
scale_color_manual(values = colors)+                
scale_fill_manual(values = colors)+
scale_shape_manual(values = shapes)+
scale_size_continuous(range=c(0,10)) +  
theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle = 45),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       )+
  labs(x=NULL,y=NULL)
p9
ggsave(paste0("heatplot_geus_ani",".pdf"), p9, width=6, height=8, units="in")
ggsave(paste0("heatplot_geus_ani",".png"), p9, width=6, height=8, units="in")


######order_level

level =c("L","R")
analysisgroup =  "ORG"
metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("Order _animal.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "Order_", level = level,colorsF = colorscat,  topN = 11 )
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="Order_animal")
analysisgroup =  "TimeGrowth"
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="Order_animal")
#物种分别热图

###数据整理 6:8
####作图
library(ggplot2)
dynamic= read.csv("Order_ani.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

shapes<-c(Seedling=22,Tillering =23,Elongation =24)
colors <-c(L="#a4cfa9",R ="#B2FA3E")

p9<- ggplot(dynamic,aes(y=reorder(variable,top, decreasing = TRUE),x=reorder(STAGES,Stageord)))+
  geom_point(aes(size=allmean,
                 fill=enrichcompartment,
                 shape=enrichstages,
                 color=enrichcompartment
                 ))+
scale_color_manual(values = colors)+                
scale_fill_manual(values = colors)+
scale_shape_manual(values = shapes)+
scale_size_continuous(range=c(0,10)) +  
theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",linewidth = 0.5),
        axis.ticks =  element_line(color = "black",linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle = 45),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       )+
  labs(x=NULL,y=NULL)
p9
ggsave(paste0("heatplot_order_ani",".pdf"), p9, width=6, height=8, units="in")
ggsave(paste0("heatplot_order_ani",".png"), p9, width=6, height=8, units="in")


######urban---filed---compried---aniaml---microbiome

library(dplyr)

taxonomy=read.csv("Genus-FIELD RA.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("Genus-URBAN RA.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "Genus")
total[is.na(total)]<-0
rownames(total)<-total$Genus
total$Genus<-NULL
write.csv(total, "total-genus.RA.csv")

metadata=read.csv("metadata-total.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("total-genus.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

colors <-c(Urban_soil="#a4cfa9",Field_soil ="#3DCBF7")
level =c("Urban_soil","Field_soil")
shapes<-c(Urban_soil=21, Field_soil =25)
analysisgroup = "Group"

PCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,label="genus")


DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="Genus_Urban_feild")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "Genus_Urban_feild", level = level,colorsF = colorscat,  topN = 21 )

#order_level

taxonomy=read.csv("Order-FILED RA.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("Order-URBAN RA.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "Order")
total[is.na(total)]<-0
rownames(total)<-total$Order
total$Order<-NULL
write.csv(total, "total-order.RA.csv")

metadata=read.csv("metadata-total.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("total-order.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

colors <-c(Urban_soil="#a4cfa9",Field_soil ="#3DCBF7")
level =c("Urban_soil","Field_soil")
shapes<-c(Urban_soil=21, Field_soil =25)
analysisgroup = "Group"
PCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,label="Order")

DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="Oder_Urban_feild")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "Order_Urban_feild", level = level,colorsF = colorscat,  topN = 10 )

####################################urban---filed---compried---function

taxonomy=read.csv("FAPROTAX _URBAN.RA.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("faprotax-FIELD .RA.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "func")
total[is.na(total)]<-0
rownames(total)<-total$func
total$func<-NULL
write.csv(total, "total-func.RA.csv")

metadata=read.csv("metadata-total.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("total-func.RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

colors <-c(Urban_soil="#a4cfa9",Field_soil ="#3DCBF7")
level =c("Urban_soil","Field_soil")
shapes<-c(Urban_soil=21, Field_soil =25)
analysisgroup = "Group"
PCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,label="func")

####6:6
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="func_Urban_feild")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "func_Urban_feild", level = level,colorsF = colorscat,  topN = 10 )

####################################urban---filed---compried---aniaml

library(dplyr)

taxonomy=read.csv("animal . Genus-filed .RA.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("animal . Genus-urban .RA.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "genus")
total[is.na(total)]<-0
write.csv(total, "total-genus.csv")

taxonomy=read.csv("animal . Order-filed .RA.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("animal . Order-urban .RA.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "order")
total[is.na(total)]<-0
write.csv(total, "total-order.csv")


##################差异分析
metadata=read.csv("metadata-total.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

colors <-c(Urban_soil="#a4cfa9",Field_soil ="#3DCBF7")
level =c("Urban_soil","Field_soil")
shapes<-c(Urban_soil=21, Field_soil =25)
analysisgroup = "Group"

##############order
taxonomy=read.csv("total-order.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="order_all_animal")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "order_all_animal", level = level,colorsF = colorscat,  topN = 10 )
PCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,label="urban_field_order")
##############genus
taxonomy=read.csv("total-genus.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="genus_all_animal")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "genus_all_animal", level = level,colorsF = colorscat,  topN = 20 )
PCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,label="urban_field_genus")


#####################
  library(rstatix)
  library(ggplot2)
  library( ggpubr )
metadata=read.csv("metadata-total.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

theme<-theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
               strip.text = element_text(size=12),
               axis.line = element_line(color = "black",linewidth = 0.25),
               axis.ticks =  element_line(color = "black",linewidth = 0.25),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(color="black",size=10),
               axis.text.x = element_text(color="black",size=10,angle = 45),
               panel.spacing.x = unit(0,"cm"),
               panel.border = element_blank(),
               panel.spacing = unit(0,"lines"))

colors <-c(Urban_soil="#a4cfa9",Field_soil ="#3DCBF7")
level =c("Urban_soil","Field_soil")
shapes<-c(Urban_soil=21, Field_soil =25)
analysisgroup = "Group"

alpha_dif = metadata
  df <- alpha_dif%>%mutate(Group=factor(Group,levels =level))
  dfn <- df
  names(dfn) <- c("Content", "Group")
  stat.test<- dfn %>% t_test(Content ~ Group)%>%
  add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")

 p1 <-ggplot(dfn,aes(Group,Content)) +
    geom_boxplot(aes(fill=Group,alpha = 0.8),alpha=1, width = 0.4,outlier.shape = NA, outlier.size = 0)+
    scale_fill_manual(values = colors)+
    scale_colour_manual(values = colors)+
    stat_pvalue_manual(stat.test,label = "p = {p}",hide.ns = F)+
    scale_size_continuous(range=c(1,3)) +  
    theme
  p1
  
  ggsave(paste0("animal_parasites_or_symbionts",analysisgroup,".pdf"), p1, width=3.6, height=3.4, units="in")
  ggsave(paste0("animal_parasites_or_symbionts",analysisgroup,".png"), p1, width=3.6, height=3.4, units="in")


####作图
library(ggplot2)
dynamic= read.csv("order_urban_field.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

colors <-c(Urban_soil="#a4cfa9",Field_soil ="#3DCBF7")
shapes<-c(Urban_soil=21, Field_soil =25)

p9<- ggplot(dynamic,aes(x=reorder(variable,top, decreasing = FALSE),y=reorder(Group,Stageord)))+
  geom_point(aes(size=allmean,
                 fill=enrichment,
                 shape=enrichment,
                 color=enrichment
                 ))+
scale_color_manual(values = colors)+                
scale_fill_manual(values = colors)+
scale_shape_manual(values = shapes)+
scale_size_continuous(range=c(0,10)) +  
theme+
labs(x=NULL,y=NULL)
p9
ggsave(paste0("heatplot_order_ani",".pdf"), p9, width=6, height=3, units="in")
ggsave(paste0("heatplot_order_ani",".png"), p9, width=6, height=3, units="in")

###################genus
dynamic= read.csv("genus_urban_field.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

colors <-c(Urban_soil="#a4cfa9",Field_soil ="#3DCBF7")
shapes<-c(Urban_soil=21, Field_soil =25)

p9<- ggplot(dynamic,aes(x=reorder(variable,top, decreasing = FALSE),y=reorder(Group,Stageord)))+
  geom_point(aes(size=allmean,
                 fill=enrichment,
                 shape=enrichment,
                 color=enrichment
                 ))+
scale_color_manual(values = colors)+                
scale_fill_manual(values = colors)+
scale_shape_manual(values = shapes)+
scale_size_continuous(range=c(0,10)) +  
theme+
labs(x=NULL,y=NULL)
p9
ggsave(paste0("heatplot_genus_ani",".pdf"), p9, width=8, height=3, units="in")
ggsave(paste0("heatplot_genus_ani",".png"), p9, width=8, height=3, units="in")


##################################################multi_area_microbiome

####################################urban---filed---compried---all

################################合并数据
library(dplyr)

taxonomy=read.csv("total-genus.RA.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("Genus RA.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "genus")
total[is.na(total)]<-0
write.csv(total, "ALL-genus.csv")

taxonomy=read.csv("total-order.RA.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("Order RA.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "order")
total[is.na(total)]<-0
write.csv(total, "ALL-order.csv")

taxonomy=read.csv("total-func.RA.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("faprotax .RA.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "func")
total[is.na(total)]<-0
write.csv(total, "ALL-function.csv")

#############################################

#####################差异分析###all--aps

  library(rstatix)
  library(ggplot2)
  library( ggpubr )
  metadata=read.csv("metadata2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

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
               panel.spacing = unit(0,"lines"))

colors <-c(Rice_L_Field="#3DCBF7",
ZH11_L_Field="#3DCBF7",
ZH11_R_Field="#3DCBF7",
ZH11_R_Urban="#a4cfa9",
ZH11_L_Urban="#a4cfa9",
Ningbo_L_Urban="#a4cfa9")
level <-c("Ningbo_L_Urban","ZH11_L_Urban","ZH11_R_Urban","ZH11_R_Field","ZH11_L_Field","Rice_L_Field")
shapes <-c(Rice_L_Field=22,
ZH11_L_Field=23,
ZH11_R_Field=25,
ZH11_R_Urban=1,
ZH11_L_Urban=21,
Ningbo_L_Urban=24)
analysisgroup = "Group2"

alpha_dif = metadata
  df <- alpha_dif%>%mutate(Group2=factor(Group2,levels =level))
  dfn <- df 
  names(dfn) <- c("Content", "Group")
  stat.test<- dfn %>% t_test(Content ~ Group)%>%
  add_significance(p.col = "p")%>% add_xy_position(fun = "max", x = "Group")

 p1 <-ggplot(dfn,aes(Group,Content)) +
    geom_boxplot(aes(fill=Group,alpha = 0.8),alpha=1, width = 0.4,outlier.shape = NA, outlier.size = 0)+
    scale_fill_manual(values = colors)+
    scale_colour_manual(values = colors)+
    stat_pvalue_manual(stat.test,label = " {p.signif}",hide.ns = F)+
    scale_size_continuous(range=c(1,3)) +  
    theme 
  p1
  
  ggsave(paste0("aps_all",analysisgroup,".pdf"), p1, width=6, height=6, units="in")
  ggsave(paste0("aps_all",analysisgroup,".png"), p1, width=6, height=6, units="in")



###############################

#order_level

metadata=read.csv("metadata2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("ALL-order.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

colors <-c(Rice_L_Field="#3DCBF7",
ZH11_L_Field="#3DCBF7",
ZH11_R_Field="#3DCBF7",
ZH11_R_Urban="#a4cfa9",
ZH11_L_Urban="#a4cfa9",
Ningbo_L_Urban="#a4cfa9")
level <-c("Ningbo_L_Urban","ZH11_L_Urban","ZH11_R_Urban","ZH11_R_Field","ZH11_L_Field","Rice_L_Field")
shapes <-c(Rice_L_Field=22,
ZH11_L_Field=23,
ZH11_R_Field=25,
ZH11_R_Urban=1,
ZH11_L_Urban=21,
Ningbo_L_Urban=24)
analysisgroup = "Group2"

CPCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,shapes = shapes,label="Order")

DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="Oder_all_Urban_feild")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "Order_all_Urban_feild", level = level,colorsF = colorscat,  topN = 10 )

####################################urban---filed---compried---function

metadata=read.csv("metadata2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("ALL-function.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0

colors <-c(Rice_L_Field="#3DCBF7",
ZH11_L_Field="#3DCBF7",
ZH11_R_Field="#3DCBF7",
ZH11_R_Urban="#a4cfa9",
ZH11_L_Urban="#a4cfa9",
Ningbo_L_Urban="#a4cfa9")
level <-c("Ningbo_L_Urban","ZH11_L_Urban","ZH11_R_Urban","ZH11_R_Field","ZH11_L_Field","Rice_L_Field")
shapes <-c(Rice_L_Field=22,
ZH11_L_Field=23,
ZH11_R_Field=25,
ZH11_R_Urban=1,
ZH11_L_Urban=21,
Ningbo_L_Urban=24)
analysisgroup = "Group2"

CPCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,shapes = shapes,label="func")

####6:6
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="func_all_Urban_feild")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "func_all_Urban_feild", level = level,colorsF = colorscat,  topN = 10 )

AFD14A nitrate reduction
animal parasites or symbionts ABDCE9
nitrogen respiration 7EC7EA
nitrogen fixation 47BFE6
phototrophy 9DCE82
ureolysis C4C4C1


####################################urban---filed---compried---aniaml

library(dplyr)

taxonomy=read.csv("total-genus.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("animal . Genus .RA2.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "genus")
total[is.na(total)]<-0
write.csv(total, "all-aps-genus.csv")

taxonomy=read.csv("total-order.csv", header=T,  comment.char="", stringsAsFactors=F)
taxonomy2=read.csv("animal . Order .RA2.csv", header=T,  comment.char="", stringsAsFactors=F)
total <- full_join(taxonomy, taxonomy2, by = "order")
total[is.na(total)]<-0
write.csv(total, "all-aps-order.csv")


##################差异分析
metadata=read.csv("metadata2.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

colors <-c(Rice_L_Field="#3DCBF7",
ZH11_L_Field="#3DCBF7",
ZH11_R_Field="#3DCBF7",
ZH11_R_Urban="#a4cfa9",
ZH11_L_Urban="#a4cfa9",
Ningbo_L_Urban="#a4cfa9")
level <-c("Ningbo_L_Urban","ZH11_L_Urban","ZH11_R_Urban","ZH11_R_Field","ZH11_L_Field","Rice_L_Field")
shapes <-c(Rice_L_Field=22,
ZH11_L_Field=23,
ZH11_R_Field=25,
ZH11_R_Urban=1,
ZH11_L_Urban=21,
Ningbo_L_Urban=24)
analysisgroup = "Group2"

##############order
taxonomy=read.csv("all-aps-order.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
DF_ANA (taxonomy = taxonomy,metadata,analysisgroup,repea=1,label="order_all_animal")
STACK_ANA(taxnonmy,metadata,analysisgroup, label = "order_all_animal", level = level,colorsF = colorscat,  topN = 10 )
CPCOA_ANA  (taxonomy,metadata,analysisgroup ,colors ,shapes = shapes,label="urban_field_order")


################################################################################FIG5
################################################################################FIG5
################################################################################FIG5
################################################################################FIG5
################################################################################FIG5

source("~/FUNCTION/visulization.R")

metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("no_ animal _OTU .RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy[is.na(taxonomy)]<-0
metadata2 <- metadata[metadata$ORG == "L", ]
metadata3 <- metadata[metadata$ORG == "R", ]
analysisgroup = "Group"

topn = 20
new_rownames <- gsub("-", "_", rownames(taxonomy))
new_rownames <- gsub("\\[|\\]", "", rownames(taxonomy))
rownames(taxonomy) <- new_rownames

analysisgroup = "Group"
RF_ANA  (taxonomy = taxonomy ,metadata2,analysisgroup = analysisgroup,topn,label="L")
RF_ANA  (taxonomy = taxonomy ,metadata3,analysisgroup = analysisgroup,topn,label="R")


#DF_ANA (taxonomy = taxonomy,metadata2,analysisgroup,repea=1,label="L_otu_all")
#DF_ANA (taxonomy = taxonomy,metadata3,analysisgroup,repea=1,label="R_otu_all")

#RF特征棒状图

COL <-  c("Order","MeanDecreaseGini","Group")
Topn  = 20
shapes<-c(Seedling=22,Tillering =23,Elongation =24)

data1=read.csv("noanimal_L_OTU_RF.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
bb_plot (data1, COL, Topn,shapes=shapes, colors=colorscat,label="NOAMotu_L_RF") 

data1=read.csv("noanimal_R_OTU_RF.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
bb_plot (data1, COL, Topn,shapes=shapes, colors=colorscat,label="NOAMotu_R_RF") 


########################微生物组变异贡献度分析
########################################################L_RDA
########################################################L_RDA
########################################################L_RDA

library(vegan)
library(ggplot2)
library( tibble)
library( rdacca.hp)

taxonomy=read.csv("all_OTU .RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata1=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata1

metadata = metadata1[metadata1$ORG == "L", ]
idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]

taxonomy[is.na(taxonomy)]<-0
taxonomy = taxonomy[, rownames(metadata)]
taxonomy = taxonomy[rowSums(taxonomy)!=0, ]

#RDA 
colors <-c(Seedling="#a4cfa9",Tillering ="#a4cfa9",Elongation ="#a4cfa9",
animal_parasites_or_symbionts="#A8D6E3",
intracellular_parasites = "#9DC4A2",
predatory_or_exoparasitic = "#F0C673")
shapes<-c(Seedling=22,Tillering =23,Elongation =24)

###############################################发育阶段影响

corrmetabo <-"Group"
sampFile = as.data.frame(metadata[, corrmetabo], row.names = row.names(metadata))
colnames(sampFile)[1] = corrmetabo
indiv_rda_all <- rda (t(taxonomy) ~ Group, data = sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)
eig2 = indiv_rda.scaling1$cont$importance
points2 = as.data.frame(indiv_rda.scaling1$sites)
points2 = cbind(metadata, points2)

p1 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = RDA2, fill= Group ,shape=Group),alpha = 0.8,size =6, stroke = 0.25,  color = "#000000") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("RDA 2 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  scale_fill_manual(values =colors )+
  scale_shape_manual(values = shapes )+
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
p1

ggsave(paste0("L_ALL_RDA_dev",".pdf"), p1, width=5.5, height=6, units="in")
ggsave(paste0("L_ALL_RDA_dev",".png"), p1, width=5.5, height=6, units="in")

########################################animal_parasites_or_symbionts
corrmetabo<-c(
"animal_parasites_or_symbionts"
)

sampFile = as.data.frame(metadata[, corrmetabo], row.names = row.names(metadata))
colnames(sampFile)[1] = "animal_parasites_or_symbionts"

indiv_rda_all <- rda(t(taxonomy)~animal_parasites_or_symbionts, sampFile)

#######################################微生物组变异
variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)
eig2 = indiv_rda.scaling1$cont$importance
points2 = as.data.frame(indiv_rda.scaling1$sites)
points2 = cbind(metadata, points2)

p2 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill= Group ,shape=Group),alpha = 0.8, size =6,  stroke = 0.25,  color = "#000000") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 1 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  scale_fill_manual(values =colors )+
  scale_shape_manual(values = shapes )+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",linewidth = 0.25),
        axis.ticks =  element_line(color = "black",linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle =0),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       )
p2

ggsave(paste0("L_ALL_RDA_Aps",".pdf"), p2, width=5.5, height=6, units="in")
ggsave(paste0("L_ALL_RDA_Aps",".png"), p2, width=5.5, height=6, units="in")

###################################################################微生物相关
library(psych)
corrr <- cbind(sampFile, t(taxonomy))

animal <- sampFile$animal_parasites_or_symbionts
corre <- sapply( corrr,function(x) cor.test (animal, x, method="spearman"))
write.csv(corre, file="L_ALL_func_micro_corr.csv")


###########################################################根际
###########################################################根际
###########################################################根际

taxonomy=read.csv("all_OTU .RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata1=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
metadata1

metadata = metadata1[metadata1$ORG == "R", ]
idx = rownames(metadata) %in% colnames(taxonomy)
metadata = metadata[idx, , drop = F]
taxonomy[is.na(taxonomy)]<-0
taxonomy = taxonomy[, rownames(metadata)]
taxonomy = taxonomy[rowSums(taxonomy)!=0, ]

#RDA CCA判断
colors <-c(Seedling="#B2FA3E",Tillering ="#B2FA3E",Elongation ="#B2FA3E",animal_parasites_or_symbionts="#A8D6E3",
intracellular_parasites = "#9DC4A2",
predatory_or_exoparasitic = "#F0C673")
shapes<-c(Seedling=22,Tillering =23,Elongation =24)


###############################################发育阶段影响
corrmetabo <-"Group"
sampFile = as.data.frame(metadata[, corrmetabo], row.names = row.names(metadata))
colnames(sampFile)[1] = corrmetabo
indiv_rda_all <- rda (t(taxonomy) ~ Group, data = sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)
eig2 = indiv_rda.scaling1$cont$importance
points2 = as.data.frame(indiv_rda.scaling1$sites)
points2 = cbind(metadata, points2)

p5 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill= Group ,shape=Group),alpha = 0.8, size =6, stroke = 0.25, color = "#000000") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 2 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  scale_fill_manual(values =colors )+
  scale_shape_manual(values = shapes )+
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
p5

ggsave(paste0("R_ALL_RDA_dev",".pdf"), p5, width=5.5, height=6, units="in")
ggsave(paste0("R_ALL_RDA_dev",".png"), p5, width=5.5, height=6, units="in")
########################################animal_parasites_or_symbionts
corrmetabo<-c(
"animal_parasites_or_symbionts"
)

sampFile = as.data.frame(metadata[, corrmetabo], row.names = row.names(metadata))
colnames(sampFile)[1] = "animal_parasites_or_symbionts"
indiv_rda_all <- rda(t(taxonomy)~animal_parasites_or_symbionts, sampFile)

variance_indiv = (indiv_rda_all$CCA$tot.chi/indiv_rda_all$tot.chi)
perm_indiv_rda = anova.cca(indiv_rda_all, permutations = 1000, parallel = 4)
p.val2 = perm_indiv_rda[1, 4]

indiv_rda.scaling1 <- summary(indiv_rda_all, scaling = 1)
eig2 = indiv_rda.scaling1$cont$importance
points2 = as.data.frame(indiv_rda.scaling1$sites)
points2 = cbind(metadata, points2)

p6 = ggplot(points2) + 
  geom_point(aes(x =  RDA1 , y = PC1, fill= Group ,shape=Group),alpha = 0.8,size =6,  stroke = 0.25, color = "#000000") + 
  labs(x = paste("RDA 1 (", format(100 * eig2[2,1], digits = 4), "%)", sep = ""), 
       y = paste("PC 1 (", format(100 * eig2[2,2], digits = 4), "%)", sep = ""), color = "black") + 
  ggtitle(paste(format(100 * variance_indiv, digits = 3), " % of variance; P = ", format(p.val2 , digits = 2), sep = "")) +
  scale_fill_manual(values =colors )+
  scale_shape_manual(values = shapes )+
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

p6

ggsave(paste0("R_ALL_RDA_Aps",".pdf"), p6, width=5.5, height=6, units="in")
ggsave(paste0("R_ALL_RDA_Aps",".png"), p6, width=5.5, height=6, units="in")


##################################################phyllosphere micorbiome Associated APS 
library(vegan)
library(ggplot2)
library(tibble)

taxonomy <- t(read.table("all_OTU_L .RA.csv", header = TRUE, sep = ",", row.names = 1))
sampFile <- t(read.table("animal _OTU_L .RA.csv", header = TRUE, sep = ",", row.names = 1))

otu_rda_all <- rda(taxonomy, sampFile)

otu_rda_all_envfit <- envfit(otu_rda_all,sampFile,permu=999)
cor_metabo <- cbind(as.data.frame(otu_rda_all_envfit$vectors$r),as.data.frame(otu_rda_all_envfit$vectors$pvals))
colnames(cor_metabo) = c("r", "p")
cor_metabo <-rownames_to_column(cor_metabo, var = "metabo")
write.csv(as.data.frame(cor_metabo), file =paste0("aps-response",".","l.csv") , row.names = FALSE)

################################################rhizoshpere micorbiome Associated APS

taxonomy <- t(read.table("all_OTU_R .RA.csv", header = TRUE, sep = ",", row.names = 1))
sampFile <- t(read.table("animal _OTU_R .RA.csv", header = TRUE, sep = ",", row.names = 1))

otu_rda_all <- rda(taxonomy, sampFile)

otu_rda_all_envfit <- envfit(otu_rda_all,sampFile,permu=999)
cor_metabo <- cbind(as.data.frame(otu_rda_all_envfit$vectors$r),as.data.frame(otu_rda_all_envfit$vectors$pvals))
colnames(cor_metabo) = c("r", "p")
cor_metabo <-rownames_to_column(cor_metabo, var = "metabo")
write.csv(as.data.frame(cor_metabo), file =paste0("aps-response",".","r.csv") , row.names = FALSE)

####################################

########################微生物相关
library(psych)
corrr <- cbind(sampFile, t(taxonomy))
animal <- sampFile$animal_parasites_or_symbionts
corre <- sapply( corrr,function(x) cor.test (animal, x, method="spearman"))

write.csv(corre, file="R_ALL_func_micro_COR.csv")


########################微生物相关
library(psych)
corrr <- cbind(sampFile, t(taxonomy))
animal <- sampFile$animal_parasites_or_symbionts
corre <- sapply( corrr,function(x) cor.test (animal, x, method="spearman"))

write.csv(corre, file="indi_func_micro_COR.csv")

########################

colorscat <-c(Bacteroidales="#a4cfa9",Clostridiales="#B2FA3E",Burkholderiales="#ABE6F4",Rhizobiales="#7FCFF5",Acidobacteriales="#3DCBF7",
Xanthomonadales="#90FF8C",Rhodospirillales="#C4C4C2",Sphingobacteriales="#FFD071",Pseudomonadales="#01FF7C",Other="#A2A2A2",
Methylophilales="#FFFB00",Caulobacterales="#BB5120",Oceanospirillales="#7E7E7E",Gemmatimonadales="#BBBBBB",Streptomycetales="#1F41EE",
Desulfuromonadales="#C96E6E",Fusobacteriales="#FCF4E7",Gaiellales="#3E51CA",Lactobacillales="#76cad9",AT425_EubC11_terrestrial_group="#1D7AEA",
Rhodobacterales="#dca759",TRA3_20="#956632",Rickettsiales="#b38b4e",Corynebacteriales="#b27e36",Flavobacteriales="#B67AE7",
Erysipelotrichales="#b29569",Myxococcales="#74c3aa",Ignavibacteriales="#FF8247",Campylobacterales="#cd903d",Bacillales = "#EE3C3C",
Unassigned="#686046",Selenomonadales="#41584c",Micrococcales="#6c4a30",Enterobacteriales="#83b072",Subgroup_3 = "#82705C",
Rhodocyclales="#672A8C",Spirochaetales="#92a159",Sphingomonadales="#89854d",Pseudonocardiales="#F781B0",
animal_parasites_or_symbionts="#A8D6E3",Solirubrobacterales = "#440769",
intracellular_parasites="#9DC4A2",
predatory_or_exoparasitic="#F0C673"
)

data1=read.csv("L_otu_aps_corr.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
p14 <- ggplot(
data1,aes(color=Order,fill=Order, x = r,y =reorder(OTU,r)))+
    geom_col( width =0.8, linewidth=0.25)+
    geom_vline(xintercept = 0,linewidth=0.5)+
    labs(x = '', y = '')+
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",linewidth = 0.5),
        axis.ticks = element_line(color = "black",linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle = 0),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       )+
    scale_color_manual(values=colorscat)+
    scale_fill_manual(values=colorscat)
p14

ggsave(paste0("L-aps_otu_cor.pdf"), p14, width=8, height=8, units="in")
ggsave(paste0("L-aps_otu_cor.png"),p14, width=8, height=8, units="in")

data1=read.csv("R_otu_aps_corr.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
p15 <- ggplot(
data1,aes(color=Order,fill=Order, x = r,y =reorder(OTU,r)))+
    geom_col( width =0.8, linewidth=0.25)+
    geom_vline(xintercept = 0, linewidth=0.5)+
    labs(x = '', y = '')+
    theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
        strip.text = element_text(size=12),
        axis.line = element_line(color = "black",linewidth = 0.5),
        axis.ticks = element_line(color = "black",linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle = 0),
        panel.spacing.x = unit(0,"cm"),
        panel.border = element_blank(),
       panel.spacing = unit(0,"lines")
       )+
    scale_color_manual(values=colorscat)+
    scale_fill_manual(values=colorscat)
p15

ggsave(paste0("R-aps_otu_cor.pdf"), p15, width=8, height=8, units="in")
ggsave(paste0("R-aps_otu_cor.png"),p15, width=8, height=8, units="in")




