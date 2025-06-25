
##############示例文件 
##############EXAMPLE sample used

source("Visu_Function_v1.R")
####data_analysis

#####################################
library(reshape2)
library(ggplot2)
library(tibble)

metadata=read.csv("metadata.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
taxonomy=read.csv("all_OTU .RA.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)

colorscat <-c("#a4cfa9","#B2FA3E","#ABE6F4","#7FCFF5","#3DCBF7",
"#90FF8C","#C4C4C2","#FFD071","#01FF7C","#A2A2A2",
"#FFFB00","#FFD071","#7E7E7E","#BBBBBB",
"#dca759","#956632","#b38b4e","#b27e36","#b29569","#74c3aa","#FF8247","#cd903d",
"#686046","#41584c","#6c4a30","#83b072","#76cad9","#92a159","#89854d")
level =c("R","L")
shapes<-c(L=21, R =21)
analysisgroup =  "ORG"
colors <-c(L="#d2dae2",R ="#B2FA3E")

box_ANA (alphadiversity=t(taxonomy),metadata,analysisgroup, col_name= "OTU_4", colors ,level,shapes,label="box_sample_")
PCOA_ANA (taxonomy,metadata,analysisgroup ,colors,shapes,label="PCOA_sample_") 
VEEN_ANA (taxonomy,metadata,analysisgroup = analysisgroup)
metadata2 <- metadata[metadata$ORG == "L", ]
analysisgroup = "Group"
colors <- c(Seedling="#a4cfa9",Tillering ="#a4cfa9",Elongation ="#a4cfa9")
level <- c("Seedling","Tillering","Elongation")
shapes <- c(Seedling=22,Tillering =23,Elongation =24)
STACK_ANA(taxonomy,metadata2,analysisgroup, label = "stack_sample_", level = level,colorsF = colorscat, topN = 20)
CPCOA_ANA  (taxonomy,metadata2,analysisgroup ,colors,shapes = shapes,label="CPCOA_sample")
DF_ANA (taxonomy = taxonomy,metadata2,analysisgroup,repea=0,label="DF_t_test")

###############################################

data1=read.csv("volcano_TEST.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
cut_off_FDR = 0.05
cut_off_log2FC = 3
COL <-  c("L_allmean","R_allmean","Pvalue","name")
color <- c("#B2FA3E", "#d2dae2", "#a4cfa9")
volcano_plot  (data1, cut_off_FDR, cut_off_log2FC, COL = COL,  color,label="volcano_sample_",topn=5)

data1=read.csv("pie_test.csv", header=T, row.names= NULL , comment.char="", stringsAsFactors=F)
topN  = 10  
COL <- c("Order","L_allmean")
analysis <- data1[data1$Sig  == 'Up', ]
pie_plot  (analysis, COL, topN,colorsw = colorscat,label="pie_test")

#########################
#RF特征棒状图
COL <- c("Order","MeanDecreaseGini","Group")
Topn <- 20
shapes <-c(Seedling=22,Tillering =23,Elongation =24)

data1=read.csv("bbplot_TEST.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
bb_plot (data1, COL, Topn,shapes=shapes, colors=colorscat,label="bbplot_sample") 

###########################all_network (need_long_time)
#analysisgroup = "TimeGrowth"
#NET_ANA (taxonomy,metadata,analysisgroup)

###########################bubble_plot
data1<-read.csv("bubble_TEST.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
colors1 <-c(L="#a4cfa9",R ="#B2FA3E")
shapes<-c(Seedling=22,Tillering =23,Elongation =24)
colors2 <- c(no="#000000",all ="#00AEFF")
COL <- c("clustering_coefficient","average_degree", "node" ,"Organ","anmi","Develop")
bubble_plot (data1,COL,colors1,shapes,colors2)

###########################shape_heatmap
dynamic= read.csv("shape_heatmap_TEST.csv", header=T, row.names=1, comment.char="", stringsAsFactors=F)
shapes<-c(Seedling=22,Tillering =23,Elongation =24)
colors <-c(L="#a4cfa9",R ="#B2FA3E")
COL <- c("variable", "top" ,"STAGES","Stageord","allmean","enrichcompartment","enrichstages")
shape_heatmap (data1=dynamic,COL,shapes,colors)

