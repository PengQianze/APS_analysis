### APS_analysis

#############################################
APS_analysis 中包括微生物组可视化的相关函数 （全部为R语言）
使用data文件夹中的数据可进行相关示例分析
microbiome-analysis-code 文件夹中包括从原始下机数据到丰度表的分析流程 （如果需要）

#############################################
APS_analysis includes related functions for microbiome visualization (all written in R language).
The data folder contains datasets that can be used for example analyses.
The microbiome-analysis-code folder contains the analysis workflow from raw sequencing data to abundance tables (if required).

##########依赖关系，完整使用需安装以下R包
########## Dependency Requirements: Full functionality requires installation of the following R packages
install.packages(reshape2)
install.packages(ggplot2)
install.packages(tibble)
install.packages( rstatix)
install.packages(vegan)
install.packages( ggpubr )
install.packages(ggrepel)
install.packages( plyr   )
install.packages(ggVennDiagram)
install.packages(scales)
install.packages(tidyr)
install.packages(phyloseq)
install.packages(NetCoMi)
install.packages(randomForest)

###########快速使用
#################Quick Start
1.拷贝sample_using_data中的数据,Visu_Function_v1.r,Sample_used.r 文件至相同文件夹
1.Copy the data from the sample_using_data folder, along with Visu_Function_v1.r and Sample_used.r files, into the same working directory

2.使用R语言运行Sample_used.r
2.Run Sample_used.r using R language

###########自定义使用
##################Customized Usage

1.完整environmental microbiome 的数据文件放于full_using文件夹
1.Full data files of environmental microbiome in the full_using folder

2.完整environmental microbiome 的脚本于full_using.r （可根据自己的需求来选择相应的脚本运行）
2.Use full_using.r for environmental microbiome analysis workflows (you may select specific scripts according to your requirements)

