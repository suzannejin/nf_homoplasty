
## Libraries & packages
library(data.table)
library(ggplot2)
library(plotly)
#library(gridExtra)
library(ggpubr)


# Import data
fil=read.csv("merged_data.csv",header=TRUE,sep=" ",row.names=NULL)

metrics=c("tc","hom","w_hom","w_hom2","ngap","ngap2")



