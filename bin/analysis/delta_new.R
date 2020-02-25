
## Libraries & packages
#install.packages('ggplot2')
library(ggplot2)
#install.packages('plotly')
library(plotly)
#install.packages('gridExtra')
library(gridExtra)
#install.packages('ggpubr')
library(ggpubr)


bucket=100
typ="raw"

ori=paste("data/delta/bucket",bucket,sep="")

tc=read.csv(file=paste(ori,"/tc.",typ,".csv",sep=""), header=TRUE, sep=",")
sp=read.csv(file=paste(ori,"/sp.",typ,".csv",sep=""), header=TRUE, sep=",")
homo=read.csv(file=paste(ori,"/homo.",typ,".csv",sep=""), header=TRUE, sep=",")
whomo=read.csv(file=paste(ori,"/whomo.",typ,".csv",sep=""), header=TRUE, sep=",")
whomo2=read.csv(file=paste(ori,"/whomo2.",typ,".csv",sep=""), header=TRUE, sep=",")
ngap=read.csv(file=paste(ori,"/ngap.",typ,".csv",sep=""), header=TRUE, sep=",")
ngap2=read.csv(file=paste(ori,"/ngap2.",typ,".csv",sep=""), header=TRUE, sep=",")



# Import data
filename=paste(metric)
fil=read.csv("delta_data.csv",header=TRUE,sep=" ",row.names=NULL)