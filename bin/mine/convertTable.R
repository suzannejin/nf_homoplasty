# https://stackoverflow.com/questions/37756389/r-convert-multiple-columns-into-single-column

#install.packages("data.table")
library(data.table)

setwd("~/Downloads")
aligner="clustalo"


## import data
tc <- read.csv(file=paste(aligner,"_TC.csv",sep=""),header=TRUE,sep=",")
hom <- read.csv(file=paste(aligner,"_HOM.csv",sep=""),header=TRUE,sep=",")
w_hom <- read.csv(file=paste(aligner,"_W_HOM.csv",sep=""),header=TRUE,sep=",")
w_hom2 <- read.csv(file=paste(aligner,"_W_HOM2.csv",sep=""),header=TRUE,sep=",")
ngap <- read.csv(file=paste(aligner,"_NGAP.csv",sep=""),header=TRUE,sep=",")
ngap2 <- read.csv(file=paste(aligner,"_NGAP2.csv",sep=""),header=TRUE,sep=",")


