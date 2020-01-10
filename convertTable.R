# https://stackoverflow.com/questions/37756389/r-convert-multiple-columns-into-single-column


install.packages("data.table")
library(data.table)

## import data
tc_homo <- read.csv(file="clustalo_ref_TC.csv", header=TRUE, sep=",")
hom_homo <- read.csv(file="clustalo_ref_HOM.csv", header=TRUE, sep=",")


## convert to 1 column dataframe
setDT(tc_homo)
fullTC=melt(tc_homo, id = 1, value.name = "tc", variable.name = "delta")

setDT(hom_homo)
fullHOM=melt(hom_homo, id = 1, value.name = "hom", variable.name = "delta")

#combine both datasdts
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
total = merge(fullTC, fullHOM, by.x=c("family", "delta"), by.y=c("family", "delta"))

## analisys
install.packages("ggplot2")
library(ggplot2)
install.packages("plotly")
library(plotly)

p_trees <- ggplot(data = total, aes(x = hom, y = tc, color = delta)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)+
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

ggplotly(p_trees)

#remove RND trees
mBed_famsa =subset(total, delta=="mBed.FAMSA")
mBed_pt =subset(total,  delta=="mBed.PartTree" )
pt_famsa =subset(total, delta=="PartTree.FAMSA")
homoplasy_nonRND = rbind(mBed_famsa, mBed_pt, pt_famsa)

## count how many points in each sector

## upp - left
length(which(homoplasy_nonRND$hom<0 & homoplasy_nonRND$tc>0))
## upp - right
length(which(homoplasy_nonRND$hom>0 & homoplasy_nonRND$tc>0))
## down - left
length(which(homoplasy_nonRND$hom<0 & homoplasy_nonRND$tc<0))
## down - right
length(which(homoplasy_nonRND$hom>0 & homoplasy_nonRND$tc<0))



