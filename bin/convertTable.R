# https://stackoverflow.com/questions/37756389/r-convert-multiple-columns-into-single-column


## Libraries & packages
library(data.table)
library(ggplot2)
library(plotly)
#library(gridExtra)
library(ggpubr)


setwd("~/homoplasy/data")

metrics=c("tc","hom","w_hom","w_hom2","ngap","ngap2")
files=c(paste(aligner,"_TC.csv",sep=""),
        paste(aligner,"_HOM.csv",sep=""),
        paste(aligner,"_W_HOM.csv",sep=""),
        paste(aligner,"_W_HOM2.csv",sep=""),
        paste(aligner,"_NGAP.csv",sep=""),
        paste(aligner,"_NGAP2.csv",sep=""))

aligner="clustalo"



## Import data
import_data=function(aligner,metrics,files){
    "
    Input  -  aligner : name of the aligner (clustalo, )
              metrics : vector with the name of the metrics (tc, hom, w_hom, w_hom2, ngap, ngap2)
              files   : vector with the name of the .csv files to be imported
    Output -  dataList: list of imported data 
    "
    dataList=list()
    n=length(metrics)
    for (i in 1:n){
        metric=metrics[i]
        fil=files[i]
        dataList[[i]]=read.csv(file=fil,header=TRUE,sep=",")
    }
    names(dataList)=metrics
    return(dataList)
}



## Convert to 1 column dataframe
to1coldataList=function(dataList,metrics){
    "
    Input  - dataList   : list of imported data 
             metrics    : vector with the name of the metrics (tc, hom, w_hom, w_hom2, ngap, ngap2)
    Output - df1colList : list of 1 column dataframe
    "
    df1colList=list()
    n=length(metrics)
    for (i in 1:n){
        data=dataList[[i]]
        nam=metrics[i]
        setDT(data)
        df1colList[[i]]=melt(data, id = 1, value.name = nam, variable.name = "tree")
    }
    names(df1colList)=metrics
    return(df1colList)
}



## Combine data
## https://stackoverflow.com/questions/6709151/how-do-i-combine-two-data-frames-based-on-two-columns
combine_datasets=function(df1colList){
    n=length(df1colList)
    mergedDf=merge(df1colList[[1]], df1colList[[2]], by.x=c("Family", "tree"), by.y=c("Family", "tree"))
    for (i in 3:n){
        mergedDf=merge(mergedDf, df1colList[[i]], by.x=c("Family", "tree"), by.y=c("Family", "tree"))
    }
    return(mergedDf)
}



# Plot tc vs metrics
gplot=function(df,xs){
    plotList=list()
    n=length(xs)
    for (i in 1:n){
        plotList[[i]]=ggplot(df,aes_string(x=xs[i],y="tc",color="tree")) +
            geom_point() +
            geom_smooth(method = "lm", fill = NA)+
            geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
            ggtitle("ALL") +
            xlab(xs[i]) + ylab("") +
            ylim(0,100)
    }
    names(plotList)=xs
    return(plotList)
}


get_delta=function(df,metrics){
    trees=df$tree;fams=df$Family
    ntree=length(trees);nfam=length(fams)
    nmetrics=length(metrics)
    deltaArray=array(0,dim=c(ntree,ntree,nmetrics,nfam),dimnames=list(trees,trees,metrics,fams))

    for (f in 1:nfam)
        for (m in 1:nmetrics){
            for (i in 1:ntree){
                for (j in 1:ntree){
                    tree1=trees[i];tree2=trees[j];metric=metrics[m];fam=fams[f]
                    df2=df[which(df$Family==fam),]
                    value1=as.numeric(df2[which(df2$tree==tree1),])[m]
                    value2=as.numeric(df2[which(df2$tree==tree2),])[m]
                    delta=value1-value2
                    deltaArray[i,j,m,f]=delta
                }
            }
        }
    }
    return(deltaArray)
}



# Run
dataList=import_data(aligner,metrics,files)
df1colList=to1coldataList(dataList,metrics)
mergedDf=combine_datasets(df1colList)




df=mergedDf
plotList=gplot(df,metrics[2:6])


## Here you can explore the plots by:
#plotList$ngap
#ggplotly(plotList$ngap)




#######
### GRID WITH ALL PLOTS
#######

figure=ggarrange(plotList$ngap,plotList$ngap2,plotList$hom,plotList$w_hom,plotList$w_hom2,common.legend=TRUE)

annotate_figure(figure,
                top = text_grob(paste("Alignment with ",toupper(aligner),sep=""), size = 14),
                left = text_grob("TC Score", rot = 90)
)

