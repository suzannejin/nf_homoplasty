

args = commandArgs(trailingOnly=TRUE)

# Read input raw table
rawtable=read.table(args[1],header=TRUE)

metrics=c("tc","sp","col","homo","whomo","whomo2","ngap","ngap2")



sumByAlnTree = function(df,aligner,tree,bucket,metrics)
{
    pos=which(df$aligner==aligner & df$tree==tree & df$bucket==bucket)
    df2=df[pos,]

    meanDf=data.frame(aligner=aligner,tree=tree,bucket=bucket)
    for (i in metrics)
    {
        tmp=data.frame(mean(get(i,df2)))
        colnames(tmp)=i
        meanDf=cbind(meanDf,tmp)
    }
    return(meanDf)
}

n=3+length(metrics)  # aligner, tree, bucket, <metrics ...>
#meanDf=data.frame(matrix(ncol=n,nrow=0),dimnames=list(c(),c("aligner","tree","bucket",metrics)),stringsAsFactors=F)
#colnames(meanDf)=c("aligner","tree","bucket",metrics)
meanDf=read.table(text="",col.names=c("aligner","tree","bucket",metrics))
for (aligner in levels(df$aligner))
{
    for (tree in levels(df$tree))
    {
        for (bucket in levels(df$bucket))
        {
            df=rawtable
            tmp=sumByAlnTree(df,aligner,tree,bucket)
            meanDf=rbind(meanDf,tmp)
        }
    }
}

