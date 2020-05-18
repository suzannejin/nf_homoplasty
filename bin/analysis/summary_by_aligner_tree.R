
# -----
# Input
# -----

args = commandArgs(trailingOnly=TRUE)   # Rscript <script> <in csv> <outdir>

# Read input raw table
rawtable=read.table(args[1],header=TRUE,sep=",")

metrics=c("tc","sp","col","homo","whomo","whomo2","ngap","ngap2")


# ------------------
# Averaging function
# ------------------

# Function to obtain the average values 
sumByAlnTree = function(df,groupfam,combs,metrics)
{
    # Parse data by family type
    if (groupfam == "1000-10K")
    {
        pos = which(df$nseq >= 1000 & df$nseq <= 10000)
        df = df[pos,]
        rownames(df) = 1:nrow(df)
    }else if (groupfam == "top20")
    {
        fams = unique(df$family)
        top20 = fams[ (length(fams)-20) : length(fams) ]
        pos = which(df$family %in% top20)
        df = df[pos,]
        rownames(df) = 1:nrow(df)
    }

    # Parse data by combination
    if (length(combs) > 0)
    {
        pos = list()
        for (key in names(combs))
        {
            value = as.character(combs[key])
            pos[key] = list(which(df[key] == value))
        }
        if (length(combs) > 2)
        {
            for (i in 1:(length(combs)-1))
            {
                name1 = names(combs[i])
                name2 = names(combs[i+1])
                if (i == 1)
                {
                    pos2 = intersect(unlist(pos[name1]),unlist(pos[name2]))
                }else{
                    pos2 = intersect(pos2,unlist(pos[name2]))
                }
                
            }
        }else if (length(combs) == 2){
            pos2 = intersect(pos[[1]],pos[[2]])
        }else if (length(combs) == 1){
            pos2 = pos[[1]]
        }
        df2=df[pos2,]
        rownames(df2) = 1:nrow(df2)
    } else {
        df2 = df
    }

    # Calculate average values
    tmp_values = c(groupfam, as.character(combs))
    tmp_keys = c("groupfam", names(combs))
    meanDf = data.frame(matrix(tmp_values,nrow=1))
    names(meanDf) = tmp_keys
    for (i in metrics)
    {
        tmp=data.frame(mean(get(i,df2), na.rm=TRUE))
        colnames(tmp)=i
        meanDf=cbind(meanDf,tmp)
    }

    return(meanDf)
}


# -------------
# Create tables
# -------------

df = rawtable
meanDf_bucket_aligner_tree=read.table(text="",col.names=c("groupfam","bucket","aligner","tree",metrics))
meanDf_bucket_aligner=read.table(text="",col.names=c("groupfam","bucket","aligner",metrics))
meanDf_bucket=read.table(text="",col.names=c("groupfam","bucket",metrics))
meanDf=read.table(text="",col.names=c("groupfam",metrics))
for (groupfam in c("all","1000-10K","top20"))
{
    for (bucket in unique(df$bucket))
    {
        for (aligner in unique(df$aligner))
        {
            for (tree in unique(df$tree))
            {
                # With bucket:aligner:tree
                df = rawtable
                combs = c("bucket"=bucket, "aligner"=aligner, "tree"=tree)
                tmp=sumByAlnTree(df,groupfam,combs,metrics)
                meanDf_bucket_aligner_tree=rbind(meanDf_bucket_aligner_tree,tmp)
            }
            # With bucket:aligner
            df = rawtable
            combs = c("bucket"=bucket, "aligner"=aligner)
            tmp=sumByAlnTree(df,groupfam,combs,metrics)
            meanDf_bucket_aligner=rbind(meanDf_bucket_aligner,tmp)
        }
        # With bucket
        df = rawtable
        combs = c("bucket"=bucket)
        tmp=sumByAlnTree(df,groupfam,combs,metrics)
        meanDf_bucket=rbind(meanDf_bucket,tmp)
    }
    df = rawtable
    tmp=sumByAlnTree(df,groupfam,c(),metrics)
    meanDf=rbind(meanDf,tmp)
}

#print(meanDf)
