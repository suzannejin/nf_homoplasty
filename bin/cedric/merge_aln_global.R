library(argparse)

parser <- ArgumentParser()
parser$add_argument("--in_aln",type="character",help="Mrdelta data table - aln")
parser$add_argument("--in_global",type="character",help="Mrdelta data table - global")
parser$add_argument("--out",type="character",help="Merged table")
args <- parser$parse_args()


t1=read.table(args$in_aln,header=TRUE,sep="\t")
t2=read.table(args$in_global,header=TRUE,sep="\t")


t=data.frame(mrdelta=c(as.character(t1$mrdelta),as.character(t2$mrdelta)),
             score=c(as.character(t1$score),as.character(t2$score)),
             metric=c(as.character(t1$metric),as.character(t2$metric)),
             typ=c(as.character(t1$typ),as.character(t2$typ)),
             aligner=c(as.character(t1$aligner),as.character(t2$aligner)),
             bucket=c(as.character(t1$bucket),as.character(t2$bucket)),
             Nratio=c(as.character(t1$Nratio),as.character(t2$Nratio)),
             unusedRatio=c(as.character(t1$unusedRatio),as.character(t2$unusedRatio)),
             minacc=c(as.character(t1$minacc),as.character(t2$minacc)),
             maxacc=c(as.character(t1$maxacc),as.character(t2$maxacc)),
             deltaacc=c(as.character(t1$deltaacc),as.character(t2$deltaacc)))

write.table(t,args$out,quote=FALSE,row.names=FALSE,sep="\t")