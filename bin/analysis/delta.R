
## Libraries & packages
#install.packages('ggplot2')
library(ggplot2)
#install.packages('plotly')
library(plotly)
#install.packages('gridExtra')
library(gridExtra)
#install.packages('ggpubr')
library(ggpubr)


# Remember to order levels tree.
# And plot the same plots but without "parttreednd2 parttreednd2size"


args = commandArgs(trailingOnly=TRUE)
bucket=args[1] # 50, 100, 200, 500, 1000
typ=args[2]  # original | normPerLen | normByLen
rm_t=args[3] # yes | no
treen=args[4]  # tree1 | tree2

# bucket=100
# typ="original"
# rm_t="yes"
# treen="tree1"

# for (bucket in buckets){for (typ in typs){for (rm_t in rm_ts){for (treen in treens){source("bin/analysis/delta.R")}}}}

metrics=c("tc","homo","whomo","whomo2","len","ngap","ngap2")
tree_names=c("codnd","dpparttreednd1","dpparttreednd2","dpparttreednd2size","fastaparttreednd","fftns1dnd","fftns1dndmem","fftns2dnd","fftns2dndmem","mafftdnd","parttreednd0","parttreednd1","parttreednd2","parttreednd2size","FAMSA","CLUSTALO-RANDOM")



# Import data
f=paste("data/delta/bucket",bucket,"/delta.",typ,".csv",sep="")
fil=read.csv(f,header=TRUE,sep=",",row.names=NULL)

fil$tree1=factor(fil$tree1,levels=tree_names[1:length(tree_names)-1])
fil$tree2=factor(fil$tree2,levels=tree_names[2:length(tree_names)])

#fil=read.csv("delta_viceversa_data.csv",header=TRUE,sep=" ",row.names=NULL)
fil$tree=paste(fil$tree1,fil$tree2,sep=".")



# Remove guide trees
rm_tree=function(df,treename)
{
    pos=c(which(df$tree1==treename))
    pos=c(pos,which(df$tree2==treename))
    df2=df[-pos,]
    return(df2)
}

if (rm_t == "yes")
{
    to_be_rm=c("parttreednd2","parttreednd2size","CLUSTALO-RANDOM")
    for (tr in to_be_rm)
    {
        fil=rm_tree(fil,tr)
    }
}


# gg color
gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



# Plot function : tc vs metric
gplot_tc_metric=function(df,xs,color,shape,title){
    c_trees_names=unique(c(as.character(unique(fil$tree1)),as.character(unique(fil$tree2))))
    c_trees_color=gg_color(length(c_trees_names))
    if (color=="tree1")
    {
        labels=c_trees_names[1:(length(c_trees_names)-1)]
        colors=c_trees_color[1:(length(c_trees_color)-1)]
    }
    else
    {
        labels=c_trees_names[2:length(c_trees_names)]
        colors=c_trees_color[2:length(c_trees_color)]
    }
    plotList=list()
    n=length(xs)
    for (i in 1:n){
        plotList[[i]]=ggplot(df,aes_string(x=xs[i],y="tc",color=color,shape=shape)) +
            geom_point() +
            #geom_smooth(method = "lm", fill = NA)+
            geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
            ggtitle(title) +
            xlab(paste("Δ",xs[i],sep="")) + ylab("")  +
            ylim(-100,100) +
            scale_color_manual(values=colors,labels=labels)
    }
    names(plotList)=xs
    return(plotList)
}



# Plot all in grid
out=paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/delta/bucket",bucket,"/delta.",typ,".combinedSeqs.allaligner.",rm_t,"rm.",treen,".bucket",bucket,".png",sep="")
if (!file.exists(out))
{
    plotList=gplot_tc_metric(fil,metrics[2:7],treen,"aligner","")
    figure=ggarrange(plotList$ngap,plotList$ngap2,plotList$len,plotList$homo,plotList$whomo,plotList$whomo2,common.legend=TRUE)

    g=annotate_figure(figure,
                    top = text_grob(paste("All delta together - ",typ," - bucket",bucket,sep=""), size = 14),
                    left = text_grob("Δ TC Score", rot = 90)
    )
    #ggsave(out)
    ggsave(out,plot=g,width=14.2,height=7.34)
}





# Plot dataset by dataset in grid
fams=levels(fil$family)
for (fam in fams)
{
    out=paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/delta/bucket",bucket,"/delta.",typ,".fam.",fam,".combinedSeqs.allaligner.",rm_t,"rm.",treen,".bucket",bucket,".png",sep="")
    if (!file.exists(out))
    {
        df2=fil[which(fil$family==fam),]
        plotList=gplot_tc_metric(df2,metrics[2:7],treen,"aligner","")
        figure=ggarrange(plotList$ngap,plotList$ngap2,plotList$len,plotList$homo,plotList$whomo,plotList$whomo2,common.legend=TRUE)
        g=annotate_figure(figure,
                    top = text_grob(paste("Family: ",fam," - ",typ," - bucket",bucket,sep=""), size = 14),
                    left = text_grob("Δ TC Score", rot = 90)
        )
        #ggsave(out)
        ggsave(out,plot=g,width=14.2,height=7.34)
    } 

}



# Plot dataset by dataset & aligner in grid
fams=levels(fil$family)
aligners=levels(fil$aligner)
for (fam in fams){
    for (aln in aligners){
        out=paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/delta/bucket",bucket,"/delta.",typ,".fam.",fam,".combinedSeqs.",aln,".",rm_t,"rm.",treen,".bucket",bucket,".png",sep="")
        if (!file.exists(out))
        {
            df2=fil[which(fil$family==fam & fil$aligner==aln),]
            plotList=gplot_tc_metric(df2,metrics[2:7],treen,"aligner","")
            figure=ggarrange(plotList$ngap,plotList$ngap2,plotList$len,plotList$homo,plotList$whomo,plotList$whomo2,common.legend=TRUE)
            g=annotate_figure(figure,
                        top = text_grob(paste("Family: ",fam," - ",typ," - bucket",bucket,sep=""), size = 14),
                        left = text_grob("Δ TC Score", rot = 90)
            )
            #ggsave(out)
            ggsave(out,plot=g,width=14.2,height=7.34)
        }
    }
}







###
##  Function to perform CORRELATION by GROUP
###
require(plyr)

func_homo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$homo)))
}

func_whomo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$whomo)))
}

func_whomo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$whomo2)))
}

func_whomo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$ngap)))
}

func_whomo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$ngap2)))
}
###



