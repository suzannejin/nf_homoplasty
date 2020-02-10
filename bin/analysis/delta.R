
## Libraries & packages
#install.packages('ggplot2')
library(ggplot2)
#install.packages('plotly')
library(plotly)
#install.packages('gridExtra')
library(gridExtra)
#install.packages('ggpubr')
library(ggpubr)


metrics=c("tc","homo","w_homo","w_homo2","ngap","ngap2")


# Import data
fil=read.csv("delta_data.csv",header=TRUE,sep=" ",row.names=NULL)

#fil=read.csv("delta_viceversa_data.csv",header=TRUE,sep=" ",row.names=NULL)
fil$tree=paste(fil$tree1,fil$tree2,sep="_")



# Remove guide trees
rm_tree=function(df,treename)
{
    pos=c(which(df$tree1==treename))
    pos=c(pos,which(df$tree2==treename))
    df2=df[-pos,]
    return(df2)
}

to_be_rm=c("parttreednd2","parttreednd2size")
for (tr in to_be_rm)
{
    fil=rm_tree(fil,tr)
}




# Plot function : tc vs metric
gplot_tc_metric=function(df,xs,color,shape,title){
    plotList=list()
    n=length(xs)
    for (i in 1:n){
        plotList[[i]]=ggplot(df,aes_string(x=xs[i],y="tc",color=color,shape=shape)) +
            geom_point() +
            #geom_smooth(method = "lm", fill = NA)+
            geom_vline(xintercept = 0) + geom_hline(yintercept = 0)+
            ggtitle(title) +
            xlab(xs[i]) + ylab("") 
    }
    names(plotList)=xs
    return(plotList)
}



# Plot all in grid
plotList=gplot_tc_metric(fil,metrics[2:6],"tree1","aligner","")
figure=ggarrange(plotList$ngap,plotList$ngap2,plotList$homo,plotList$w_homo,plotList$w_homo2,common.legend=TRUE)

annotate_figure(figure,
                top = text_grob(paste("All delta together",sep=""), size = 14),
                left = text_grob("TC Score", rot = 90)
)
ggsave("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/plots_12trees_above1000/tc_metric.grid.above1000.3aligner.12tree.png") 



# Plot dataset by dataset in grid
fams=levels(fil$family)
for (fam in fams){
    df2=fil[which(fil$family==fam),]
    plotList=gplot_tc_metric(df2,metrics[2:6],"tree1","aligner","")
    figure=ggarrange(plotList$ngap,plotList$ngap2,plotList$homo,plotList$w_homo,plotList$w_homo2,common.legend=TRUE)
    annotate_figure(figure,
                top = text_grob(paste("Family: ",fam,sep=""), size = 14),
                left = text_grob("TC Score", rot = 90)
    )
    out=paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/plots_12trees_above1000/tc_metric.grid.",fam,".3aligner.12tree.png",sep="")
    ggsave(out) 
}



# Plot dataset by dataset & aligner in grid
fams=levels(fil$family)
aligners=levels(fil$aligner)
for (fam in fams){
    for (aln in aligners){
        df2=fil[which(fil$family==fam & fil$aligner==aln),]
        plotList=gplot_tc_metric(df2,metrics[2:6],"tree1","aligner","")
        figure=ggarrange(plotList$ngap,plotList$ngap2,plotList$homo,plotList$w_homo,plotList$w_homo2,common.legend=TRUE)
        annotate_figure(figure,
                    top = text_grob(paste("Family: ",fam,sep=""), size = 14),
                    left = text_grob("TC Score", rot = 90)
        )
        out=paste("/users/cn/sjin/projects/homoplasy/nf_homoplasty/plots/plots_12trees_above1000/tc_metric.grid.",fam,".",aln,"12tree.png",sep="")
        ggsave(out) 
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

func_w_homo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$w_homo)))
}

func_w_homo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$w_homo2)))
}

func_w_homo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$ngap)))
}

func_w_homo=function(xx)
{
    return(data.frame(COR = cor(xx$tc, xx$ngap2)))
}
###



