library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)

ori="/users/cn/sjin/projects/homoplasy/nf_homoplasty"
foldname="/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta_quantest/global/ylim"  # Output folder

#ylim_min=-10
#ylim_min=-5
ylim_min=-2
ylim_max=2.5
#ylim_max=15
xlim_min=0
xlim_max=1

buckets=c(50,100,200,500,1000)
metrics=c("homo","whomo","whomo2","len","ngap","ngap2")


args = commandArgs(trailingOnly=TRUE)
score=args[1]  # tc | sp
typ=args[2]    # original | normPerLen | normByLen
fil=args[3]


# Read file
#df=read.table("data/mrdelta/out.global.tsv",header=TRUE,sep="\t")
df=read.table(fil,header=TRUE,sep="\t")
pos=which(df$score==score&df$typ==typ)
df2=df[pos,]


# Functions
get_color=function(metrics)
{
    colorset=hue_pal()(length(metrics))
    names(colorset)=metrics
    return(colorset)
}


plot_vsmrdelta=function(df2)
{
    colorset=get_color(metrics)   
    g_deltaacc=ggplot(df2[which(df2$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(~bucket,nrow=1)+
               geom_point(size=0.6)+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    g_Nratio=ggplot(df2[which(df2$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(~bucket,nrow=1)+
               geom_point(size=0.6)+geom_line()+ylab("Nratio")+ylim(0.4,0.8)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    g_usedRatio=ggplot(df2[which(df2$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(~bucket,nrow=1)+
               geom_point(size=0.6)+geom_line()+ylab("usedRatio")+ylim(0,1)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    figure=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,legend="bottom",nrow=3,font.label=list(size=12))
    g=annotate_figure(figure,top = text_grob(paste(toupper(score)," - ",typ,sep=""), size = 14))
    return(g)
}


plot_vsunused=function(df2)
{
    colorset=get_color(metrics)
    g_deltaacc=ggplot(df2[which(df2$typ=="deltaacc"),],aes(x=unusedRatio,y=value,color=metric))+facet_wrap(~bucket,nrow=1)+
               geom_point(size=0.6)+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    g_Nratio=ggplot(df2[which(df2$typ=="Nratio"),],aes(x=unusedRatio,y=value,color=metric))+facet_wrap(~bucket,nrow=1)+
               geom_point(size=0.6)+geom_line()+ylab("Nratio")+ylim(0.4,1)+xlim(xlim_min,xlim_max)+scale_color_manual(values=colorset)
    figure=ggarrange(g_deltaacc,g_Nratio,common.legend=TRUE,legend="bottom",nrow=2,font.label=list(size=12))
    g=annotate_figure(figure,top = text_grob(paste(toupper(score)," - ",typ,sep=""), size = 14))
    return(g)
}



## =========================================
## Plot Nratio/deltaacc/usedRatio vs mrdelta

df_plot=data.frame(value=c(as.numeric(as.character(df2$Nratio)),
                           as.numeric(as.character(df2$deltaacc)),
                           1-as.numeric(as.character(df2$unusedRatio))),
                    mrdelta=rep(as.numeric(as.character(df2$mrdelta)),3),
                    metric=rep(as.character(df2$metric),3),
                    bucket=rep(as.factor(as.numeric(as.character(df2$bucket))),3),
                    typ=c(rep("Nratio",nrow(df2)),rep("deltaacc",nrow(df2)),rep("usedRatio",nrow(df2)))
                   )
g=plot_vsmrdelta(df_plot)
ggsave(paste(foldname,"/mrdelta.",score,".",typ,".global.png",sep=""),plot=g,width=14.5,height=8)




## ==============================
## Nratio/deltaacc vs unusedRatio

df_plot=data.frame(value=c(as.numeric(as.character(df2$Nratio)),as.numeric(as.character(df2$deltaacc))),
                   mrdelta=rep(as.numeric(as.character(df2$mrdelta)),2),
                   score=rep(as.character(df2$score),2),
                   metric=rep(as.character(df2$metric),2),
                   unusedRatio=rep(as.numeric(as.character(df2$unusedRatio)),2),
                   bucket=rep(as.factor(as.numeric(as.character(df2$bucket))),2),
                   typ=c(rep("Nratio",nrow(df2)),rep("deltaacc",nrow(df2)))
                   )

g=plot_vsunused(df_plot)
ggsave(paste(foldname,"/mrdelta.unused.",score,".",typ,".global.png",sep=""),plot=g,width=14,height=7.2)




