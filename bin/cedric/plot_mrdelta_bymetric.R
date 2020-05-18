
# Load libraries
library(argparse)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i","--input",type="character",help="Mrdelta data table. Eg: data/mrdelta/out.aln.tsv")
parser$add_argument("-o","--output",type="character",help="Out png")
#parser$add_argument("-s","--score",type="character",nargs='+')
parser$add_argument("-t","--typ",type="character",default="original",help="Eg. original, normPerLen, normByLen")
parser$add_argument("-b","--bucket",type="integer",nargs='+')
parser$add_argument("-a","--aligner",type="character",nargs='+')
parser$add_argument("-T","--title",type="character")
parser$add_argument("--ymax_tc_delta",type="double",default=12.5)
parser$add_argument("--ymin_tc_delta",type="double",default=-2.5)
parser$add_argument("--ymax_quan_delta",type="double",default=1.25)
parser$add_argument("--ymin_quan_delta",type="double",default=-0.25)
parser$add_argument("--ymax_tc_Nratio",type="double",default=0.8)
parser$add_argument("--ymin_tc_Nratio",type="double",default=0.4)
parser$add_argument("--ymax_quan_Nratio",type="double",default=0.7)
parser$add_argument("--ymin_quan_Nratio",type="double",default=0.4)
parser$add_argument("--ymax_unused",type="double",default=1)
parser$add_argument("--ymin_unused",type="double",default=0)
parser$add_argument("--xmax",type="double",default=1)
parser$add_argument("--xmin",type="double",default=0)
args <- parser$parse_args()


#metrics=c("homo","whomo","whomo2","len","ngap","ngap2")
#score1=c("tc","sp")
#score2=c("quantest30","quantest50","quantest150","quantest1000")
metrics=c("quantest30","quantest50","quantest150","quantest1000")
score1=c("tc")


#### LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
table=read.table(args$input,header=TRUE,sep="\t")

pos=c()
for (typ in args$typ)
{
    for (bucket in args$bucket)
    {
        for (aligner in args$aligner)
        {
            pos=c(pos,which(table$typ==typ & table$bucket==bucket & table$aligner==aligner))
        }
    }
}
tmp=table[pos,]

# TC|SP table
pos=c()
for (score in score1)
{
    pos=c(pos,which(tmp$score==score))
}
table_tc=tmp[pos,]

# # Quantest2 table
# pos=c()
# for (score in score2)
# {
#     pos=c(pos,which(tmp$score==score))
# }
# table_quantest=tmp[pos,]
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




#### FUNCTIONS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
get_color=function(metrics)
{
    colorset=hue_pal()(length(metrics))
    names(colorset)=metrics
    return(colorset)
}


plot_mrdelta_byscore=function(df1,df2)
{
    colorset=get_color(metrics) 

    # TC|SP  
    g_deltaacc_tc=ggplot(df1[which(df1$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+
                facet_wrap(~score,nrow=1)+geom_point(size=0.6)+geom_line()+ylab("Deltaacc")+
                xlim(args$xmin,args$xmax)+ylim(args$ymin_tc_delta,args$ymax_tc_delta)+
                scale_color_manual(values=colorset)
    g_Nratio_tc=ggplot(df1[which(df1$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+
                facet_wrap(~score,nrow=1)+geom_point(size=0.6)+geom_line()+ylab("Nratio")+
                xlim(args$xmin,args$xmax)+ylim(args$ymin_tc_Nratio,args$ymax_tc_Nratio)+
                scale_color_manual(values=colorset)
    g_usedRatio_tc=ggplot(df1[which(df1$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+
                facet_wrap(~score,nrow=1)+geom_point(size=0.6)+geom_line()+ylab("usedRatio")+
                xlim(args$xmin,args$xmax)+ylim(args$ymin_unused,args$ymax_unused)+scale_color_manual(values=colorset)

    # Quantest2
    # g_deltaacc_quan=ggplot(df2[which(df2$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+
    #             facet_wrap(~score,nrow=1)+geom_point(size=0.6)+geom_line()+ylab("Deltaacc")+
    #             xlim(args$xmin,args$xmax)+ylim(args$ymin_quan_delta,args$ymax_quan_delta)+
    #             scale_color_manual(values=colorset)
    # g_Nratio_quan=ggplot(df2[which(df2$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+
    #             facet_wrap(~score,nrow=1)+geom_point(size=0.6)+geom_line()+ylab("Nratio")+
    #             xlim(args$xmin,args$xmax)+ylim(args$ymin_quan_Nratio,args$ymax_quan_Nratio)+
    #             scale_color_manual(values=colorset)
    # g_usedRatio_quan=ggplot(df2[which(df2$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+
    #             facet_wrap(~score,nrow=1)+geom_point(size=0.6)+geom_line()+ylab("usedRatio")+
    #             xlim(args$xmin,args$xmax)+ylim(args$ymin_unused,args$ymax_unused)+scale_color_manual(values=colorset)

    # figure_tc=ggarrange(g_deltaacc_tc,g_Nratio_tc,g_usedRatio_tc,
    #             legend="none",nrow=3,font.label=list(size=12))
    # figure_quan=ggarrange(g_deltaacc_quan,g_Nratio_quan,g_usedRatio_quan,
    #             common.legend=TRUE,legend="bottom",nrow=3,font.label=list(size=12))
    # figure=ggarrange(figure_tc,figure_quan,common.legend=TRUE,legend="bottom",nrow=1,font.label=list(size=12))
    figure=ggarrange(g_deltaacc_tc,g_deltaacc_quan,
            g_Nratio_tc,g_Nratio_quan,
            g_usedRatio_tc,g_usedRatio_quan,
            common.legend=TRUE,legend="bottom",nrow=3,ncol=2,font.label=list(size=12))
    g=annotate_figure(figure,top = text_grob(args$title, size = 14))
    return(g)
}



plot_mrdelta_byscore=function(df1) #,df2)
{
    colorset=get_color(metrics) 

    # DeltaAcc
    g_deltaacc=list()
    for (score in score1)
    {
        pos=which(df1$score==score)
        df=df1[pos,]
        g=ggplot(df[which(df$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+
                geom_point(size=0.6)+geom_line()+ylab("Deltaacc")+labs(title=score)+
                xlim(args$xmin,args$xmax)+ylim(args$ymin_tc_delta,args$ymax_tc_delta)+
                scale_color_manual(values=colorset)
        g_deltaacc[[score]]=g
    }
    # for (score in score2)
    # {
    #     pos=which(df2$score==score)
    #     df=df2[pos,]
    #     g=ggplot(df[which(df$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+
    #             geom_point(size=0.6)+geom_line()+ylab("Deltaacc")+labs(title=score)+
    #             xlim(args$xmin,args$xmax)+ylim(args$ymin_quan_delta,args$ymax_quan_delta)+
    #             scale_color_manual(values=colorset)
    #     g_deltaacc[[score]]=g
    # }

    # Nratio
    g_Nratio=list()
    for (score in score1)
    {
        pos=which(df1$score==score)
        df=df1[pos,]
        g=ggplot(df[which(df$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+
                geom_point(size=0.6)+geom_line()+ylab("Nratio")+labs(title=score)+
                xlim(args$xmin,args$xmax)+ylim(args$ymin_tc_Nratio,args$ymax_tc_Nratio)+
                scale_color_manual(values=colorset) + ylab("Pratio")
        g_Nratio[[score]]=g
    }
    # for(score in score2)
    # {
    #     pos=which(df2$score==score)
    #     df=df2[pos,]
    #     g=ggplot(df[which(df$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+
    #             geom_point(size=0.6)+geom_line()+ylab("Nratio")+labs(title=score)+
    #             xlim(args$xmin,args$xmax)+ylim(args$ymin_quan_Nratio,args$ymax_quan_Nratio)+
    #             scale_color_manual(values=colorset)
    #     g_Nratio[[score]]=g
    # }

    # UsedRatio
    g_usedRatio=list()
    for (score in score1)
    {
        pos=which(df1$score==score)
        df=df1[pos,]
        g=ggplot(df[which(df$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+
                geom_point(size=0.6)+geom_line()+ylab("usedRatio")+labs(title=score)+
                xlim(args$xmin,args$xmax)+ylim(args$ymin_unused,args$ymax_unused)+scale_color_manual(values=colorset)
        g_usedRatio[[score]]=g
    }
    # for (score in score2)
    # {
    #     pos=which(df2$score==score)
    #     df=df2[pos,]
    #     g=ggplot(df[which(df$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+
    #             geom_point(size=0.6)+geom_line()+ylab("usedRatio")+labs(title=score)+
    #             xlim(args$xmin,args$xmax)+ylim(args$ymin_unused,args$ymax_unused)+scale_color_manual(values=colorset)
    #     g_usedRatio[[score]]=g
    # }



    # Merge figures
    figure_deltaacc=ggarrange(plotlist=g_deltaacc,font.label=list(size=12),nrow=1,legend="none")
    figure_nratio=ggarrange(plotlist=g_Nratio,font.label=list(size=12),nrow=1,legend="none")
    figure_usedratio=ggarrange(plotlist=g_usedRatio,
                     common.legend=TRUE,legend="right",font.label=list(size=12),nrow=1)
    figure=ggarrange(figure_deltaacc,figure_nratio,figure_usedratio,nrow=1,
                    common.legend=TRUE,legend="right",font.label=list(size=12))
    g=annotate_figure(figure,top = text_grob(args$title, size = 14))
    return(g)
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



#### Plot Nratio/deltaacc/usedRatio vs mrdelta
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df_tc=data.frame(value=c(as.numeric(as.character(table_tc$Nratio)),
                           as.numeric(as.character(table_tc$deltaacc)),
                           1-as.numeric(as.character(table_tc$unusedRatio))),
                    mrdelta=rep(as.numeric(as.character(table_tc$mrdelta)),3),
                    metric=rep(as.character(table_tc$metric),3),
                    score=rep(as.character(table_tc$score),3),
                    typ=c(rep("Nratio",nrow(table_tc)),rep("deltaacc",nrow(table_tc)),rep("usedRatio",nrow(table_tc)))
                   )
# df_quantest=data.frame(value=c(as.numeric(as.character(table_quantest$Nratio)),
#                            as.numeric(as.character(table_quantest$deltaacc)),
#                            1-as.numeric(as.character(table_quantest$unusedRatio))),
#                     mrdelta=rep(as.numeric(as.character(table_quantest$mrdelta)),3),
#                     metric=rep(as.character(table_quantest$metric),3),
#                     score=rep(as.character(table_quantest$score),3),
#                     typ=c(rep("Nratio",nrow(table_quantest)),rep("deltaacc",nrow(table_quantest)),rep("usedRatio",nrow(table_quantest)))
#                    )

df_tc$metric=factor(df_tc$metric,level=metrics)

g=plot_mrdelta_byscore(df_tc) #,df_quantest)
ggsave(args$output,plot=g,width=16,height=8)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<