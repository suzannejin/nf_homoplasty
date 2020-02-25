library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)



metrics=c("homo","whomo","whomo2","ngap","ngap2",
            "homoPerLen","whomoPerLen","whomo2PerLen","ngapPerLen","ngap2PerLen",
            "homoByLen","whomoByLen","whomo2ByLen","ngapByLen","ngap2ByLen",
            "homoPerSeq","whomoPerSeq","whomo2PerSeq","ngapPerSeq","ngap2PerSeq",
            "homoBySeq","whomoBySeq","whomo2BySeq","ngapBySeq","ngap2BySeq",
            "homoPerLenSeq","whomoPerLenSeq","whomo2PerLenSeq","ngapPerLenSeq","ngap2PerLenSeq",
            "homoByLenSeq","whomoByLenSeq","whomo2ByLenSeq","ngapByLenSeq","ngap2ByLenSeq")
metrics_original=c("homo","whomo","whomo2","ngap","ngap2")
metrics_perlen=c("homoPerLen","whomoPerLen","whomo2PerLen","ngapPerLen","ngap2PerLen")
metrics_bylen=c("homoByLen","whomoByLen","whomo2ByLen","ngapByLen","ngap2ByLen")
metrics_perseq=c("homoPerSeq","whomoPerSeq","whomo2PerSeq","ngapPerSeq","ngap2PerSeq")
metrics_byseq=c("homoBySeq","whomoBySeq","whomo2BySeq","ngapBySeq","ngap2BySeq")
metrics_perlenseq=c("homoPerLenSeq","whomoPerLenSeq","whomo2PerLenSeq","ngapPerLenSeq","ngap2PerLenSeq")
metrics_bylenseq=c("homoByLenSeq","whomoByLenSeq","whomo2ByLenSeq","ngapByLenSeq","ngap2ByLenSeq")

ylim_min=-5
ylim_max=15
foldname="global/bucket100_1000"



f_global_100=read.table("/users/cn/sjin/projects/homoplasy/nf_homoplasty/data/mrdelta/out.bucket100.global.tsv",header=TRUE)
f_global_1000=read.table("/users/cn/sjin/projects/homoplasy/nf_homoplasty/data/mrdelta/out.bucket1000.global.tsv",header=TRUE)

f_global_100$metric=factor(f_global_100$metric,levels=metrics,ordered=TRUE)
f_global_1000$metric=factor(f_global_1000$metric,levels=metrics,ordered=TRUE)



df_global_100=data.frame(value=c(f_global_100$Nratio,f_global_100$deltaacc,1-f_global_100$unusedRatio),
                         mrdelta=rep(f_global_100$mrdelta,3),
                         score=rep(as.character(f_global_100$score),3),
                         metric=rep(as.character(f_global_100$metric),3),
                         bucket=rep("bucket100",3*nrow(f_global_100)),
                         typ=c(rep("Nratio",nrow(f_global_100)),rep("deltaacc",nrow(f_global_100)),rep("usedRatio",nrow(f_global_100)))
                        )
df_global_1000=data.frame(value=c(f_global_1000$Nratio,f_global_1000$deltaacc,1-f_global_1000$unusedRatio),
                         mrdelta=rep(f_global_1000$mrdelta,3),
                         score=rep(as.character(f_global_1000$score),3),
                         metric=rep(as.character(f_global_1000$metric),3),
                         bucket=rep("bucket1000",3*nrow(f_global_1000)),
                         typ=c(rep("Nratio",nrow(f_global_1000)),rep("deltaacc",nrow(f_global_1000)),rep("usedRatio",nrow(f_global_1000)))
                        )
df_global=rbind(df_global_100,df_global_1000)


get_color=function(metrics)
{
    colorset=hue_pal()(length(metrics))
    names(colorset)=metrics
    return(colorset)
}


# Original metrics
df_global_original=df_global[which(df_global$metric %in% metrics_original),]
# Metrics normalized per length
df_global_perlen=df_global[which(df_global$metric %in% metrics_perlen),]
# Metrics normalized by length
df_global_bylen=df_global[which(df_global$metric %in% metrics_bylen),]
# Metrics normalized per nseq
df_global_perseq=df_global[which(df_global$metric %in% metrics_perseq),]
# Metrics normalized by nseq
df_global_byseq=df_global[which(df_global$metric %in% metrics_byseq),]
# Metrics normalized per length & nseq
df_global_perlenseq=df_global[which(df_global$metric %in% metrics_perlenseq),]
# Metrics normalized by length & nseq
df_global_bylenseq=df_global[which(df_global$metric %in% metrics_bylenseq),]


# Original metrics
colorset=get_color(metrics_original)
g_deltaacc=ggplot(df_global_original[which(df_global_original$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+scale_color_manual(values=colorset)
g_Nratio=ggplot(df_global_original[which(df_global_original$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+scale_color_manual(values=colorset)
g_usedRatio=ggplot(df_global_original[which(df_global_original$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("usedRatio")+ylim(0,1)+scale_color_manual(values=colorset)
g=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,nrow=3)
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.originalmetrics.global.png",sep=""),plot=g,width=17.5,height=10.5)
# Using metrics normalized per length
colorset=get_color(metrics_perlen)
g_deltaacc=ggplot(df_global_perlen[which(df_global_perlen$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+scale_color_manual(values=colorset)
g_Nratio=ggplot(df_global_perlen[which(df_global_perlen$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+scale_color_manual(values=colorset)
g_usedRatio=ggplot(df_global_perlen[which(df_global_perlen$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("usedRatio")+ylim(0,1)+scale_color_manual(values=colorset)
g=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,nrow=3)
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.metricperlen.global.png",sep=""),plot=g,width=17.5,height=10.5)
# Using metrics multiplied by length
colorset=get_color(metrics_bylen)
g_deltaacc=ggplot(df_global_bylen[which(df_global_bylen$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+scale_color_manual(values=colorset)
g_Nratio=ggplot(df_global_bylen[which(df_global_bylen$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+scale_color_manual(values=colorset)
g_usedRatio=ggplot(df_global_bylen[which(df_global_bylen$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("usedRatio")+ylim(0,1)+scale_color_manual(values=colorset)
g=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,nrow=3)
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.metricbylen.global.png",sep=""),plot=g,width=17.5,height=10.5)
# Using metrics normalized per nseq
colorset=get_color(metrics_perseq)
g_deltaacc=ggplot(df_global_perseq[which(df_global_perseq$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+scale_color_manual(values=colorset)
g_Nratio=ggplot(df_global_perseq[which(df_global_perseq$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+scale_color_manual(values=colorset)
g_usedRatio=ggplot(df_global_perseq[which(df_global_perseq$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("usedRatio")+ylim(0,1)+scale_color_manual(values=colorset)
g=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,nrow=3)
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.metricperseq.global.png",sep=""),plot=g,width=17.5,height=10.5)
# Using metrics multiplied by nseq
colorset=get_color(metrics_byseq)
g_deltaacc=ggplot(df_global_byseq[which(df_global_byseq$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+scale_color_manual(values=colorset)
g_Nratio=ggplot(df_global_byseq[which(df_global_byseq$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+scale_color_manual(values=colorset)
g_usedRatio=ggplot(df_global_byseq[which(df_global_byseq$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("usedRatio")+ylim(0,1)+scale_color_manual(values=colorset)
g=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,nrow=3)
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.metricbyseq.global.png",sep=""),plot=g,width=17.5,height=10.5)
# Using metrics normalized per length & nseq
colorset=get_color(metrics_perlenseq)
g_deltaacc=ggplot(df_global_perlenseq[which(df_global_perlenseq$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+scale_color_manual(values=colorset)
g_Nratio=ggplot(df_global_perlenseq[which(df_global_perlenseq$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+scale_color_manual(values=colorset)
g_usedRatio=ggplot(df_global_perlenseq[which(df_global_perlenseq$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("usedRatio")+ylim(0,1)+scale_color_manual(values=colorset)
g=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,nrow=3)
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.metricperlenseq.global.png",sep=""),plot=g,width=17.5,height=10.5)
# Using metrics multiplied by length & nseq
colorset=get_color(metrics_bylenseq)
g_deltaacc=ggplot(df_global_bylenseq[which(df_global_bylenseq$typ=="deltaacc"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Deltaacc")+ylim(ylim_min,ylim_max)+scale_color_manual(values=colorset)
g_Nratio=ggplot(df_global_bylenseq[which(df_global_bylenseq$typ=="Nratio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("Nratio")+ylim(0,1)+scale_color_manual(values=colorset)
g_usedRatio=ggplot(df_global_bylenseq[which(df_global_bylenseq$typ=="usedRatio"),],aes(x=mrdelta,y=value,color=metric))+facet_wrap(score~bucket,nrow=1)+geom_point()+geom_line()+ylab("usedRatio")+ylim(0,1)+scale_color_manual(values=colorset)
g=ggarrange(g_deltaacc,g_Nratio,g_usedRatio,common.legend=TRUE,nrow=3)
ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.metricbylenseq.global.png",sep=""),plot=g,width=17.5,height=10.5)









# # All metrics
# df_global_sp=df_global[which(df_global$score=="sp"),]
# df_global_tc=df_global[which(df_global$score=="tc"),]
# # Original metrics
# df_global_sp_original=df_global_sp[which(df_global_sp$metric %in% metrics_original),]
# df_global_tc_original=df_global_tc[which(df_global_tc$metric %in% metrics_original),]
# # Metrics normalized per length
# df_global_sp_perlen=df_global_sp[which(df_global_sp$metric %in% metrics_perlen),]
# df_global_tc_perlen=df_global_tc[which(df_global_tc$metric %in% metrics_perlen),]
# # Metrics normalized by length
# df_global_sp_bylen=df_global_sp[which(df_global_sp$metric %in% metrics_bylen),]
# df_global_tc_bylen=df_global_tc[which(df_global_tc$metric %in% metrics_bylen),]
# # Metrics normalized per nseq
# df_global_sp_perseq=df_global_sp[which(df_global_sp$metric %in% metrics_perseq),]
# df_global_tc_perseq=df_global_tc[which(df_global_tc$metric %in% metrics_perseq),]
# # Metrics normalized by nseq
# df_global_sp_byseq=df_global_sp[which(df_global_sp$metric %in% metrics_byseq),]
# df_global_tc_byseq=df_global_tc[which(df_global_tc$metric %in% metrics_byseq),]
# # Metrics normalized per length & nseq
# df_global_sp_perlenseq=df_global_sp[which(df_global_sp$metric %in% metrics_perlenseq),]
# df_global_tc_perlenseq=df_global_tc[which(df_global_tc$metric %in% metrics_perlenseq),]
# # Metrics normalized by length & nseq
# df_global_sp_bylenseq=df_global_sp[which(df_global_sp$metric %in% metrics_bylenseq),]
# df_global_tc_bylenseq=df_global_tc[which(df_global_tc$metric %in% metrics_bylenseq),]


# # Plot
# # Original metrics
# g_global_sp_original=ggplot(df_global_sp_original,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_tc_original=ggplot(df_global_tc_original,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_ratio_sp_original
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.sp.originalmetrics.global.png",sep=""))
# g_global_ratio_tc_original
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.tc.originalmetrics.global.png",sep=""))
# # Using metrics normalized per length
# g_global_sp_perlen=ggplot(df_global_sp_perlen,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_tc_perlen=ggplot(df_global_tc_perlen,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_sp_perlen
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.sp.metricperlen.global.png",sep=""))
# g_global_tc_perlen
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.tc.metricperlen.global.png",sep=""))
# # Using metrics multiplied by length
# g_global_sp_bylen=ggplot(df_global_sp_bylen,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_tc_bylen=ggplot(df_global_tc_bylen,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_sp_bylen
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.sp.metricbylen.global.png",sep=""))
# g_global_tc_bylen
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.tc.metricbylen.global.png",sep=""))
# # Using metrics normalized per nseq
# g_global_sp_perseq=ggplot(df_global_sp_perseq,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_tc_perseq=ggplot(df_global_tc_perseq,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_sp_perseq
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.sp.metricperseq.global.png",sep=""))
# g_global_tc_perseq
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.tc.metricperseq.global.png",sep=""))
# # Using metrics multiplied by nseq
# g_global_sp_byseq=ggplot(df_global_sp_byseq,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_tc_byseq=ggplot(df_global_tc_byseq,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_sp_byseq
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.sp.metricbyseq.global.png",sep=""))
# g_global_tc_byseq
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.tc.metricbyseq.global.png",sep=""))
# # Using metrics normalized per length & nseq
# g_global_sp_perlenseq=ggplot(df_global_sp_perlenseq,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_tc_perlenseq=ggplot(df_global_tc_perlenseq,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_sp_perlenseq
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.sp.metricperlenseq.global.png",sep=""))
# g_global_tc_perlenseq
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.tc.metricperlenseq.global.png",sep=""))
# # Using metrics multiplied by length & nseq
# g_global_sp_bylenseq=ggplot(df_global_sp_bylenseq,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_tc_bylenseq=ggplot(df_global_tc_bylenseq,aes(x=mrdelta,y=value,color=metric))+
#                                     facet_wrap(bucket~typ)+geom_point()+geom_line()+
#                                     scale_y_continuous("deltaacc",limits=c(ylim_min,ylim_max),sec.axis=sec_axis(~ rescale(.,to=c(0,1),from=c(ylim_min,ylim_max)),name="Ratio"))
# g_global_sp_bylenseq
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.sp.metricbylenseq.global.png",sep=""))
# g_global_tc_bylenseq
# ggsave(paste("/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/plots/mrdelta/",foldname,"/mrdelta.tc.metricbylenseq.global.png",sep=""))
