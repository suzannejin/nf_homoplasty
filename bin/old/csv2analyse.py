#!/usr/bin/env python3

__description__='''

Python version of the Cedric's script for the analysis of homoplasy.

'''

def csv2pandas(fil):
    import pandas as pd
    df=pd.read_csv(fil,header=[0,1],index_col=0)
    return(df)

def pandas2headers(df):
    # Get headers: families, aligners, trees
    families=list(df.index)
    tmp=[x[0] for x in df["nseq"].values.tolist()]
    nseqs={}
    for i in range(len(families)):
        fam=families[i]
        n=tmp[i]
        nseqs[fam]=n
    aligners_rep=[x[0] for x in list(df.columns)][1:]  # Note that they are repeated (just as in the csv file)
    trees_rep=[x[1] for x in list(df.columns)][1:]
    aligners,trees=[],[]
    for x in aligners_rep:
        if x not in aligners:
            aligners.append(x)
    for x in trees_rep:
        if x not in trees:
            trees.append(x)
    return(families,nseqs,aligners,trees)

def csv2analyse(metricdf,scoredf,families,nseqs,aligners,trees,minseq,maxseq,mrdelta):
    import numpy as np
    import pandas as pd

    # Initialize numbers
    nfam=len(families)
    naln=len(aligners)
    unused=pd.DataFrame(np.zeros( (naln,nfam)), index=aligners, columns=families )
    total=pd.DataFrame(np.zeros( (naln,nfam)) ,index=aligners, columns=families )
    positiv=pd.DataFrame(np.zeros( (naln,nfam)) ,index=aligners, columns=families )
    negativ=pd.DataFrame(np.zeros( (naln,nfam)) ,index=aligners, columns=families )
    pick=pd.DataFrame(np.zeros( (naln,4)) ,index=aligners, columns=["minS","maxS","minN","maxN"] )

    # Compute
    for aligner in aligners:
        for family in families:
            for tr1 in trees:
                for tr2 in trees:
                    total.at[aligner,family]+=1  # total +1
                    use=0
                    # Metrics
                    h1=metricdf[aligner][tr1][family]
                    h2=metricdf[aligner][tr2][family]
                    # Score
                    s1=scoredf[aligner][tr1][family]
                    s2=scoredf[aligner][tr2][family]
                    if h1=="NA" or h2=="NA" or s1=="NA" or s2=="NA":
                        unused.at[aligner,family]+=1
                    elif nseqs[family]<minseq or nseqs[family]>maxseq:
                        unused.at[aligner,family]+=1
                    elif h1==h2 and h1==0:  # If both metrics are 0
                        unused.at[aligner,family]+=1
                    elif abs(h1-h2)/(h1+h2)<mrdelta:
                        unused.at[aligner,family]+=1
                    elif abs(h1-h2)<=0.001 and abs(h1-h2)<=0.001:
                        unused.at[aligner,family]+=1
                    else:
                        d1=h1-h2
                        d2=s1-s2
                        if (d2>=0 and d1>=0) or (d2<=0 and d1<=0):
                            use=1
                            positiv.at[aligner,family]+=1
                        elif (d2>0 and d1<0) or (d2<0 and d1>0):
                            use=1
                            negativ.at[aligner,family]+=1
                    if use>0:
                        if h1>h2:
                            pick.at[aligner,"minS"]+=s2
                            pick.at[aligner,"maxS"]+=s1
                        else:
                            pick.at[aligner,"minS"]+=s1
                            pick.at[aligner,"maxS"]+=s2
                        pick.at[aligner,"minN"]+=1
                        pick.at[aligner,"maxN"]+=1
        
        N=negativ.loc[aligner,family]
        t=total.loc[aligner,family]
        if t!=0:
            Nratio=N/t
        else:
            Nratio=0
        u=unused.loc[aligner,family]
        if u!=0:
            unusedRatio=u/t
        else:
            unusedRatio=0
        sys.stdout.write("FAM::{} [{}]: {:.3f} ({}/{}) -- UNUSED: {:.3f} ({}/{})\n".format(family,aligner,Nratio,N,t,unusedRatio,u,t))




if __name__ == '__main__':

    import sys
    import argparse

    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("-dir",type=str,help="Directory where the input (tc,sp,homo,whomo,whomo2,ngap,ngap2).csv files are stored.")
    app.add_argument("-score",type=str,choices=["tc","sp"])
    app.add_argument("-metrics",type=str,choices=["homo","whomo","whomo2","ngap","ngap2"])
    app.add_argument("-norm",type=str,choices=["PerLen","PerSeq","PerLenSeq","ByLen","BySeq","ByLenSeq"],default=None,help="Data normalization. \
                    Divided by length, number of sequence, or both. Or multiplied by Length, number of sequence, or both.")
    app.add_argument("-mrdelta",type=float,default=0)
    app.add_argument("-minseq",type=int,default=0)
    app.add_argument("-maxseq",type=int,default=100000)
    args=app.parse_args()

    # Read csv to panda Data Frames
    if args.norm:
        metricfile=args.dir+"/"+args.metrics+".norm"+args.norm+".csv"
    else:
        metricfile=args.dir+"/"+args.metrics+".csv"
    scorefile=args.dir+"/"+args.score+".csv"
    metricdf=csv2pandas(metricfile)
    scoredf=csv2pandas(scorefile)
    # Headers & index
    families,nseqs,aligners,trees=pandas2headers(metricdf)

    csv2analyse(metricdf,scoredf,families,nseqs,aligners,trees,args.minseq,args.maxseq,args.mrdelta)



