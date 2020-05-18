import sys
import pandas as pd
import argparse


app = argparse.ArgumentParser()
app.add_argument('folder',type=str,help="Folder where the {tc,sp,homo,...}.csv files are stored")
app.add_argument('bucket',type=int,help="Bucket size")
app.add_argument('nseqfile',type=str,help="path to combinedSeqs_nseq file")
app.add_argument('typ',type=str,choices=["original","normPerLen","normByLen"],default=None)
args=app.parse_args()



metrics=["tc","sp","col","homo","whomo","whomo2","len","ngap","ngap2"]

files={}
if args.typ:
    for metric in metrics:
        files[metric]=args.folder+"/"+metric+"."+args.typ+".csv"
else:
    for metric in metrics:
        files[metric]=args.folder+"/"+metric+".csv"



# Read nseq
fam_nseq={}
with open(args.nseqfile) as f:
    for line in f:
        line=line.strip("\n")
        fields=line.split("\t")
        fam_nseq[fields[0]]=fields[1]

# Read data frame
dfs={}
for name,fil in files.items():
    dfs[name]=pd.read_csv(fil,header=[0,1],index_col=0)
    dfs[name]=dfs[name].fillna("NA")

# Get headers
families=list(dfs["sp"].index)
aligners_rep=[x[0] for x in dfs["sp"].columns][1:]
trees_rep=[x[1] for x in dfs["sp"].columns][1:]
aligners,trees=[],[]
for x in aligners_rep:
    if x not in aligners:
        aligners.append(x)
for x in trees_rep:
    if x not in trees:
        trees.append(x)

# # Prepare data frame
# o.write("family aligner tree1 tree2 tc sp homo whomo whomo2 ngap ngap2\n")
# for fam in families:
#     for aligner in aligners:
#         for tree in trees:
#             tree1=tree.split(".")[0]
#             tree2=tree.split(".")[1]
#             o.write("{} {} {} {}".format(fam,aligner,tree1,tree2))
#             for metric in metrics:
#                 df=dfs[metric]
#                 value=df[aligner].loc[fam,tree]
#                 o.write(" "+str(value))
#             o.write("\n")
            

# Prepare data frame
df={}
df["family"]=[item for item in families for i in range(len(trees_rep))]
df["nseq"]=[fam_nseq[item] for item in families for i in range(len(trees_rep))]
df["mode"]=["reg_align"]*len(families)*len(trees_rep)
df["bucket"]=[str(args.bucket)]*len(families)*len(trees_rep)
df["aligner"]=aligners_rep*len(families)
tree_col=trees_rep*len(families)
df["tree1"]=[x.split(".")[0] for x in tree_col]
df["tree2"]=[x.split(".")[1] for x in tree_col]
for metric in metrics:
    df[metric]=[]
    for fam in families:
        df[metric]+=list(dfs[metric].loc[fam,:])[1:]
l=[]
for x in ["family","nseq","mode","bucket","aligner","tree1","tree2"]+metrics:
    l.append(df[x])
lt=map(list, zip(*l))

# Print
print(",".join(["family","nseq","mode","bucket","aligner","tree1","tree2"]+metrics))
for i in lt:
    tmp=",".join([str(x) for x in i])
    print(tmp)


