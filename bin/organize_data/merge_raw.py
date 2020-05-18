##
## Merge the files from different bucket sizes to one, in order to be used in csv2analyse.pl 
##

import collections
import sys
import argparse

app = argparse.ArgumentParser()
app.add_argument('--bucket',type=int,nargs='+',help="Bucket size(s)")
app.add_argument('--dir',type=str,default="data/raw")
app.add_argument('--metric',type=str)
app.add_argument('--out',type=str)
args=app.parse_args()


# Read files
files=collections.OrderedDict()
for bucket in args.bucket:
    #fil=open(args.dir+"/"+args.metric+".csv")
    fil=open(args.dir+"/bucket"+str(bucket)+"/"+args.metric+".csv")
    files[bucket]=fil.readlines()
    fil.close()


out=open(args.out,"w")

# Add bucket size header
line1=files[args.bucket[0]][0].strip("\n").split(",")
ntimes=len(line1)-2  # N times bucket size  // minus Family & nseq columns
line="Family,nseq"
for bucket in args.bucket:
    tmp=(","+str(bucket))
    line+=tmp*ntimes
out.write(line+"\n")

# Write content
nrow=len(files[args.bucket[0]])
for i in range(nrow):
    # Add first bucket file content
    fil=files[args.bucket[0]]
    line=fil[i].replace("\n","")
    out.write(line)
    #sys.stdout.write(line)
    # Add rest bucket file content
    nbuckets=len(args.bucket)  
    for j in range(1,nbuckets):
        fil=files[args.bucket[j]]
        line=fil[i].replace("\n","").split(",")[2:]   # Exclude family,nseq columns
        line=",".join(line)
        out.write(","+line)
        #sys.stdout.write(","+line)
        #sys.stdout.write("\n")
        #sys.stdout.write("\n")
    #sys.stdout.write("\n")
    #if i >5:
        #exit(1)
    out.write("\n")
out.close()
        
