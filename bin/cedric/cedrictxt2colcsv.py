import collections
import operator
import os.path
import re
import sys

fil=sys.argv[1]


class Fam:
    def __init__(self,dic,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint):
        self.dic=dic  # FAM, BUCKET, ALIGNER
        self.Nratio=Nratio
        self.Nnumber=Nnumber
        self.usedNumber=usedNumber
        self.unusedRatio=unusedRatio
        self.unusedNumber=unusedNumber
        self.totalpoint=totalpoint

class Global:
    def __init__(self,dic,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint,nfam,minacc,maxacc,deltaacc):
        self.dic=dic  # BUCKET | ALIGNER | TREE
        self.Nratio=Nratio
        self.Nnumber=Nnumber
        self.usedNumber=usedNumber
        self.unusedRatio=unusedRatio
        self.unusedNumber=unusedNumber
        self.totalpoint=totalpoint
        self.nfam=nfam
        self.minacc=minacc
        self.maxacc=maxacc
        self.deltaacc=deltaacc


fams=[]
glob={}
with open(fil) as f:
    for line in f:
        line=line.strip("\n")
        fields=line.split(" ")

        if fields[0]=="----Process":
            continue

        # Get N / total used ratio
        Nratio=float(fields[1])

        # Get number of N
        tmp=fields[2][1:-1].split("/")
        Nnumber=int(tmp[0])
        usedNumber=int(tmp[1])

        # Get unused ratio
        unusedRatio=float(fields[5])
        tmp=fields[6][1:-1].split("/")
        unusedNumber=int(tmp[0])
        totalpoint=int(tmp[1])


        if line[0:3]!="ALL":

            # Get typ. Eg. [FAM:BUCKET:ALIGNER]:[seatoxin:50:CLUSTALO]
            tmp=line.split("::")[0]
            print(line)
            match=re.findall("(\[[A-Za-z:0-9_-]+\])+",tmp) 
            print(match)
            typ=match[0][1:-1].split(":")
            val=match[1][1:-1].split(":")
            dic=collections.OrderedDict()
            for n in range(len(typ)):
                t=typ[n]
                v=val[n]
                dic[t]=v

            if "FAM" in dic.keys():
                tmp=Fam(dic,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint)
                fams.append(tmp)
            else:
                # Get nfam
                nfam=fields[7][1:]

                # Get acc
                minacc=float(fields[10])
                maxacc=float(fields[13])
                deltaacc=float(fields[16][:-1])

                tmp=Global(dic,Nratio,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint,nfam,minacc,maxacc,deltaacc)
                if match[0] not in glob:
                    glob[match[0]]=[]
                glob[match[0]].append(tmp)

        else:
            # Get nfam
            nfam=fields[7][1:]

            # Get acc
            minacc=float(fields[10])
            maxacc=float(fields[13])
            deltaacc=float(fields[16][:-1])

            tmp=Global("ALL",Nratio,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint,nfam,minacc,maxacc,deltaacc)
            if match[0] not in glob:
                glob["ALL"]=[]
            glob["ALL"].append(tmp)


# Write output
pre=os.path.splitext(fil)[0]
o=open(pre+"_fam.tsv","w")
o.write("family\taligner\tNratio\tNnumber\tusedNumber\tunusedRatio\tunusedNumber\ttotal\n")
for c in sorted(fams, key=lambda x: (x.aligner, x.fam)):
    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(c.fam, c.aligner, c.Nratio, c.Nnumber, c.usedNumber, c.unusedRatio, c.unusedNumber, c.totalpoint))
o.close()
o=open(pre+"_aln.tsv","w")
o.write("aligner\tNratio\tunusedRatio\tnfam\tminacc\tmaxacc\tdeltaacc\n")
for c in alns:
    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(c.aligner, c.Nratio, c.unusedRatio, c.nfam, c.minacc, c.maxacc, c.deltaacc))
o.close()
o=open(pre+"_global.tsv","w")
o.write("Nratio\tunusedRatio\tminacc\tmaxacc\tdeltaacc\n")
o.write("{}\t{}\t{}\t{}\t{}\n".format(glob.Nratio, glob.unusedRatio, glob.minacc, glob.maxacc, glob.deltaacc))
o.close()