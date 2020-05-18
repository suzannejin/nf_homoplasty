folder=$1
typ=$2  # fam, aln, global
out=${folder}/out.${typ}.tsv

cat ${folder}/*.${typ}.tsv | head -n 1 > $out
for fil in $(echo ${folder}/*.${typ}.tsv); do
    awk 'NR!=1{print}' $fil >> $out
done