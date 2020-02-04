t_coffee -reg -reg_method clustalo_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.aln

echo "t_coffee -reg -reg_method clustalo_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.aln" > ${id}.test

## Homoplasy
cat ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.homoplasy
## homo
awk 'NR == 1 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.homo
## weight homo
awk 'NR == 2 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo
## weight homo 2
awk 'NR == 3 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo2
## alignment length
awk 'NR == 4 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.len
## number gaps
awk 'NR == 5 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.ngap
## square number gaps
awk 'NR == 6 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.ngap2

