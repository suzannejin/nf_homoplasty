t_coffee -reg -reg_method clustalo_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.aln

## homo
awk 'NR == 1 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.homo
## weight homo
awk 'NR == 2 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo
