t_coffee -dpa -dpa_method t_coffee_msa \
         -dpa_tree ${guide_tree} \
         -seq ${seqs} \
         -dpa_nseq ${bucket_size} \
         -n_core=1 \
         -outfile ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.aln
