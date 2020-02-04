t_coffee -reg -reg_method mafft_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -outfile ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.aln


