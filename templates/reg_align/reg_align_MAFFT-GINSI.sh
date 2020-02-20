# export DUMP_ALN_BUCKETS=1
export NO_MAFFT_BINARIES=1

t_coffee -reg -reg_method mafftginsi_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -n_core=1 \
         -outfile ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.aln