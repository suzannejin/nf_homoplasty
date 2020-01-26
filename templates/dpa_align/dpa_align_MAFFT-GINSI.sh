# export DUMP_ALN_BUCKETS=1
export NO_MAFFT_BINARIES=1
t_coffee -reg -reg_method mafftginsi_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -n_core=1 \
         -outfile ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.aln

if test -f "alndump..1.aln_bucket"; then
  mv alndump..1.aln_bucket ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.parent.aln
fi
## homo
awk 'NR == 1 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.homo
## weight homo
awk 'NR == 2 {print \$2}' ${id}.homoplasy > ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.w_homo

