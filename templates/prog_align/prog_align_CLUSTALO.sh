clustalo --infile=${seqs} \
         --guidetree-in=${guide_tree} \
         --outfmt=fa \
         -o ${id}.prog_align.NA.${align_method}.with.${tree_method}.tree.aln
