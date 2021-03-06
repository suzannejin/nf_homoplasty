bucket=$1
outputfolder=$2

declare -a aligner=(CLUSTALO MAFFT-FFTNS1) # MAFFT-GINSI MAFFT-SPARSECORE)
declare -a foldname=(CO FFTNS1) # GINSI SPARSE)

##declare -a familyName=(seatoxin hip scorptoxin cyt3 rnasemam bowman toxin ghf11 TNF sti Stap_Strp_toxin profilin ricin ghf22 ChtBD ins trfl slectin phoslip ltn il8 az kringle cryst DEATH cah mmp rub ghf10 tgfb sodcu KAS DMRL_synthase tms GEL kunitz Sulfotransfer mofe Ald_Xan_dh_2 ghf5 phc aadh annexin serpin cytb asp oxidored_q6 hpr hormone_rec hr tim glob ace cys ghf1 sodfe peroxidase uce flav HMG_box OTCace msb icd proteasome cyclo LIM HLH ldh subt int lyase_1 gpdh egf blm gluts myb_DNA-binding tRNA-synt_2b biotin_lipoyl hom ghf13 aldosered hla Rhodanese PDZ blmb rhv p450 adh aat rrm Acetyltransf sdr zf-CCHH rvp)
##declare -a familyName=(seatoxin)
#familyList=/users/cn/sjin/projects/homoplasy/homfam_info/combinedSeqs_nseqs
familyList=/users/cn/sjin/projects/homoplasy/homfam_info/bigger1000_small2big

#declare -a tree=(codnd dpparttreednd1 dpparttreednd2 dpparttreednd2size fastaparttreednd fftns1dnd fftns1dndmem fftns2dnd fftns2dndmem mafftdnd parttreednd0 parttreednd1 parttreednd2 parttreednd2size FAMSA CLUSTALO-RANDOM)
declare -a tree=(codnd FAMSA mafftdnd dpparttreednd1 fastaparttreednd fftns1dnd parttreednd0 parttreednd2 CLUSTALO-RANDOM)

declare -a tcs=(tc sp col)
declare -a homos=(homo w_homo w_homo2 len ngap ngap2)

declare -a n_quantest2=(30 50 100 150 1000)



[[ ! -d $outputfolder ]] && mkdir -p $outputfolder

print_tc () {

  typ=$1  # tc / sp / col
  TYP=${typ^^}

  # Print aligners header
  printf "Family,nseq"
  for a in ${aligner[@]}
  do
    for z in ${tree[@]} 
    do 
        printf ","${a}
    done
  done
  printf "\n"

  # Print trees header
  printf "Family,nseq"
  for a in ${aligner[@]}
  do
    for z in ${tree[@]} 
    do 
        printf ","${z}
    done
  done
  printf "\n"

  # Print values
  ##for i in ${familyName[@]}
  for i in $(cut -f1 $familyList)
  do
    printf "${i},"
    awk -v i=$i '{if($1==i){printf $2}}' $familyList
          n=0
          for x in ${aligner[@]}
      do
          f=results_fullTree_${foldname[$n]}_bucket${bucket}
          for z in ${tree[@]}
          do
            fil=${f}/individual_scores/${i}.dpa_align.${bucket}.${x}.${z}.${typ}
            if [[ -s $fil ]]; then
              awk 'NR==1{printf ","$1}' ${fil}
            else
              printf ",NA"
            fi
          done
          n=$(( $n+1 ))
      done
      printf "\n"
  done

}

print_homo () {

typ=$1  # homo / whomo / whomo2 / len / ngap / ngap2
TYP=${typ^^}

# Print aligners header
printf "Family,nseq"
for a in ${aligner[@]}
do
for z in ${tree[@]} 
do 
    printf ","${a}
done
done
printf "\n"

# Print trees header
printf "Family,nseq"
for a in ${aligner[@]}
do
for z in ${tree[@]} 
do 
    printf ","${z}
done
done
printf "\n"

# Print values
##for i in ${familyName[@]}
for i in $(cut -f1 $familyList)
do
  printf "${i},"
  awk -v i=$i '{if($1==i){printf $2}}' $familyList
  
  n=0
  for x in ${aligner[@]}
  do
      f=results_fullTree_${foldname[$n]}_bucket${bucket}
      for z in ${tree[@]}
      do
        fil=${f}/alignments/${i}.dpa_${bucket}.${x}.with.${z}.tree.${typ}
        
        if [[ -s $fil ]]; then
          if [[ $typ == "ngap2" ]]; then 
            awk 'NR==1{printf ","$1}' ${fil} | tr -d "l"
          else
            awk 'NR==1{printf ","$1}' ${fil}
          fi
        else
          printf ",NA"
        fi
      done
      n=$(( $n+1 ))
  done
  printf "\n"
done


}

print_quantest () {

nquan=$1

# Print aligners header
printf "Family,nseq"
for a in ${aligner[@]}
do
for z in ${tree[@]} 
do 
    printf ","${a}
done
done
printf "\n"

# Print trees header
printf "Family,nseq"
for a in ${aligner[@]}
do
for z in ${tree[@]} 
do 
    printf ","${z}
done
done
printf "\n"

# Print values
##for i in ${familyName[@]}
for i in $(cut -f1 $familyList)
do
  printf "${i},"
  awk -v i=$i '{if($1==i){printf $2}}' $familyList
  
  n=0
  for x in ${aligner[@]}
  do
      for z in ${tree[@]}
      do
        fil=/users/cn/sjin/projects/homoplasy/quantest2/out/n${nquan}/${i}/${i}.dpa_${bucket}.${x}.with.${z}.tree.random${nquan}.quantest2
        if [[ -s $fil ]]; then
            awk 'NR==1{printf ","$1}' ${fil}
        else
          printf ",NA"
        fi
      done
      n=$(( $n+1 ))
  done
  printf "\n"
done

}


for n in ${n_quantest2[@]}
do
  printf "\t+print quantest${n} bucket${bucket}\n"
  print_quantest $n > ${outputfolder}/quantest${n}.csv
done

for i in ${tcs[@]}
do
  printf "\t+print raw $i bucket${bucket}\n"
  print_tc $i > ${outputfolder}/${i}.csv
done


for i in ${homos[@]}
do
  printf "\t+print raw $i bucket${bucket}\n"
  print_homo $i > ${outputfolder}/${i}.csv
done
mv ${outputfolder}/w_homo.csv ${outputfolder}/whomo.csv
mv ${outputfolder}/w_homo2.csv ${outputfolder}/whomo2.csv