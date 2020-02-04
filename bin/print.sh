



declare -a familyName=(ltn     il8     az      kringle cryst   DEATH   cah     mmp     rub     ghf10   tgfb    sodcu   KAS     DMRL_synthase   tms     GEL     kunitz  Sulfotransfer   mofe    Ald_Xan_dh_2    ghf5    phc     aadh    annexin serpin  cytb    asp     oxidored_q6     hpr     hormone_rec     hr      tim     glob    ace     cys     ghf1    sodfe   peroxidase      uce     flav    HMG_box OTCace  msb     icd     proteasome      cyclo   LIM     HLH     ldh     subt    int     lyase_1 gpdh    egf     blm     gluts   myb_DNA-binding tRNA-synt_2b    biotin_lipoyl   hom     ghf13   aldosered       hla     Rhodanese       PDZ     blmb    rhv     p450    adh     aat     rrm     Acetyltransf    sdr     zf-CCHH rvp)
##declare -a familyName=(seatoxin)

declare -a aligner=(MAFFT-GINSI)  ## CLUSTALO, MAFFT-FFTNS1, MAFFT-GINSI, MAFFT-SPARSECORE, UPP

declare -a tree=(codnd dpparttreednd1 dpparttreednd2 dpparttreednd2size fastaparttreednd fftns1dnd fftns1dndmem fftns2dnd fftns2dndmem mafftdnd parttreednd0 parttreednd1 parttreednd2 parttreednd2size)
##declare -a tree=(codnd dpparttreednd1 dpparttreednd2)


declare -a tcs=(tc sp col)
declare -a homos=(homo w_homo w_homo2 len ngap ngap2)

f=results_fullTree_GINSI

print_tc () {

  typ=$1  # tc / sp / col
  TYP=${typ^^}

  printf "\t\t########################\n"
  printf "\t\t######### $TYP ###########\n"
  printf "\t\t########################\n"

  printf "Family "

  for z in ${tree[@]} 
  do 
    printf ${z}" "
  done
  printf "\n"

  for i in ${familyName[@]}
  do
    printf "${i} "
          
          for x in ${aligner[@]}
      do
          for z in ${tree[@]}
          do
            fil=${f}/individual_scores/${i}.dpa_align.1000.${x}.${z}.${typ}
            if [[ -s $fil ]]; then
              awk 'NR==1{printf $1" "}' ${fil}
            else
              printf "NA "
            fi
          done
      done
      printf "\n"
  done

}


print_homo () {

typ=$1  # homo / w_homo / w_homo2 / len / ngap / ngap2
TYP=${typ^^}

printf "\t\t########################\n"
printf "\t\t######### $TYP ###########\n"
printf "\t\t########################\n"

printf "Family "

for z in ${tree[@]} 
do 
  printf ${z}" "
done
printf "\n"

for i in ${familyName[@]}
do
  printf "${i} "
  
  for x in ${aligner[@]}
    do
        for z in ${tree[@]}
        do
          fil=${f}/alignments/${i}.dpa_1000.${x}.with.${z}.tree.${typ}
          if [[ -s $fil ]]; then
            awk 'NR==1{printf $1" "}' ${fil}
            [[ $typ == "ngap2" ]] && awk 'NR==1{printf $1" "}' ${fil} | tr -d "l"
          else
            printf "NA "
          fi
        done
    done
  printf "\n"
done


}


for i in ${tcs[@]}
do
  print_tc $i
  printf "\n"
done

for i in ${homos[@]}
do
  print_homo $i
  printf "\n"
done

