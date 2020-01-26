


declare -a familyName=(ltn     il8     az      kringle cryst   DEATH   cah     mmp     rub     ghf10   tgfb    sodcu   KAS     DMRL_synthase   tms     GEL     kunitz  Sulfotransfer   mofe    Ald_Xan_dh_2    ghf5    phc     aadh    annexin serpin  cytb    asp     oxidored_q6     hpr     hormone_rec     hr      tim     glob    ace     cys     ghf1    sodfe   peroxidase      uce     flav    HMG_box OTCace  msb     icd     proteasome      cyclo   LIM     HLH     ldh     subt    int     lyase_1 gpdh    egf     blm     gluts   myb_DNA-binding tRNA-synt_2b    biotin_lipoyl   hom     ghf13   aldosered       hla     Rhodanese       PDZ     blmb    rhv     p450    adh     aat     rrm     Acetyltransf    sdr     zf-CCHH rvp)

declare -a aligner=(CLUSTALO)

declare -a tree=(codnd dpparttreednd1 dpparttreednd2 dpparttreednd2size fastaparttreednd fftns1dnd fftns1dndmem fftns2dnd fftns2dndmem mafftdnd parttreednd0 parttreednd1 parttreednd2 parttreednd2size)


printf "\t\t########################\n"
printf "\t\t######### TC ###########\n"
printf "\t\t########################\n"

printf "Family \t"

for z in ${tree[@]} 
do 
  printf ${z}"\t"
done
printf "\n"

for i in ${familyName[@]}
do
  printf "${i}\t"
	
	for x in ${aligner[@]}
    do
        for z in ${tree[@]}
        do
			    cat ${i}.dpa_align.1000.${x}.${z}.tc
       	done
    done
    printf "\n"
done

printf "\t\t########################\n"
printf "\t\t######### W_HOMO ###########\n"
printf "\t\t########################\n"

printf "Family \t"

for z in ${tree[@]} 
do 
  printf ${z}"\t"
done
printf "\n"

for i in ${familyName[@]}
do
  printf "${i}\t"
  
  for x in ${aligner[@]}
    do
        for z in ${tree[@]}
        do
          cat ${i}.dpa_1000.${x}.with.${z}.tree.w_homo 
        done
    done
    printf "\n"
done


