
declare -a flavour=(reg_align)
declare -a aligner=(CLUSTALO MAFFT-FFTNS1 MAFFT-GINSI)
declare -a folders=(results_fullTree_CO results_fullTree_FFTNS1 results_fullTree_GINSI)

#declare -a familyName=(seatoxin hip scorptoxin cyt3 rnasemam bowman toxin ghf11 TNF sti Stap_Strp_toxin profilin ricin ghf22 ChtBD ins trfl slectin phoslip ltn il8 az kringle cryst DEATH cah mmp rub ghf10 tgfb sodcu KAS DMRL_synthase tms GEL kunitz Sulfotransfer mofe Ald_Xan_dh_2 ghf5 phc aadh annexin serpin cytb asp oxidored_q6 hpr hormone_rec hr tim glob ace cys ghf1 sodfe peroxidase uce flav HMG_box OTCace msb icd proteasome cyclo LIM HLH ldh subt int lyase_1 gpdh egf blm gluts myb_DNA-binding tRNA-synt_2b biotin_lipoyl hom ghf13 aldosered hla Rhodanese PDZ blmb rhv p450 adh aat rrm Acetyltransf sdr zf-CCHH rvp)
declare -a familyName=(seatoxin)
familyList=/users/cn/sjin/projects/homoplasy/homfam_info/combinedSeqs_nseqs

declare -a tree=(codnd dpparttreednd1 dpparttreednd2 dpparttreednd2size fastaparttreednd fftns1dnd fftns1dndmem fftns2dnd fftns2dndmem mafftdnd parttreednd0 parttreednd1 parttreednd2 parttreednd2size FAMSA CLUSTALO-RANDOM)
##declare -a tree=(codnd dpparttreednd1 dpparttreednd2)

declare -a buckets=(50 100 200 500 1000)

declare -a tcs=(tc sp col)
declare -a homos=(homo w_homo w_homo2 len ngap ngap2)


echo "family;mode;bucket;aligner;tree;nseq;tc;sp;col;homo;whomo;whomo2;len;ngap;ngap2"
#printf "family aligner tree bucket nseq tc sp col homo whomo whomo2 len ngap ngap2\n"

for fam in ${familyName[@]}; do
    nseq=$(grep $fam $familyList | awk '{print $2}')
    for mode in ${flavour[@]}; do
    i=0
    for bucket in ${buckets[@]}; do
    for aln in ${aligner[@]}; do
        fold=${folders[$i]}
        for tr in ${tree[@]}; do
            printf "$fam;$mode;$bucket;$aln;$tr;$nseq;"
            for tc in ${tcs[@]}; do
                fil=${fold}_bucket${bucket}/individual_scores/${fam}.dpa_align.${bucket}.${aln}.${tr}.${tc}
                if [[ -s $fil ]]; then
                    awk 'NR==1{printf $1";"}' ${fil}
                else
                    printf "NA;"
                fi
            done
            for hom in ${homos[@]}; do
                fil=${fold}_bucket${bucket}/alignments/${fam}.dpa_${bucket}.${aln}.with.${tr}.tree.${hom}
                if [[ -s $fil ]]; then
                    if [[ $hom == "ngap2" ]]; then 
                        awk 'NR==1{printf $1"\n"}' ${fil} | tr -d "l"
                    else
                        awk 'NR==1{printf $1";"}' ${fil}
                    fi
                else
                    if [[ $hom == "ngap2" ]]; then 
                        printf "NA\n"
                    else
                        printf "NA;"
                    fi
                fi
            done
        done
        i=$((i+1))
    done
    done
    done
done

