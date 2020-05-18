bucket=$1
#declare -a mrdeltas=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1)
declare -a mrdeltas=(0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15)
#declare -a scores=(tc sp quantest30 quantest50 quantest150 quantest1000)
#declare -a metrics=(homo whomo whomo2 len ngap ngap2)
declare -a scores=(tc)
declare -a metrics=(quantest30 quantest50 quantest150 quantest1000)
reverse="yes"
declare -a typs=(original) # normPerLen normByLen)
#declare -a metrics=(whomo whomo2 whomoPerLen whomo2PerLen whomoByLen whomo2ByLen whomoPerSeq whomo2PerSeq whomoBySeq whomo2BySeq whomoPerLenSeq whomo2PerLenSeq whomoByLenSeq whomo2ByLenSeq)
declare -a aligners=(CLUSTALO MAFFT-FFTNS1) # MAFFT-GINSI MAFFT-SPARSECORE)
#declare -a familyName=(ltn     il8     az      kringle cryst   DEATH   cah     mmp     rub     ghf10   tgfb    sodcu   KAS     DMRL_synthase   tms     GEL     kunitz  Sulfotransfer   mofe    Ald_Xan_dh_2    ghf5    phc     aadh    annexin serpin  cytb    asp     oxidored_q6     hpr     hormone_rec     hr      tim     glob    ace     cys     ghf1    sodfe   peroxidase      uce     flav    HMG_box OTCace  msb     icd     proteasome      cyclo   LIM     HLH     ldh     subt    int     lyase_1 gpdh    egf     blm     gluts   myb_DNA-binding tRNA-synt_2b    biotin_lipoyl   hom     ghf13   aldosered       hla     Rhodanese       PDZ     blmb    rhv     p450    adh     aat     rrm     Acetyltransf    sdr     zf-CCHH rvp)
#familyList=$(cat /users/cn/sjin/projects/homoplasy/homfam_info/combinedSeqs)
familyList=/users/cn/sjin/projects/homoplasy/homfam_info/bigger1000_small2big

csv2analyse=bin/cedric/csv2analyse.pl
cedrictxt2colcsv=bin/cedric/cedrictxt2colcsv.py
mergeAlnGlobal=bin/cedric/merge_aln_global.R


# Create files
compute_homoplasy () {
    #ori=data/raw/bucket${bucket}  # Folder where the raw data files are stored
    #outori=data/mrdelta
    ori=data/raw_tc_quantest_over1000/bucket${bucket}
    #outori=data/mrdelta_tc_quantest_over1000
    outori=data/mrdelta_plot_deltatc_deltaquantest
    for mrdelta in ${mrdeltas[@]}; do
        folder=${outori}/bucket${bucket}/mrdelta${mrdelta}
        [[ ! -d $folder ]] && mkdir -p $folder
        for typ in ${typs[@]}; do
        for score in ${scores[@]}; do
            for metric in ${metrics[@]}; do
                echo "perl $csv2analyse -mrdelta $mrdelta -score $score -metrics $metric -norm $typ -dir $ori -reverse $reverse > ${folder}/${score}_${metric}_${typ}.txt"
                perl $csv2analyse -mrdelta $mrdelta -score $score -metrics $metric -norm $typ -dir $ori -reverse $reverse > ${folder}/${score}_${metric}_${typ}.txt
                python $cedrictxt2colcsv ${folder}/${score}_${metric}_${typ}.txt
            done
        done
        done
    done
}

# Merge data 
# This will be used to plot
merge_data () {
    echo "+mrdelta merge $bucket"
    #ori=data/mrdelta/bucket${bucket}
    #outori=data/mrdelta  # Where mrdelta merged output tables are stored
    #ori=data/mrdelta_tc_quantest_over1000/bucket${bucket}
    #outori=data/mrdelta_tc_quantest_over1000
    ori=data/mrdelta_plot_deltatc_deltaquantest/bucket${bucket}
    outori=data/mrdelta_plot_deltatc_deltaquantest
    outfam=${outori}/out.bucket${bucket}.fam.tsv
    outaln=${outori}/out.bucket${bucket}.aln.tsv
    outglobal=${outori}/out.bucket${bucket}.global.tsv
    printf "mrdelta\tscore\tmetric\ttyp\tbucket\tfamily\taligner\tNratio\tNnumber\tusedNumber\tunusedRatio\tunusedNumber\ttotal\n" > $outfam
    printf "mrdelta\tscore\tmetric\ttyp\tbucket\taligner\tNratio\tunusedRatio\tnfam\tminacc\tmaxacc\tdeltaacc\n" > $outaln
    printf "mrdelta\tscore\tmetric\ttyp\tbucket\taligner\tNratio\tunusedRatio\tminacc\tmaxacc\tdeltaacc\n" > $outglobal
    echo -e "\t+mrdelta merge global $bucket"
    for mrdelta in ${mrdeltas[@]}; do
    for typ in ${typs[@]}; do
        for score in ${scores[@]}; do
            for metric in ${metrics[@]}; do
                awk -v mrdelta=$mrdelta -v score=$score -v metric=$metric -v typ=$typ -v bucket=$bucket 'NR!=1{print mrdelta"\t"score"\t"metric"\t"typ"\t"bucket"\tglobal\t"$0}' ${ori}/mrdelta${mrdelta}/${score}_${metric}_${typ}_global.tsv >> $outglobal 
            done
        done
    done
    done
    echo -e "\t+mrdelta merge aln $bucket"
    for mrdelta in ${mrdeltas[@]}; do
    for typ in ${typs[@]}; do
        for score in ${scores[@]}; do
            for aligner in ${aligners[@]}; do
                for metric in ${metrics[@]}; do
                awk -v mrdelta=$mrdelta -v score=$score -v metric=$metric -v typ=$typ -v aligner=$aligner -v bucket=$bucket 'NR!=1&&$1==aligner{print mrdelta"\t"score"\t"metric"\t"typ"\t"bucket"\t"$0}' ${ori}/mrdelta${mrdelta}/${score}_${metric}_${typ}_aln.tsv >> $outaln
                done
            done
        done
    done
    done
    echo -e "\t+mrdelta merge fam $bucket"
    for mrdelta in ${mrdeltas[@]}; do
    for typ in ${typs[@]}; do
        for score in ${scores[@]}; do
            for metric in ${metrics[@]}; do
                awk -v mrdelta=$mrdelta -v score=$score -v metric=$metric -v typ=$typ -v bucket=$bucket 'NR!=1{print mrdelta"\t"score"\t"metric"\t"typ"\t"bucket"\t"$0}' ${ori}/mrdelta${mrdelta}/${score}_${metric}_${typ}_fam.tsv >> $outfam
            done
        done
    done
    done

}


compute_homoplasy
merge_data

