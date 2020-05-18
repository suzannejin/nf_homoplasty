#!/usr/bin/env nextflow


/*
MRdelta all vs all

score:  [tc, sp, col, sspa]
metric: [homo, whomo, whomo2, len, ngap, ngap2, seqid]  - negative
metric: [tc, sp, col, sspa]  - positive


score:  [homo, whomo, whomo2, len, ngap, ngap2, seqid]
metric: [tc, sp, col, sspa]  - negative
metric: [homo, whomo, whomo2, len, ngap, ngap2, seqid]  - positive
*/


//params.score = ['tc', 'sp', 'col','sspa1000']
//params.metric = ['tc', 'sp', 'col','sspa1000']
params.metric= ['homo', 'whomo', 'whomo2', 'len', 'ngap', 'ngap2','dTSseqidscore1000']
params.score = ['homo', 'whomo', 'whomo2', 'len', 'ngap', 'ngap2','dTSseqidscore1000']


//params.metric = ["sspa30","sspa40","sspa50","sspa60","sspa70","sspa80","sspa90","sspa100","sspa150","sspa200","sspa250","sspa300","sspa500","sspa1000"]
//params.metric = ['homo', 'whomo', 'whomo2', 'len', 'ngap', 'ngap2']
//params.metric = ['dTSseqidscore10','dTSseqidscore100','dTSseqidscore200','dTSseqidscore500','dTSseqidscore1000','dTMseqid10','dTMseqidscore10','dTMseqid100','dTMseqidscore100']
//params.metric = ['sTMseqid10','sTMseqid100']
//params.metric = ["sspa30","sspa40","sspa50","sspa60","sspa70","sspa80","sspa90","sspa100","sspa150","sspa200","sspa250","sspa300","sspa500","sspa1000", "rdsspa30","rdsspa40","rdsspa50","rdsspa60","rdsspa70","rdsspa80","rdsspa90","rdsspa100","rdsspa150","rdsspa200","rdsspa250","rdsspa300","rdsspa500","rdsspa1000","seqid10","seqid50","seqid100"]
//params.metric = ['homo','whomo','len','ngap2']
params.mrdelta = [0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0]
params.maintain = "bucket aligner family"
params.deltaby = "tree"
params.norm = " "
params.corr = "-corr positive" 


//params.csv = "/users/cn/sjin/projects/homoplasy/seqid/tc_quantest2_seqid.csv"
//params.output = "/users/cn/sjin/projects/homoplasy/nf_homoplasty/data/nf_mrdelta_all/deltaby_${params.deltaby}_positive"
//params.csv = "/users/cn/sjin/projects/homoplasy/results_quantest2/quantest2_with_informativeAln_FAMSA_reformat.csv"
//params.output = "/users/cn/sjin/projects/homoplasy/nf_homoplasty/data/nf_mrdelta_check_famsa_vs_mafft_informativeAln/deltaby_${params.deltaby}"
params.csv = "/users/cn/sjin/projects/homoplasy/final_results/raw/medium/tc_homo_quantest2_seqid.csv"
params.output = "/users/cn/sjin/projects/homoplasy/nf_homoplasty/data/MRdelta/medium/deltaby_${params.deltaby}"


log.info """\
         Running MRdelta - python version"
         ======================================="
         Input csv                                          : ${params.csv}
         Evaluation:
            Scores                                          : ${params.score}
            Metrics                                         : ${params.metric}
            Mrdelta                                         : ${params.mrdelta}
            Between                                         : ${params.deltaby}
         Output directory                                   : ${params.output}

         """
         .stripIndent()




process mrdelta {
    
    cache false

    tag "${score}_${metric}_mrdelta${mrdelta}"
    publishDir "${params.output}/mrdelta${mrdelta}", mode: 'copy', overwrite: true

    input:
        path(csv) from params.csv
        each score from params.score
        each metric from params.metric
        each mrdelta from params.mrdelta
        val(maintain) from params.maintain
        val(deltaby) from params.deltaby
        val(norm) from params.norm
        val(corr) from params.corr
    
    output:
        file("${score}_${metric}_mrdelta${mrdelta}*.tsv") into mrdeltaOut

    when:
        if(score == metric){false}else{true}

    script:
        """
        python /users/cn/sjin/projects/homoplasy/nf_homoplasty/bin/cedric/csv2analyse.py -csv ${csv} \
            -score ${score} -metric ${metric} -mrdelta ${mrdelta} -outdir . -maintain ${maintain} -deltaby ${deltaby} ${norm} ${corr}

        """
}