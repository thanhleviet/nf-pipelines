#!/usr/bin/env nextflow

/*
================================================================================
=                     BACTERIAL | A S S E M B L Y | P I P E L I N E            =
================================================================================
@Author
Thanh Le Viet <lethanhx2k@gmail.com>
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/thanhleviet/bacteria-assembly
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/thanhleviet/bacteria-assembly/README.md
--------------------------------------------------------------------------------
 @Licence
 https://github.com/thanhleviet/bacteria-assembly/LICENSE
--------------------------------------------------------------------------------
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

params.input = "/home/ubuntu/data/seq_ref_lab_2018/May30"
params.output = "results_may_30"
params.pattern = "*_L001_R{1,2}_001.fastq.gz"
params.out_dir = "/home/ubuntu/data/seq_ref_lab_2018/results/" + params.output
params.arbicate_db = ["resfinder", "card", "vfdb", "plasmidfinder", "ncbi"]

params.threads = 8
params.ram = 15

read_pairs = params.input + "/" + params.pattern
threads = params.threads

Channel
    .fromFilePairs(read_pairs, flat: true)
    .ifEmpty { exit 1, "Read pairs could not be found: ${read_pairs}" }
    .set { read_pairs }


process shovill {
    publishDir "${params.out_dir}/shovill", mode: "copy", pattern: "${dataset_id}"

    tag { dataset_id }

    input:
    set dataset_id, file(forward), file(reverse) from read_pairs

    output:
    file("${dataset_id}")
    set dataset_id, file("${dataset_id}_shovill_contigs.fa") into (abricate_ch, mlst_ch, sendsketch)

    """
    shovill --outdir ${dataset_id} \
    --R1 $forward --R2 $reverse \
    --cpus $threads \
    --trim

    cp ${dataset_id}/contigs.fa ${dataset_id}_shovill_contigs.fa
    mv ${dataset_id}/contigs.fa ${dataset_id}/${dataset_id}_shovill_contigs.fa
    """
    }

process abricating {
    publishDir "${params.out_dir}/annotation", mode: "copy"

    tag { "${dataset_id} ${db}"}

    input:
    set dataset_id, file("${dataset_id}_shovill_contigs.fa") from abricate_ch
    each db from params.arbicate_db
    output:
    file("*.tsv") into abricate_output

    script:

    """
    abricate --db ${db} ${dataset_id}_shovill_contigs.fa > ${dataset_id}_abricate_${db}.tsv
    """

}


process mlst {
    publishDir "${params.out_dir}/annotation", mode: "copy"

    tag { "${dataset_id}"}

    input:
    set dataset_id, file("${dataset_id}_shovill_contigs.fa") from mlst_ch
    output:
    file("*.tsv")

    script:

    """
    mlst ${dataset_id}_shovill_contigs.fa > ${dataset_id}_mlst.tsv
    """

}


process sendsketch {
    publishDir "${params.out_dir}/species_id", mode: "copy"

    tag {dataset_id}

    input:
    set dataset_id, file(contig) from sendsketch

    output:
    file "${dataset_id}.tsv" into report

    """
    sendsketch.sh in=${contig} color=f minbases=1000 printall > ${dataset_id}.tsv
    """
}
 // address=http://192.168.23.115:3071/sketch

process summary {
    publishDir "${params.out_dir}/species_id", mode: "copy"

    input:
    file tsv from report.collect()

    output:
    file("summary.tsv") into summary_python

    shell:
    '''
    echo "SampleID\tWKID\tKID\tANI\tComplt\tContam\tContam2\tuContam\tScore\tDepth\tDepth2\tVolume\tMatches\tUnique\tUnique2\tUnique3\tnoHit\tLength\tTaxID\tImgID\tgBases\tgKmers\tgSize\tgSeqs\trDiv\tqDiv\trSize\tqSize\tcHits\ttaxName\tseqName\ttaxonomy" > summary.tsv
    awk 'FNR==4 {print FILENAME, "\t",$0}' !{tsv} | sed 's/_shovill_contigs.tsv//g' >> summary.tsv
    '''

}

process filter_summary_python {
    publishDir "${params.out_dir}/species_id", mode:"copy"

    input:
    file("summary.tsv") from summary_python

    output:
    file("filter_summary_python.tsv")

    script:

    """
    #!/usr/bin/env python
    with open(\"filter_summary_python.tsv\",\"w\") as fo:
        with open(\"summary.tsv\",\"r\") as fi:
            next(fi)
            for line in fi:
                fields = line.split("\t")
                fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\\n".format(fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[28], fields[29]))
    """
}

workflow.onComplete {
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command Line:     $workflow.commandLine"
    log.info "Container:        $workflow.container"
    log.info "Duration:     $workflow.duration"
    log.info "Output Directory: $params.out_dir"


    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
    //Need to config mail (see nextflow.config)
    sendMail(to: 'lethanhx2k@gmail.com', subject: 'Assembly finished!', body: msg)
}


