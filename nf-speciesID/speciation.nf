#!/usr/bin/env nextflow
params.reads = "*_L001_R{1,2}_001.fastq.gz"
params.input = "/home/ubuntu/data/seq_ref_lab_2018/May12"
params.read_pairs = params.input + "/" + params.reads
params.centrifuge_db = "/home/ubuntu/data/centrifuge/library/p_compressed"
params.krona = "/home/ubuntu/data/krona/taxonomy"
params.output = "/home/ubuntu/data/seq_ref_lab_2018/results"
params.threads = 8

Channel
    .fromFilePairs(params.read_pairs, flat: true)
    .ifEmpty{ exit 1, "Found no input reads, did you specify --read_pairs? I got: '${params.read_pairs}'"}
    .into {input_centrifuge; input_genome_size}

centrifuge_db = file(params.centrifuge_db)

threads = params.threads

process centrifuge {
    tag {sample_id}

    publishDir "${params.output}/speciation", mode: 'copy'

    errorStrategy 'terminate'

    input:
    set sample_id, file(forward), file(reverse) from input_centrifuge
    val centrifuge_db
    val params.krona

    output:
    file ("*")
    set sample_id, file("${sample_id}_speciation.csv"), file ("${sample_id}_genome_size.tsv") into output_speciation
    script:
    """
    genome_size.py --R $forward > ${sample_id}_genome_size.tsv

    centrifuge \
        --time \
        --threads ${threads} \
        -x ${centrifuge_db} \
        -1 ${forward} \
        -2 ${reverse} \
        -S ${sample_id}.centrifuge.tsv \
        --report-file ${sample_id}.centrifuge_report.tsv
    # Kraken like report
    centrifuge-kreport -x ${centrifuge_db} ${sample_id}.centrifuge.tsv > ${sample_id}_kraken_like.tsv

    getSpecies.py ${sample_id}.centrifuge_report.tsv "${params.output}/speciation/species_reads"

    # Visualization with Krona
    cat ${sample_id}.centrifuge.tsv | cut -f 1,3 > ${sample_id}.krona

    ktImportTaxonomy -tax ${params.krona} -o ${sample_id}.krona.html ${sample_id}.krona
    """
}


process summary {
    tag {sample_id}

    input:
    set sample_id, file("${sample_id}_speciation.csv"), file ("${sample_id}_genome_size.tsv") from output_speciation

    output:
    file("${sample_id}.csv") into output_summary

    script:
    """
    #!/usr/bin/env Rscript

    if (!"pacman" %in% row.names(installed.packages())) {
        chooseCRANmirror(ind=110)
        install.packages("pacman")
    }

    library(pacman)

    p_load(readr, dplyr)

    species <- read_csv("${sample_id}_speciation.csv", col_names = FALSE)

    genome_size <- read_csv("${sample_id}_genome_size.tsv", col_names = FALSE)

    all <- bind_cols(species, genome_size)
    write_csv(all, "${sample_id}.csv", col_names = FALSE)
    """

}

process cat_summary {
    publishDir "${params.output}/speciation", mode: 'copy'

    input:
    file sum from output_summary.collect()

    output:
    file "report_species.csv"

    """
    echo "sample_id,species,genome_size,coverage" > report_species.csv
    cat $sum >> report_species.csv
    """
}


workflow.onComplete {
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command Line:     $workflow.commandLine"
    log.info "Container:        $workflow.container"
    log.info "Duration:     $workflow.duration"
    log.info "Output Directory: $params.output"
}
