#!/usr/bin/env nextflow

// General configuration variables
params.input = "/home/ubuntu/data/linh/pipelines/batch2"
params.read_pattern = "*_L001_R{1,2}_001.fastq.gz"
params.output = ""
params.help = false
params.read_pairs = params.input + "/" +params.read_pattern
params.threads = 4
params.ref = ""

params.snippy = false
params.phenix = true


threads = params.threads
ref = file(params.ref)
phenix_config = file("$baseDir/config/phenix.yaml")

Channel
    .fromFilePairs(params.read_pairs, flat: true)
    .ifEmpty { exit 1, "Read pairs could not be found: ${params.read_pairs}" }
    .into { read_pairs_snippy; read_pairs_phenix }

if (params.snippy) {
    process snippying {
    publishDir params.output, mode: "copy"

    tag { dataset_id }

    input:
    set dataset_id, file(forward), file(reverse) from read_pairs_snippy
    file ref

    output:
    file "${dataset_id}"

    """
    /home/ubuntu/data/linh/pipelines/snippy/bin/snippy --cpus $threads --prefix ${dataset_id} --outdir ${dataset_id} --ref $ref --pe1 $forward --pe2 $reverse

    qualimap bamqc -bam ${dataset_id}/${dataset_id}.bam -outdir stats
    mv stats/genome_results.txt ${dataset_id}/${dataset_id}_stats_cov.txt
    rm -r stats
    """
    }
}

if (params.phenix) {
    process phenixing {
    publishDir params.output, mode: "move"

    tag { dataset_id }

    input:
    set dataset_id, file(forward), file(reverse) from read_pairs_phenix
    file phenix_config
    file(ref)

    output:
    file "${dataset_id}"

    """
    phenix.py prepare_reference -r $ref \
    --mapper bwa \
    --variant gatk
    
    phenix.py run_snp_pipeline \
    -r1 $forward \
    -r2 $reverse \
    -r ${ref} \
    -c ${phenix_config} \
    --keep-temp \
    --json \
    --sample-name ${dataset_id} \
    -o ${dataset_id}

    phenix.py vcf2fasta \
    -i ${dataset_id}/${dataset_id}.filtered.vcf \
    -o ${dataset_id}/${dataset_id}_all.fasta \
    --reference ${params.ref} \
    --regex filtered
    
    seqkit grep -n -i -v -p "reference" ${dataset_id}/${dataset_id}_all.fasta > ${dataset_id}/${dataset_id}.fasta
    rm ${dataset_id}/${dataset_id}_all.fasta

    qualimap bamqc -bam ${dataset_id}/${dataset_id}.bam -outdir stats
    mv stats/genome_results.txt ${dataset_id}/${dataset_id}_stats_cov.txt
    rm -r stats
    """
    }
}


workflow.onComplete {
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command Line:     $workflow.commandLine"
    log.info "Container:        $workflow.container"
    log.info "Duration:     $workflow.duration"
    log.info "Output Directory: $params.output"
}
