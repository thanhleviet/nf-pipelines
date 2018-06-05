#!/usr/bin/env bash
nextflow snps.nf \
    --input "/home/ubuntu/data/seq_ref_lab_2018/results/results_may_12/speciation/species_reads/Escherichia_coli" \
    --output "/home/ubuntu/data/seq_ref_lab_2018/results/results_may_12/E.coli_snps" \
    --ref "/home/ubuntu/data2/genomes/NC_000913.fa" \
    -with-report phenix.html \
    -bg \
    -resume
