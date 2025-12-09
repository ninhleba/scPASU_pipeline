#!/bin/bash

dir='/Users/ninhle/Desktop/Research/scPASU_pipeline_runs/Ureter10_scPASU_run/outputs/differentiation_stage_cellranger_peakcount/'
fprefix=${dir}WUI_APA_genes_UTR3_shortening_ranges

echo -e "SeqName\tseqlen" > ${fprefix}_seqlen.txt
awk '{ OFS="\t" } 
/^>/ {
    if (header != "") print header, seq_length;
    header = substr($0, 2);
    seq_length = 0;
    next;
} {
    # Remove all spaces before calculating length
    gsub(/ /, "", $0)
    seq_length += length($0);
} END {
    if (header != "") print header, seq_length;
}' ${fprefix}.fa >> ${fprefix}_seqlen.txt
