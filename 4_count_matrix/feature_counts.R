#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(Rsubread)
library(argparser)

parser<-arg_parser(name="feature_counts.R",description="Count reads from bam file for peaks")

parser<-add_argument(
  parser,
  arg='--bam_file',
  short = '-b',
  type="character",
  default='file.bam',
  help="Enter the path to the reference TU-assigned peak file from step 3a. Format: path/to/ref_file")

parser<-add_argument(
  parser,
  arg='--ref_file',
  short = '-r',
  type="character",
  default='peak_universe_updated.saf',
  help="Enter the path to the reference TU-assigned peak file from step 3a. Format: path/to/ref_file. Default: peak_universe_updated.saf")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='scPASU',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--cores',
  short = '-c',
  default=30,
  type='numeric',
  help="Number of cores to use when processing samples in parallel. Defulat = 30")

parser<-add_argument(
  parser,
  arg='--is_GTF',
  short = '-i',
  type="character",
  default='no',
  help="Whether or not the reference TU-assigned peak file is GTF format. Default: no")

args <- parse_args(parser)

bam_file=args$bam_file
ref=args$ref_file
out=args$outdir
prefix=args$file_prefix
cores=args$cores
isGTF=args$is_GTF

if(isGTF=='yes')
{
  counts<-featureCounts(bam_file, annot.ext=ref, isGTFAnnotationFile=TRUE, GTF.attrType='gene_name', strandSpecific=1, largestOverlap = T, nthreads = cores)
  fname<-paste0(out,prefix,'_gene_count.rds')
}else{
  counts<-featureCounts(bam_file, annot.ext = ref, isGTFAnnotationFile = FALSE,strandSpecific=1, largestOverlap = T, nthreads = cores)
  fname<-paste0(out,prefix,'_peak_count.rds')
}

saveRDS(counts, fname)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
