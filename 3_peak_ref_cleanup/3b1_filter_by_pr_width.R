#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(argparser)

parser<-arg_parser(name="3b1_filter_by_pr_width.R",description="Filter peaks by polyA site width")

parser<-add_argument(
  parser,
  arg='--ref_file',
  short = '-r',
  type="character",
  default='peak_universe_updated.txt',
  help="Enter the path to the reference TU-assigned peak file from step 3a. Format: path/to/ref_file. Default: peak_universe_updated.txt")

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
  arg='--max_width',
  short = '-m',
  default=1000,
  type='numeric',
  help="Maximum polyA site width to retain. Default = 1000")

args <- parse_args(parser)

ref_file=args$ref_file
out=args$outdir
fprefix=args$file_prefix
max_width=as.numeric(args$max_width)

# Read inputs
ref<-fread(ref_file)
ref_filtered<-ref[which(ref$pr_width<=max_width),]

write.table(ref_filtered,paste0(out,fprefix,'_pr_filtered_peak_ref.txt'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
