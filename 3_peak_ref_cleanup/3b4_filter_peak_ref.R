#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(argparser)

parser<-arg_parser(name="3b3_filter_peak_ref.R",description="Filter peaks on multi-peak TU by coverage")

parser<-add_argument(
  parser,
  arg='--ref_file',
  short = '-r',
  type="character",
  default='peak_universe_updated.txt',
  help="Enter the path to the reference TU-assigned peak file from step 3a. Format: path/to/ref_file. Default: peak_universe_updated.txt")

parser<-add_argument(
  parser,
  arg='--pcov_file',
  short = '-p',
  type="character",
  default='peak_count_updated.rds',
  help="Enter the path to the peak count file with peak coverage added from step 3b2. Format: path/to/pcov_file. Default: peak_count_updated.rds")

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
  arg='--min_cov',
  short = '-m',
  default=10,
  type='numeric',
  help="Minimum coverage to retain peaks from multi-peak TUs. Default = 10")

args <- parse_args(parser)

ref_file=args$ref_file
cov_mat_file=args$pcov_file
out=args$outdir
fprefix=args$file_prefix
mincov=as.numeric(args$min_cov)

ref<-fread(ref_file)
cov_mat<-readRDS(cov_mat_file)
if (length(which(is.na(cov_mat$pcov_pct))) != 0) {
  cov_mat<-cov_mat[-which(is.na(cov_mat$pcov_pct)),]
}
sub<-cov_mat[which(cov_mat$pcov_pct>=mincov),]

select<-match(sub$peakID,ref$peakID)
ref_filtered<-ref[select,]

write.table(ref_filtered,paste0(out,fprefix,'_cov_filtered_peak_ref.txt'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
