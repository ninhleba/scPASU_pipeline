#!/usr/bin/env Rscript
options(scipen = 999, digits = 20)

library(dplyr)
library(data.table)
library(argparser)

parser<-arg_parser(name="hexamer_check_bed_input.R",description="Prepare bed input for hexamer checking")

parser<-add_argument(
  parser,
  arg='--ref_file',
  short = '-r',
  type="character",
  default='peak_universe_updated.txt',
  help="Enter the path to the reference TU-assigned peak file from the previous step. Format: path/to/ref_file. Default: peak_universe_updated.txt")

parser<-add_argument(
  parser,
  arg='--extn',
  short = '-e',
  type="character",
  default=NA,
  help="Enter direction and number of base pair to extend, separated by a comma. For example, 5,100")

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

args <- parse_args(parser)

ref_file=args$ref_file
out=args$outdir
fprefix=args$file_prefix
extn=args$extn

# Read inputs
ref<-fread(ref_file)
ref <- ref[,c('chr','pr_start','pr_end','final_annotation','polya_count','strand')]

if (is.na(extn) == FALSE) {
  extn <- strsplit(extn, split = ',') %>% unlist
  direction <- extn[1]
  ext_len <- extn[2] %>% as.numeric
  cat(paste('Extend search space by',ext_len,'bp',direction,'prime \n'))
  if (direction == '3') {
    ref$pr_end <- ifelse(ref$strand=='+',ref$pr_end+ext_len,ref$pr_end)
    ref$pr_start <- ifelse(ref$strand=='-',ref$pr_start-ext_len,ref$pr_start)
  } else if (direction == '5') {
    ref$pr_start <- ifelse(ref$strand=='+',ref$pr_start-ext_len,ref$pr_start)
    ref$pr_end <- ifelse(ref$strand=='-',ref$pr_end+ext_len,ref$pr_end)
  } else {stop('Direction has to be either 3 or 5 \n')}
}

write.table(ref,paste0(out,fprefix,'_hexamer_check_input.bed'),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')