#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(magrittr)
library(argparser)

parser<-arg_parser(name="3b2_plot_peak_cov.R",description="Plot the peak coverage computed from step 3b2.")

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
  arg='--hist_args',
  short = '-a',
  type="character",
  default='11,8,100,5000',
  help="Enter the arguments for histograms, separated by comma. Order: width, height, x_lim, y_lim")

parser<-add_argument(
  parser,
  arg='--script_dir',
  short = '-s',
  type="character",
  default='./',
  help="Enter the script directory path. Format: path/to/script/dir/ Default= ./")

args <- parse_args(parser)

pcov_file=args$pcov_file
out=args$outdir
fprefix=args$file_prefix
hist_args=args$hist_args %>% strsplit(split = ',') %>% unlist %>% as.numeric
w=hist_args[1]
h=hist_args[2]
x=hist_args[3]
y=hist_args[4]
script_dir=args$script_dir

source(paste0(script_dir,'scPASU_functions.R'))

pcov_mat<-readRDS(pcov_file)

dat<-pcov_mat$pcov_pct

options(bitmapType='cairo')
pdf(paste0(out,fprefix,'_peak_cov_pct.pdf'),width=w, height=h)
get_histogram(dat=dat,binw=1,col='royal blue',fill='salmon', x_lim=c(0,100),y_lim=c(0,5000),title='Peak read count percentage',x_lab='peak_read_pct')
dev.off()

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')

