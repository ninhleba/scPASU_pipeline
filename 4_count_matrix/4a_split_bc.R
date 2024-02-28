#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

parser<-arg_parser(name="4a_split_bc.R",description="Split each barcode into a separate file")

parser<-add_argument(
  parser,
  arg='--bc_dir',
  short = '-b',
  type="character",
  default='./',
  help="Directory of the desired barcodes of each sample. Default=./")

parser<-add_argument(
  parser,
  arg='--samples',
  short = '-s',
  type="character",
  default=NA,
  help="Sample to split barcodes. Sample names separated by comma. For example: Sample1,Sample2,Sample3")

parser<-add_argument(
  parser,
  arg='--subset',
  short = '-c',
  type="character",
  default='all',
  help="Compartment of cluster of interest. Default=all")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Directory to store split barcodes. Default=./")

parser<-add_argument(
  parser,
  arg='--num_folder',
  short = '-n',
  type="numeric",
  default=NA,
  help="Number of folders to evenly sort barcodes into for parallel purposes.")

args <- parse_args(parser)

bc_dir <- args$bc_dir
samples <- args$samples %>% strsplit(split=',') %>% unlist
subset <- args$subset
outdir <- args$outdir
num_fold <- args$num_fold

for (s in samples) {
  file_name = paste0(bc_dir,s,'_',subset,'_barcodes.tsv')
  bc_file <- fread(file_name, header = FALSE)
  bc <- bc_file$V1
  out <- paste0(outdir,'/',s,'/')
  if(!dir.exists(out)) {dir.create(out, recursive = TRUE)}
  lapply(1:length(bc), function(x){write.table(as.data.frame(bc[x]),paste0(out,bc[x],'.tsv'),col.names=FALSE,quote=FALSE,row.names = FALSE)})
  
  if (is.na(num_fold) == FALSE) {
    setwd(out)
    all_bc_files <- list.files(out) %>% unlist
    split_length <- length(all_bc_files) %/% num_fold
    for (i in 1:num_fold) {
      subfold <- paste0('folder',i)
      if(!dir.exists(subfold)) {dir.create(subfold, recursive = TRUE)}
      start <- (i - 1) * split_length + 1
      end <- i * split_length
      if (i == num_fold) {end <- length(all_bc_files)}
      for (j in start:end) {file.rename(from=all_bc_files[j], to=paste0(subfold,'/',all_bc_files[j]))}
    }
  }
}


