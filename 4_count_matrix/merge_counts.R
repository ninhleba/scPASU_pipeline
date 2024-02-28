#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

parser<-arg_parser(name="merge_counts.R",description="Merge counts of each cell into a matrix.")

parser<-add_argument(
  parser,
  arg='--counts_dir',
  short = '-d',
  type="character",
  default='./',
  help="Enter the path where the counts of all the cells are. Format: path/to/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='scPASU',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

args <- parse_args(parser)

counts_path<-args$counts_dir
file_prefix<-args$file_prefix
outdir<- args$outdir

if(!dir.exists(outdir))
{dir.create(outdir, recursive = TRUE)}

setwd(counts_path)

f<-list.files(counts_path,pattern = '.rds',full.names = TRUE)
dat<-lapply(f,readRDS)

# The outputs of SubRead are stored as list so dat is list of lists - extracting features from one of the countr matrix
counts<-rownames(dat[[1]]$counts)

 for(i in 1:length(dat))
 {
  cnt<-dat[[i]]$counts
  counts<-cbind(counts,cnt)
  
  }

colnames(counts)<-gsub('.bam','',colnames(counts))

counts<-as.data.frame(counts)
counts<-counts[,-1]

write.table(counts,paste0(outdir,file_prefix,'_counts.txt'),row.names=TRUE,col.names=TRUE,sep='\t',quote=FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
 
