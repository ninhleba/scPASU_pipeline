#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(argparser)

parser<-arg_parser(name="create_bed_file.R",description="Create bed file from the reference peak file.")

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
  arg='--name_col',
  short = '-n',
  type="character",
  default='final_annotation',
  help="Enter the column name of the annotation column. Default: final_annotation")

parser<-add_argument(
  parser,
  arg='--chrs',
  short = '-c',
  type="character",
  default='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY',
  help="Enter the chromosomes in the desired order, separated by comma.")

args <- parse_args(parser)

ref_file=args$ref_file
out=args$outdir
fprefix=args$file_prefix
name_col=args$name_col
chrs <- args$chrs %>% strsplit(split = ',') %>% unlist

if(!dir.exists(out))
{dir.create(out, recursive = TRUE)}

ref<-fread(ref_file)
ref<-as.data.frame(ref)

cols<-c('chr','start','end',name_col)
bed<-ref[,cols]
bed$score<-rep(0,nrow(bed))
bed$strand<-ref$strand

write.table(bed, paste0(out,fprefix,'.bed'),sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)

bed_plus <- bed[which(bed$strand == '+'),]
bed_minus <- bed[which(bed$strand == '-'),]

write.table(bed_plus, paste0(out,fprefix,'_plus.bed'),sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(bed_minus, paste0(out,fprefix,'_minus.bed'),sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
