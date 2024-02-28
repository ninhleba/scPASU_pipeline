library(goldmine)
library(rtracklayer)
library(dplyr)
library(data.table)
library(gtools)
library(argparser)

parser<-arg_parser(name="create_turef.R",description="Create a TU reference from gtf annotation")

parser<-add_argument(
  parser,
  arg='--script_dir',
  short = '-d',
  type="character",
  default='./',
  help="Enter the path to dir of the scPASU_functions.R file. Format: path/to/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--gtf',
  short = '-g',
  type="character",
  default='genes.gtf',
  help="Enter the path to the cellranger gtf gene file. Format: path/to/gtf_file. Default: ./genes.gtf")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--chrs',
  short='-c',
  type='character',
  default=NA,
  help='Only process genes on these chromosomes. Default: 24 canonical chromosomes of human genome')

parser<-add_argument(
  parser,
  arg='--flank',
  short = '-f',
  default=5000,
  type='numeric',
  help="Flank region length. Default = 5000")

parser<-add_argument(
  parser,
  arg='--use_genesymbol',
  short = '-u',
  flag=TRUE,
  help="Flag set to collapse transcripts by gene symbol instead of Ensembl ID.")

parser<-add_argument(
  parser,
  arg='--final_annotation',
  short = '-a',
  type="character",
  default='gene_symbol',
  help="gene_symbol means using gene symbol for final TU annotation; gene.id means using whatever ID in gene.id column for final TU annotation, which depends on the flag --use_genesymbol. Default= gene_symbol")

parser<-add_argument(
  parser,
  arg='--isoform_rm',
  short = '-i',
  type="character",
  default=NA,
  help="Isoforms to remove from gtf annotation file. Format: Ensembl IDs separated by comma")

args <- parse_args(parser)

script_dir <- args$script_dir
gtf_file <- args$gtf
outdir <- args$outdir
chrs <- args$chrs
flank <- args$flank %>% as.numeric
use.genesymbol <-args$use_genesymbol
final_annotation <- args$final_annotation
isoform.rm <- args$isoform_rm

if (is.na(chrs) == TRUE) {
  chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
         'chr20','chr21','chr22','chrX','chrY')
} else {
  chrs <- chrs %>% strsplit(split=',') %>% unlist
}

if (is.na(isoform.rm) == FALSE) {
  cat('Remove these isoforms from gtf annotation file:',isoform.rm,'\n')
  isoform.rm <- isoform.rm %>% strsplit(split=',') %>% unlist
} else {isoform.rm <- NULL}

source(paste0(script_dir,'scPASU_functions.R'))

if(!dir.exists(outdir))
{dir.create(outdir, recursive = TRUE)}

### Create genes table ###
g<-rtracklayer::import(gtf_file)
gtf=as.data.frame(g)
cat('Reformat gtf annotation file...\n')
genes <- reformat_gtf(gtf, chrs = chrs, isoform.rm = isoform.rm)
cat('Generate TU reference...\n')
red_ens <- create_TU_from_gtf(genes,flank=flank,outdir=outdir,use.genesymbol=use.genesymbol,save=TRUE,final_tu_annot=final_annotation)

# Print session info
sessionInfo()
