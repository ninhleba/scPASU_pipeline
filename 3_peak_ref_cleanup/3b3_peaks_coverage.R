#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(argparser)

parser<-arg_parser(name="3b2_peaks_coverage.R",description="Calculate for the coverage of peaks within each TU") 

parser<-add_argument(
  parser,
  arg='--ref_file',
  short = '-r',
  type="character",
  default='peak_universe_updated.txt',
  help="Enter the path to the reference TU-assigned peak file from step 3a. Format: path/to/ref_file. Default: peak_universe_updated.txt")

parser<-add_argument(
  parser,
  arg='--counts_file',
  short = '-c',
  type="character",
  default='peak_count.rds',
  help="Enter the path to the peak count file from feature_counts.R. Format: path/to/counts_file. Default: peak_count.rds")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

args <- parse_args(parser)

ref_file=args$ref_file
counts_file=args$counts_file
out<-args$outdir

#Read counts
counts<-readRDS(counts_file)
counts_mat<-counts$counts %>% as.data.frame()
colnames(counts_mat)<-'counts'

# Add annotation cols
tu<-strsplit(rownames(counts_mat),split=':') %>% sapply(.,'[[',1)
gene<-strsplit(rownames(counts_mat),split=':') %>% sapply(.,'[[',2)
peak<-strsplit(rownames(counts_mat),split=':') %>% sapply(.,'[[',3)

counts_mat$tu<-tu
counts_mat$gene<-gene
counts_mat$peak<-peak
counts_mat$final_annotation<-rownames(counts_mat)

# Read ref 
ref<-fread(ref_file)

# Subset and sort ref 
idx<-match(rownames(counts_mat),ref$final_annotation)
ref<-ref[idx,]

# Sanity check
identical(ref$final_annotation,rownames(counts_mat))

# Add Peak coordinates as well
counts_mat$chr<-ref$chr
counts_mat$start<-ref$start
counts_mat$end<-ref$end
counts_mat$strand<-ref$strand

# Finally add peakID which will help match retained peaks later on (note: this col is labelled as 'peak' in ref table)
counts_mat$peakID<-ref$peak 

# Create separate counts matrix for P0 genes
p0<-which(counts_mat$peak=='P0')
counts_p0<-counts_mat[p0,]
counts_p0$pcov_pct<-100


# Remove P0 genes from ref
ref2<-ref[-p0,]
tu<-strsplit(ref2$tu_anno,split=':') %>% sapply(.,'[[',1)
tu<-unique(tu)

# Proceed with multi TU genes
counts_mat2<-counts_mat[-p0,]

updated_mat<-data.frame()

for(i in 1:length(tu))
{
r<-which(counts_mat2$tu==tu[i])
sub<-counts_mat2[r,]
cov<-sum(sub$counts)
pcov_pct<-(sub$counts/cov)*100
sub$pcov_pct<-pcov_pct

if(nrow(updated_mat)==0)
{
  updated_mat<-sub
}else{
  updated_mat<-rbind(updated_mat,sub)
}
}

# Merge matrices

updated_mat<-rbind(counts_p0,updated_mat)


#out<-dir(counts_file)
fname<-gsub('.rds','_updated.rds',basename(counts_file))
saveRDS(updated_mat,paste0(out,fname))

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
