#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(reshape)
library(stringr)
library(data.table)
library(edgeR)
library(DEXSeq)
library(Seurat)

indir <- '~/Desktop/SC296_scPASU_run/outputs/cancercell/5a_merged_per_sample_counts/'
script_dir <- '~/Desktop/SC296_scPASU_run/Scripts/cancercell/5_APA_testing/'
fprefix <- 'cancercell'
counts_file <- 'cancercell_counts.txt'
meta_file <- 'cancercell_meta.txt'
peak_ref_file <- '~/Desktop/SC296_scPASU_run/outputs/cancercell/3h_fragmented_peaks_to_merge/cancercell_final_peak_universe_updated.txt'
outdir <- '~/Desktop/SC296_scPASU_run/outputs/cancercell/5b_APA_testing/iAP_G0_v_iAP_G4/'

source(paste0(script_dir,'scPASU_functions_for_differential_testing.R'))

# Read merged counts
merged_counts<-read.delim(paste0(indir,counts_file))
colnames(merged_counts) <- gsub('\\.','-',colnames(merged_counts))

# Read meta data [save meta data from Seurat object - ensure the set of cells match that of the count matrix]
meta<-read.table(paste0(indir,meta_file))
meta$cell <- substr(meta$cell,4,nchar(meta$cell))
row.names(meta) <- meta$cell
meta <- subset(meta, cell_subtype == 'Cancer_cell') ### Subset meta data by the cells in your compartment only

# Peak ref is jtu$join and final_annotation is the peak column 
peak_ref<-read.table(paste0(peak_ref_file),header=TRUE,sep='\t')

# Replace peak name column with final annotation as this is more meaningful name 
cat('Use final annotations as the peak names \n')
col<-which(colnames(peak_ref)=='final_annotation')
peak_ref$peak<-peak_ref[,col]
peak_ref<-peak_ref[,-col]

## Remove P0 peaks ##
cat('Remove P0 peaks \n')
plist<-strsplit(peak_ref$peak,split=':') %>% sapply(.,'[[',3)
rem<-which(plist=='P0')
peak_ref<-peak_ref[-rem,]

### DEXSeq test ###

apa <- stratify_matrix(merged_counts=merged_counts, meta=meta, vars=c('genotype'),cutoff_pct = 0, min_cell_per_group = 0)
if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)}

###### G0 v. G4 ######

inputs <- create_test_inputs(test='APA',apa,ident1=c('iAP_G0'),ident2=c('iAP_G4'),min_cell_expr_pct=10,
                             expr.thres=1,pseudo.count=1,replicate='random',nrep=3,p=0.7)

colnames(inputs) <- c(paste0('iAP_G0_rep',1:3),paste0('iAP_G4_rep',1:3))

### DEXSeq test ###
DEXseq_res <- APA_DEXseq_test(ident1='iAP_G0',ident2='iAP_G4',inputs,min_peak=2,ncpu=4,dispersion.plot.save=TRUE,peak_ref=peak_ref,outdir=outdir)

### T Test ###
DEXseq_res <- APA_TTest(DEXseq_res,ident1='iAP_G0',ident2='iAP_G4')

## Check significance ##
DEXseq_res <- callsig(DEXseq_res,ident1='iAP_G0',ident2='iAP_G4',delta_frac_thres=0.1, padj_thres=0.01,
                      l2fc_frac_thres=log2(1.5),mean_frac_thres=0.05)

saveRDS(DEXseq_res,paste0(outdir,'iAP_G0_v_iAP_G4_APAtest.res.rds'))

### Generate raw counts files and UCSC bedGraph tracks for significant APA genes
## Pairwise
pairwise_res <- list()
i = 1
for (comp in c('iAP_G0_v_iAP_G4')) {
  DEXseq_res <- readRDS(paste0(outdir,comp,'_APAtest.res.rds'))
  pairwise_res[[i]] <- DEXseq_res$res
  i = i + 1
}

names(pairwise_res) <- c('iAP_G0_v_iAP_G4')

res <- pairwise_res$iAP_G0_v_iAP_G4
table(res$int_sig)

bedtracks_list <- list()
i = 1
j = 1
for (comp in c('iAP_G0_v_iAP_G4')) {
  ident1 <- sapply(str_split(comp,pattern='_v_',n = 2),`[`,1)
  ident2 <- sapply(str_split(comp,pattern='_v_',n = 2),`[`,2)
  bedtracks_list[[i]] <- as.data.frame(pairwise_res[[j]])[,c('chr','start','end',paste0('mean_',ident1,'_frac'),'strand','int_sig')]
  names(bedtracks_list)[i] <- paste0(comp,'-',ident1)
  i = i + 1
  bedtracks_list[[i]] <- as.data.frame(pairwise_res[[j]])[,c('chr','start','end',paste0('mean_',ident2,'_frac'),'strand','int_sig')]
  names(bedtracks_list)[i] <- paste0(comp,'-',ident2)
  i = i + 1
  j = j + 1
}

# Create bedGraph tracks
for (i in 1:length(bedtracks_list)) {
  bedGraph <- bedtracks_list[[i]] 
  ident <- names(bedtracks_list[i])
  plus <- bedGraph[which(bedGraph$strand == '+'),]
  plus <- plus[,!(names(plus) %in% c('strand','int_sig'))]
  minus <- bedGraph[which(bedGraph$strand == '-'),]
  minus <- minus[,!(names(minus) %in% c('strand','int_sig'))]
  write.table(plus,paste0(outdir,ident,'_plus.bedGraph'),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  write.table(minus,paste0(outdir,ident,'_minus.bedGraph'),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}

for (comp in c('iAP_G0_v_iAP_G4')) { 
  DEXseq_res <- readRDS(paste0(outdir,comp,'_APAtest.res.rds'))
  raw.count <- DEXseq_res$counts.raw
  write.table(raw.count,paste0(outdir,comp,'_raw.count.txt'),col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
  res <- pairwise_res[[comp]]
  res <- res[,-c('dexseq_sig','ttest_sig','all_tests_sig')]
  write.table(res,paste0(outdir,comp,'_res.txt'),col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
}
