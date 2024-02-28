#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(Seurat)
library(stringr)

# Local test
counts_dir <- '~/Desktop/SC296_scPASU_run/outputs/cancercell/4d_merged_PeakMat/'
fprefix <- 'cancercell'
seurat_obj_path <- '~/Desktop/SC296_scPASU_run/iTAP scRNA (Ting)/SC93_epi_anno.RDS'
outdir <- '~/Desktop/SC296_scPASU_run/outputs/cancercell/5a_merged_per_sample_counts/'

if(!dir.exists(outdir))
{dir.create(outdir, recursive = TRUE)}

f<-list.files(counts_dir,full.names = TRUE)
f <- list.files(f,full.names = TRUE)

peak_counts<-lapply(f,read.table,header=TRUE,sep='\t')

# Modify cell barcodes if provide Seurat object
if (is.na(seurat_obj_path)==FALSE) {
  cat('Modifying cell barcodes to match Seurat when integrated \n')
  z <- readRDS(seurat_obj_path)
  meta <- z@meta.data
  write.table(meta,paste0(outdir,fprefix,'_meta.txt'),col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')
  meta$sample_idx <-as.numeric(sub(".*_([0-9]+)$", "\\1", row.names(meta)))
  sample_idx <- unique(meta$sample_idx)
  names(sample_idx) <- meta[match(sample_idx,meta$sample_idx),'sample']
  
  for (i in 1:length(peak_counts)) {
    for (s in names(sample_idx)) {
      if (all(str_detect(colnames(peak_counts[[i]]),paste0('^',s))) == TRUE) {
        colnames(peak_counts[[i]]) <- sub(".*_(.*)", "\\1",colnames(peak_counts[[i]]))
        colnames(peak_counts[[i]]) <- paste0(colnames(peak_counts[[i]]),'_',sample_idx[names(sample_idx)==s])
        }
      }
  }
}

# Merge samples
merged_counts<-do.call(cbind,peak_counts)
colnames(merged_counts) <- gsub('\\.','-',colnames(merged_counts))

stopifnot(identical(sort(colnames(merged_counts)),sort(substr(row.names(subset(meta,cell_subtype == 'Cancer_cell')),4,
                                                              nchar(row.names(subset(meta,cell_subtype == 'Cancer_cell')))))))

write.table(merged_counts,paste0(outdir,fprefix,'_counts.txt'),row.names = TRUE, col.names = TRUE, quote = FALSE, sep='\t')
