#https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

outdir <- "~/Desktop/SC296_scPASU_run/outputs/cancercell/DEG_pseudobulk/iAP_G0_v_iAP_G4/"
script_dir <- '~/Desktop/SC296_scPASU_run/Scripts/cancercell/5_APA_testing/'
seurat_obj <- '~/Desktop/SC296_scPASU_run/iTAP scRNA (Ting)/SC93_epi_anno.RDS'

source(paste0(script_dir,'scPASU_functions_for_differential_testing.R'))
if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
z <- readRDS(seurat_obj)
counts <- z@assays$RNA@counts
meta <- z@meta.data
meta <- subset(meta, cell_subtype == 'Cancer_cell')

counts_by_genotype <- stratify_matrix(merged_counts = counts, meta = meta, vars = 'genotype',
                                     cutoff_pct = 1, min_cell_per_group = 10)

##### iAP_G0 v. iAP_G4
inputs <- create_test_inputs(test = 'APA', all_groups = counts_by_genotype, ident1=c('iAP_G0'),
                             ident2=c('iAP_G4'), APA.feature.filter = FALSE, 
                             replicate = 'random', nrep = 3 ,p = 0.7)

colnames(inputs) <- c(paste0('iAP_G0_rep',1:3),paste0('iAP_G4_rep',1:3))

sig_genes <- DEG_DESeq2_pseudobulk(inputs = inputs, comp = c('iAP_G4','iAP_G0'), test.used = 'LRT',
                                   padj_thres=0.01,l2fc_frac_thres=log2(1.5),outdir)

