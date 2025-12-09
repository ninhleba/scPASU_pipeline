library(dplyr)
library(goldmine)
library("readxl")

outdir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5b_APA_testing/differentiation_stage_cellranger_peakcount/'
peak_ref_file <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3h_fragmented_peaks_to_merge/Ureter10_urothelial_final_peak_universe_updated.txt'

setwd(outdir)

ref <- read_excel('wui_apa_gene_UTR3_shorterning_AT.xls')
ref <- read_excel('wui_apa_gene_UTR3_lengthening_AT.xls')
ref1 <- ref[is.na(ref$Delete),]
ref1$gene %>% unique() %>% length()
ref2 <- split(ref1,ref1$gene)
ref2 <- lapply(ref2, function(x) {
  chr <- unique(x$chr)
  gene <- unique(x$gene)
  score <- 1
  strand <- unique(x$strand)
  order <- order(x$pr_start, decreasing = F)
  if (!identical(x,x[order,])) {
    cat(gene,'needs sorting \n')
    x <- x[order,]
    }
  ranges <- data.frame(chr = chr, start = x$pr_start[1], end = x$pr_end[nrow(x)], gene = paste0(gene,'_UTR3'), score = score, strand = strand)
  return(ranges)
})

ref2 <- do.call('rbind',ref2)
write.table(ref2,'WUI_APA_genes_UTR3_shorterning_ranges.bed',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

# bedtools getfasta -fi $ref_genome -bed $bed_input -s -nameOnly -fo $fa_output
