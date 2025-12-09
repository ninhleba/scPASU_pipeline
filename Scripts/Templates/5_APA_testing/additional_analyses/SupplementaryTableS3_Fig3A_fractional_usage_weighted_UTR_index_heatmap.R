library(ggplot2)
library(ggrepel)
library(dplyr)
library(goldmine)
library(gtools)
library(parallel)
library(RColorBrewer)

outdir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5b_APA_testing/differentiation_stage_cellranger_peakcount/'
peak_ref_file <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3h_fragmented_peaks_to_merge/Ureter10_urothelial_final_peak_universe_updated.txt'
turef_file <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3f_merge_two_prongs/genes_newflankupdated.rds'

turef <- readRDS(turef_file)
peak_ref <- read.delim(peak_ref_file)
peak_ref <- peak_ref %>% mutate(key = paste(chr,start,end,strand,sep='_'))
which(peak_ref$key %>% duplicated())
setwd(outdir)

files <- list.files()
for (file in grep('res.txt$', files, value = TRUE)) {
  assign(gsub('_res.txt','',file),read.delim(file))
}

unique(subset(basal_v_intermediate, int_sig == TRUE)$gene_name) %>% length() #315
unique(subset(intermediate_v_umbrella, int_sig == TRUE)$gene_name) %>% length() #525
unique(subset(umbrella_v_basal, int_sig == TRUE)$gene_name) %>% length() #622

union_apa_genes <- subset(basal_v_intermediate, int_sig == TRUE)$gene_name
union_apa_genes <- c(union_apa_genes,subset(intermediate_v_umbrella, int_sig == TRUE)$gene_name)
union_apa_genes <- c(union_apa_genes,subset(umbrella_v_basal, int_sig == TRUE)$gene_name)
union_apa_genes <- unique(union_apa_genes)

union_apa_peaks <- subset(basal_v_intermediate, int_sig == TRUE)$peak
union_apa_peaks <- c(union_apa_peaks,subset(intermediate_v_umbrella, int_sig == TRUE)$peak)
union_apa_peaks <- c(union_apa_peaks,subset(umbrella_v_basal, int_sig == TRUE)$peak)
union_apa_peaks <- unique(union_apa_peaks)

for (compartment in c('basal','intermediate','umbrella')) {
  res_plus_file <- paste0(compartment,'_merged_plus.bedGraph')
  res_minus_file <- paste0(compartment,'_merged_minus.bedGraph')
  
  res_plus <- read.delim(res_plus_file, header = F)
  res_plus$strand <- '+'
  res_minus <- read.delim(res_minus_file, header = F)
  res_minus$strand <- '-'
  res <- rbind(res_plus, res_minus)
  colnames(res) <- c('chr','start','end',paste0(compartment,'_fractional_usage'),'strand')
  res <- res %>% mutate(key = paste(chr,start,end,strand,sep='_'))
  res$final_annotation <- peak_ref$final_annotation[match(res$key,peak_ref$key)]
  res$key <- NULL
  
  assign(compartment,res)
  rm(res_plus_file,res_minus_file,res_plus,res_minus,res)
}

fractional_usage_df <- cbind(basal[,c('final_annotation','basal_fractional_usage')],intermediate[,c('final_annotation','intermediate_fractional_usage')],
                             umbrella[,c('final_annotation','umbrella_fractional_usage')])
row.names(fractional_usage_df) <- fractional_usage_df$final_annotation
fractional_usage_df$gene <- peak_ref$gene[match(fractional_usage_df$final_annotation,peak_ref$final_annotation)]
fractional_usage_df <- subset(fractional_usage_df, final_annotation %in% union_apa_peaks)  
fractional_usage_df[,c(1,3,5,7)] <- NULL

# Inputs for gProfiler

write.table(as.data.frame(unique(subset(basal_v_intermediate, int_sig == TRUE)$gene_name)), 'b_v_i_apa_genes.txt',row.names = F, col.names = F, quote = F)
write.table(as.data.frame(unique(subset(intermediate_v_umbrella, int_sig == TRUE)$gene_name)), 'i_v_u_apa_genes.txt',row.names = F, col.names = F, quote = F)
write.table(as.data.frame(unique(subset(umbrella_v_basal, int_sig == TRUE)$gene_name)), 'u_v_b_apa_genes.txt',row.names = F, col.names = F, quote = F)
write.table(as.data.frame(union_apa_genes),'union_apa_genes.txt',row.names = F, col.names = F, quote = F)

### Weighted usage index

#### Filter for genes with UTR3 and flank peaks only 
## Filter for genes with only UTR3 and flank peaks. Annotate what UTR3 peaks belong to the same UTR3 using the UTR3 slot from turef
peakref_spl <- split(peak_ref, peak_ref$gene)
utr3flankonly <- mclapply(peakref_spl, function(x) {
  if (all(x$final_classification %in% c('UTR3','Flank'))) {
    strand <- unique(x$strand)
    gene <- unique(x$gene)

    if (gene %in% turef$utr3$name) {
      utr3_gene <- turef$utr3[which(turef$utr3$name == gene)]
      if (strand == '+') {
        utr3_gene$mostdownstream <- 1:length(utr3_gene)==length(utr3_gene)
      } else if (strand == '-') {
        utr3_gene$mostdownstream <- 1:length(utr3_gene)==1
      }
    }

    utr3_peaks <- subset(x, final_classification == 'UTR3')
    num_utr3_peaks <- nrow(utr3_peaks)
    
    if (num_utr3_peaks == 0) {
      x$shared_UTR3 <- NA
      x$mostdownstreamUTR3 <- NA
      return(x)
    } else if (num_utr3_peaks == 1) {
      utr3_peaks$shared_UTR3 <- 'None'
      utr3_peaks <- makeGRanges(utr3_peaks,strand=T)
      ovl <- findOverlaps(utr3_peaks,utr3_gene) %>% as.data.frame()
      utr3_peaks$mostdownstreamUTR3 <- any(utr3_gene$mostdownstream[ovl$subjectHits])
      x$shared_UTR3 <- utr3_peaks$shared_UTR3[match(x$final_annotation,utr3_peaks$final_annotation)]
      x$mostdownstreamUTR3 <- utr3_peaks$mostdownstreamUTR3[match(x$final_annotation,utr3_peaks$final_annotation)]
      
      return(x)
    } else {
      utr3_peaks$shared_UTR3 <- ''
      utr3_peaks$mostdownstreamUTR3 <- ''
      utr3_peaks <- makeGRanges(utr3_peaks,strand=T)

      ovl <- findOverlaps(utr3_peaks,utr3_gene) %>% as.data.frame()
      ovl_spl <- split(ovl, ovl$subjectHits)
      
      ovl_spl <- lapply(ovl_spl, function(x) {
        x <- x[!duplicated(x$queryHits),]
        if (nrow(x) > 1) {return(x)}
      })
      
      ovl_spl <- ovl_spl[!unlist(lapply(ovl_spl,is.null))]
      if (length(ovl_spl) != 0) {
        l <- list()
        for (i in 1:length(ovl_spl)) {
          l[[i]] <- ovl_spl[[i]]
        }
        
        l <- l[!duplicated(lapply(l,function(x) {sort(x$queryHits)}))]
        sorted_order <- lapply(l,nrow) %>% unlist() %>% order(decreasing = T)
        l <- l[sorted_order]
    
        for (i in 1:length(l)) {
          utr3_idx <- l[[i]]$subjectHits %>% unique
          utr3_id <- utr3_gene$utr[utr3_idx]
          mostdownstream <- utr3_gene$mostdownstream[utr3_idx]
          utr3_peaks$shared_UTR3[l[[i]]$queryHits] <- ifelse(utr3_peaks$shared_UTR3[l[[i]]$queryHits] == '',
                                                             utr3_id,
                                                             paste(utr3_peaks$shared_UTR3[l[[i]]$queryHits],utr3_id,sep = ','))
          utr3_peaks$mostdownstreamUTR3[l[[i]]$queryHits] <- ifelse(utr3_peaks$mostdownstreamUTR3[l[[i]]$queryHits] == '',
                                                                    mostdownstream,
                                                                    paste(utr3_peaks$mostdownstreamUTR3[l[[i]]$queryHits],mostdownstream,sep=','))
        }
      }
      
      utr3_peaks <- as.data.frame(utr3_peaks)
      utr3_peaks$width <- NULL
      colnames(utr3_peaks)[1] <- 'chr'
      
      num_loners <- length(which(utr3_peaks$shared_UTR3 == ''))
      if (num_loners != 0) {
        utr3_peaks$shared_UTR3[which(utr3_peaks$shared_UTR3 == '')] <- 'None'
        utr3_peaks$mostdownstreamUTR3[which(utr3_peaks$shared_UTR3 == 'None')] <- NA
      }
      
      x$shared_UTR3 <- utr3_peaks$shared_UTR3[match(x$final_annotation,utr3_peaks$final_annotation)]
      x$mostdownstreamUTR3 <- utr3_peaks$mostdownstreamUTR3[match(x$final_annotation,utr3_peaks$final_annotation)]
      return(x)
    }
  }
}, mc.cores = 4)

utr3flankonly <- utr3flankonly[!unlist(lapply(utr3flankonly,is.null))]

## UTR3 peaks have to share at least one UTR3. If flank peaks are also there, the shared UTR3 has to be most downstream
wui_eligbile <- lapply(utr3flankonly, function(x) {
    gene <- unique(x$gene)
    utr3_peaks <- subset(x, final_classification == 'UTR3')
    flank_peaks <- subset(x, final_classification == 'Flank')
    num_utr3_peaks <- nrow(utr3_peaks)
    num_flank_peaks <- nrow(flank_peaks)
    if (num_utr3_peaks == 0 | (num_utr3_peaks == 1 & num_flank_peaks == 0)) {
      return(x)
    } else if (num_utr3_peaks == 1 & num_flank_peaks != 0) {
      if (utr3_peaks$mostdownstreamUTR3 == TRUE) {return(x)}
    } else if (num_utr3_peaks > 1) {
      all_utr3s <- paste0(utr3_peaks$shared_UTR3,collapse = ',') %>% strsplit(split = ',') %>% unlist %>% table
      if (!('None' %in% names(all_utr3s)) & (max(all_utr3s) == num_utr3_peaks)) {
        if (num_flank_peaks == 0) {
          return(x)
        } else {
          shared_utr3 <- names(all_utr3s)[which(all_utr3s == max(all_utr3s))]
          if (length(shared_utr3) > 1) {cat(shared_utr3,'\n')}
          all_utr3s <- paste0(utr3_peaks$shared_UTR3,collapse = ',') %>% strsplit(split = ',') %>% unlist
          mostdownstreamUTR3 <- paste0(utr3_peaks$mostdownstreamUTR3,collapse = ',') %>% strsplit(split = ',') %>% unlist %>% as.logical()
          names(mostdownstreamUTR3) <- all_utr3s
          mostdownstreamUTR3_status <- mostdownstreamUTR3[which(names(mostdownstreamUTR3) == shared_utr3)] %>% unique
          if (mostdownstreamUTR3_status) {return(x)}
        }
     }
  }
})

wui_ineligbile <- utr3flankonly[names(which(unlist(lapply(wui_eligbile,is.null))))]
wui_ineligbile <- do.call('rbind',wui_ineligbile)

wui_eligbile <- wui_eligbile[!unlist(lapply(wui_eligbile,is.null))]
wui_eligbile <- do.call('rbind',wui_eligbile)
wui_eligbile <- wui_eligbile[!grepl(':P0$',wui_eligbile$final_annotation),]
wui_eligbile$key <- NULL
write.table(wui_eligbile,'genes_wui_eligible.txt',sep = '\t',row.names = F, col.names = T, quote = F)
wui_eligbile <- read.delim('genes_wui_eligible.txt')

# Filter for wui eligible gene
fractional_usage_df <- cbind(basal[,c('final_annotation','basal_fractional_usage')],intermediate[,c('final_annotation','intermediate_fractional_usage')],
                             umbrella[,c('final_annotation','umbrella_fractional_usage')])
row.names(fractional_usage_df) <- fractional_usage_df$final_annotation
fractional_usage_df$gene <- peak_ref$gene[match(fractional_usage_df$final_annotation,peak_ref$final_annotation)]
fractional_usage_df <- subset(fractional_usage_df, gene %in% union_apa_genes)
fractional_usage_df[,c(1,3,5)] <- NULL
fractional_usage_df <- fractional_usage_df[row.names(fractional_usage_df) %in% wui_eligbile$final_annotation,]
fractional_usage_df$gene %>% unique %>% length #297

# positional rank
peakref_subset <- subset(peak_ref, final_annotation %in% row.names(fractional_usage_df))
peakref_subset_spl <- split(peakref_subset, peakref_subset$gene)
peakref_subset_spl <- lapply(peakref_subset_spl, function(x) {
  if(unique(x$strand) == '+') {
    x <- x[order(x$pr_start, decreasing = F),]
    region_length <- x$pr_start[nrow(x)] - x$pr_start[1]
    x$positional_rank <- (x$pr_start - x$pr_start[1])/region_length
  } else if (unique(x$strand) == '-') {
    x <- x[order(x$pr_end, decreasing = T),]
    region_length <- x$pr_end[nrow(x)] - x$pr_end[1]
    x$positional_rank <- (x$pr_end - x$pr_end[1])/region_length
  }
  return(x)
})
peakref_subset_spl <- do.call('rbind',peakref_subset_spl)
write.table(peakref_subset_spl,'apa_genes_wui_eligible_positional_ranks.txt',row.names = F,col.names = T,sep = '\t',quote = F)

fractional_usage_df$positional_rank <- peakref_subset_spl$positional_rank[match(row.names(fractional_usage_df),peakref_subset_spl$final_annotation)]

wui <- split(fractional_usage_df, fractional_usage_df$gene)
wui <- lapply(wui, function(x) {
  x <- x[,c(1,2,3)] * x$positional_rank
  weighted_UTR_index <- colSums(x)
  return(weighted_UTR_index)
})

wui <- do.call('rbind',wui)

wui[wui == 0] <- NA
missing_wui <- apply(wui,1,function(x) {if(any(is.na(x))) {return(unname(x[!is.na(x)]))}}) 
missing_wui <- do.call('rbind',missing_wui[!(lapply(missing_wui,is.null) %>% unlist)])
full_wui <- wui[!(row.names(wui) %in% row.names(missing_wui)),]

mp <- pheatmap::pheatmap(missing_wui,
                         clustering_distance_rows = 'correlation',
                         clustering_method="ward.D2",
                         cluster_cols = FALSE, cluster_rows = T,
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         main = 'U10 Urothelial Basal, Intermediate, Umbrella Weighted UTR index \n WUI-eligible APA genes with missing WUI')

mp_gm <- cutree(mp$tree_row,k = 2)

fp <- pheatmap::pheatmap(full_wui,
                         clustering_distance_rows = 'correlation',
                         clustering_method="ward.D2",
                         cluster_cols = FALSE, cluster_rows = T,
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         main = 'U10 Urothelial Basal, Intermediate, Umbrella Weighted UTR index \n WUI-eligible APA genes with full WUI')

fp_gm <- cutree(fp$tree_row,k=2)

shortening_gm_fp <- names(fp_gm)[which(fp_gm==1)]
shortening_gm_fp <- shortening_gm_fp[match(fp$tree_row$labels[fp$tree_row$order],shortening_gm_fp)]
shortening_gm_fp <- shortening_gm_fp[!is.na(shortening_gm_fp)]

shortening_gm_mp <- names(mp_gm)[which(mp_gm==2)]
shortening_gm_mp <- shortening_gm_mp[match(mp$tree_row$labels[mp$tree_row$order],shortening_gm_mp)]
shortening_gm_mp <- shortening_gm_mp[!is.na(shortening_gm_mp)]

shortening_gm <- c(shortening_gm_fp,shortening_gm_mp)
write.table(as.data.frame(shortening_gm),'WUI_APA_genes_UTR3_shortening.txt',row.names = F, col.names = F, quote = F)

colfunc <- colorRampPalette(c('#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603'))

sp <- pheatmap::pheatmap(wui[shortening_gm,],
                         color = colfunc(50),
                         na_col = 'white',
                         clustering_distance_rows = 'correlation',
                         clustering_method="ward.D2",
                         cluster_cols = F, cluster_rows = F,
                         show_colnames = F,
                         show_rownames = T,
                         main = 'WUI-eligible APA genes: Shortening')

ggsave('WUI_APA_genes_UTR3_shortening.png',
       height = 16.4, width = 9, units = 'in', limitsize = F, sp)

lengthening_gm_fp <- names(fp_gm)[which(fp_gm==2)]
lengthening_gm_fp <- lengthening_gm_fp[match(fp$tree_row$labels[fp$tree_row$order],lengthening_gm_fp)]
lengthening_gm_fp <- lengthening_gm_fp[!is.na(lengthening_gm_fp)]

lengthening_gm_mp <- names(mp_gm)[which(mp_gm==1)]
lengthening_gm_mp <- lengthening_gm_mp[match(mp$tree_row$labels[mp$tree_row$order],lengthening_gm_mp)]
lengthening_gm_mp <- lengthening_gm_mp[!is.na(lengthening_gm_mp)]

lengthening_gm <- c(lengthening_gm_fp,lengthening_gm_mp)
write.table(as.data.frame(lengthening_gm),'WUI_APA_genes_UTR3_lengthening.txt',row.names = F, col.names = F, quote = F)

lp <- pheatmap::pheatmap(wui[lengthening_gm,],
                         color = colfunc(50),
                         na_col = 'white',
                         clustering_distance_rows = 'correlation',
                         clustering_method="ward.D2",
                         cluster_cols = F, cluster_rows = F,
                         show_colnames = F,
                         show_rownames = T,
                         main = 'WUI-eligible APA genes: Lengthening')

ggsave('WUI_APA_genes_UTR3_lengthening.png',
       height = 13.3, width = 9, units = 'in', limitsize = F, lp)

sessionInfo()
