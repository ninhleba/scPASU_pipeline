library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

comp <- list(c('basal','intermediate'),c('intermediate','umbrella'),c('umbrella','basal'))
apa_dir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5b_APA_testing/differentiation_stage_cellranger_peakcount/'
deg_dir <- '~/Desktop/Ureter10_scPASU_run/outputs/urothelial/5b_DEG_testing/differentiation_stage_cellranger_genecount/'
df <- matrix(ncol = 3, nrow = 0) %>% as.data.frame()
colnames(df) <- c('Comparison','Type','Num')

padj_thres = 0.01
l2fc_frac_thres = abs(log2(1.5))

for (i in 1:length(comp)) {
  apa <- fread(paste0(apa_dir,paste0(comp[[i]], collapse = '_v_'),'_res.txt'))
  apa$gene_symbol <- strsplit(apa$peak,split=':') %>% sapply(.,'[[',2)
  apa_sig <- apa[which(apa$int_sig == TRUE),]
  
  deg <- fread(paste0(deg_dir,paste0(comp[[i]], collapse = '_vs_'),'_LRT_all_genes.txt'))
  deg_sig <- deg %>% filter(padj < padj_thres & abs(log2FoldChange) >= l2fc_frac_thres)
  
  int_num <- intersect(unique(apa_sig$gene_symbol),deg_sig$gene) %>% length
  apa_num <- length(unique(apa_sig$gene_symbol))
  deg_num <- length(deg_sig$gene)
  total = apa_num + deg_num - int_num
  
  df <- rbind(df,data.frame(Comparison = rep(paste0(comp[[i]], collapse = ' v. '), times = 3),
                            Type = c('Differentially expressed genes','Overlap','Significant APA genes'),
                            Num = c(deg_num - int_num, int_num, apa_num - int_num)))
}

#colors <- c("#FC8D62","#8DA0CB","#66C2A5")
colors <- c("#abd9e9","#fdae61","#a50026")
names(colors) <- c('Differentially expressed genes','Overlap','Significant APA genes')

g <- ggplot(df, aes(x = Comparison, y = Num, fill = Type)) + 
  geom_bar(stat = "identity", color = 'black') + scale_fill_manual(values = colors) + 
  scale_y_continuous(breaks = seq(0,5000,500)) + 
  theme_classic()
g + theme(text = element_text(size = 20))     
ggsave(paste0(apa_dir,'DEG_APA_ovl_barplot.png'),width = 12, height = 8, unit = 'in')

for (i in 1:length(comp)) {
  df <- matrix(ncol = 4, nrow = 0) %>% as.data.frame()
  colnames(df) <- c('gene','DEG','APA','Both')
  
  apa <- fread(paste0(apa_dir,paste0(comp[[i]], collapse = '_v_'),'_res.txt'))
  apa$gene_symbol <- strsplit(apa$peak,split=':') %>% sapply(.,'[[',2)
  apa_sig <- apa[which(apa$int_sig == TRUE),]
  
  deg <- fread(paste0(deg_dir,paste0(comp[[i]], collapse = '_vs_'),'_LRT_all_genes.txt'))
  deg_sig <- deg %>% filter(padj < padj_thres & abs(log2FoldChange) >= l2fc_frac_thres)
  
  apa_gene <- unique(apa_sig$gene_name)
  deg_gene <- deg_sig$gene
  int_gene <- intersect(apa_gene,deg_gene)
  
  apa_gene_only <- setdiff(apa_gene,int_gene)
  deg_gene_only <- setdiff(deg_gene,int_gene)
  
  df <- rbind(df,data.frame(gene = c(int_gene,deg_gene_only,apa_gene_only),
                            DEG = c(rep(TRUE, times = length(int_gene)), 
                                    rep(TRUE, times = length(deg_gene_only)),
                                    rep(FALSE, times = length(apa_gene_only))),
                            APA = c(rep(TRUE, times = length(int_gene)),
                                    rep(FALSE, times = length(deg_gene_only)),
                                    rep(TRUE, times = length(apa_gene_only)))))
  
  df$Both <- ifelse(df$DEG == TRUE & df$APA == TRUE, TRUE, FALSE)
  
  write.table(df,paste0(apa_dir,paste0(comp[[i]], collapse = '_vs_'),'_DEG_APA.txt'),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
}
