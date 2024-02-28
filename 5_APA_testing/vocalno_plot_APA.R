#https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

library(ggplot2)
library(dplyr)
library(ggrepel)

setwd('~/Desktop/SC296_scPASU_run/outputs/cancercell/5b_APA_testing/iAP_G0_v_iAP_G4/')
options(scipen = 999, digits = 1000)

files <- list.files()
for (file in grep('res.txt$', files, value = TRUE)) {
  assign(gsub('_res.txt','',file),read.delim(file))
}

volcano_plot <- function(res,foldchange.col,ident1,ident2,colors=c("blue", "red", "black"),max.overlaps = 10) {
  res <- res[,c('tu','peak','padj',foldchange.col,'int_sig')]
  colnames(res) <- c('tu','peak','padj','log2fc','int_sig')
  res$diffexpressed <- ifelse(res$int_sig==FALSE,"NO","YES")
  res$diffexpressed <- ifelse(res$int_sig==TRUE & res$log2fc > 0,
                              ident1,res$diffexpressed)
  res$diffexpressed <- ifelse(res$int_sig==TRUE & res$log2fc < 0,
                              ident2,res$diffexpressed)
  
  spl <- split(res,res$tu)
  opp.tu <- lapply(spl,function(x) {
    sig.peaks <- subset(x, diffexpressed != 'NO')
    if (length(unique(sig.peaks$diffexpressed)) > 1) {
      return(sig.peaks$peak)
    }
  })
  
  features.to.label <- opp.tu %>% unlist() %>% unname()
  res$delabel <- NA
  res$delabel <- features.to.label[match(res$peak,features.to.label)]
  
  names(colors) <- c(ident1, ident2, "NO")
  
  res$diffexpressed <- factor(res$diffexpressed, levels = c(ident1,ident2,'NO'))
  
  g <- ggplot(data=res, aes(x=log2fc, y=-log10(padj), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(max.overlaps = max.overlaps) +
    scale_color_manual(values=colors) +
    geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="red") +
    geom_hline(yintercept=-log10(0.01), col="red")
  
  return(g)
}


g <- volcano_plot(iAP_G0_v_iAP_G4,foldchange.col='l2_iAP_G0_over_iAP_G4_frac',
                     ident1='iAP_G0',ident2='iAP_G4',colors = c("red","blue","black"), max.overlaps = 40)

g

ggsave('iAP_G0_v_iAP_G4_volcano_plot.png',width = 16, height = 8, units = 'in', bg = 'white',g)


