#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)
dir=args[1]
strand=args[2]
chr=args[3]
min_polya=as.numeric(args[4])

prdir<-paste0(dir,'3b_pr_dir/')
peak_ref_dir<-paste0(dir,'3a_peak_ref_table/')
count_dir<-paste0(dir,'2_peaks_count_dir/')
out<-peak_ref_dir
out2<-paste0(dir,'4_polya_supported_peak_ref/')

pr_file<-paste0(prdir,'pr_',chr,'_sorted_',strand,'.bed')
peak_ref_file<-paste0(peak_ref_dir,'peak_ref_',chr,'_',strand,'.bed')
count_file<-paste0(count_dir,'peaks_',chr,'_count_sorted_',strand,'.txt')

pr<-fread(pr_file)
colnames(pr)<-c('peakID','start','end')

peak_ref<-fread(peak_ref_file)

colnames(peak_ref)<-c('peak_chr','peak_start','peak_end','peakID','peak_strand','pr_chr','pr_start','pr_end','readID','A_len','pr_strand','polya')

count<-fread(count_file)
colnames(count)<-c('count','peakID')

# Remove polya coordinate columns
peak_ref<-peak_ref[,-c(7,8,9)]

peak_ref$pr_start<-0
peak_ref$pr_end<-0

# Add inferred PR coordinates 
indx<-match(pr$peakID,peak_ref$peakID)
peak_ref$pr_start[indx]<-pr$start
peak_ref$pr_end[indx]<-pr$end

# Add peak and PR width
peak_ref$peak_width=peak_ref$peak_end-peak_ref$peak_start
peak_ref$pr_width=peak_ref$pr_end-peak_ref$pr_start

#Add polya count
peak_ref$polya_count<-0
indx2<-match(count$peakID,peak_ref$peakID)
peak_ref$polya_count[indx2]<-count$count

#save individual chr output
fname<-basename(peak_ref_file) %>% gsub('peak_ref','peak_ref_updated',.)

write.table(peak_ref,file = paste0(out,fname),sep = "\t",row.names = FALSE,quote=FALSE)

# Filter by min polyA criterion

keep<-which(peak_ref$polya_count>=min_polya)
peak_ref2<-peak_ref[keep,]

fname2<-gsub('.bed',paste0('_',min_polya,'polya.bed'),fname)
write.table(peak_ref2,paste0(out2,fname2),sep='\t',row.names=FALSE,quote=FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
