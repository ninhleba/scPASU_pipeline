#!/usr/bin/env Rscript
options(scipen = 999, digits = 20)

library(goldmine)
library(rtracklayer)
library(dplyr)
library(data.table)
library(gtools)
library(argparser)

parser<-arg_parser(name="3a_assign_tu.R",description="Assign TU to peaks")

parser<-add_argument(
  parser,
  arg='--script_dir',
  short = '-d',
  type="character",
  default='./',
  help="Enter the path to dir of the scPASU_functions.R file. Format: path/to/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--peaksdir',
  short = '-p',
  type="character",
  default='./',
  help="Enter the path to the final peak directory from step 2. Format: path/to/peaks_dir/")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='scPASU',
  help="Enter file prefix to mark files for current run")

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
  arg='--peaks',
  short = '-b',
  type="character",
  default=NA,
  help="Enter the path to the merged polyA peak files from step 2 or the peak file from the last filtering step. Format: path/to/peaks_file")

parser<-add_argument(
  parser,
  arg='--turef',
  short = '-t',
  type="character",
  default=NA,
  help="Enter the path to the TU reference file. Format: path/to/tu_ref_file")

parser<-add_argument(
  parser,
  arg='--extn',
  short = '-e',
  type="character",
  default=NA,
  help="Enter direction and number of base pair to extend, separated by a comma. For example, 5,100")

parser<-add_argument(
  parser,
  arg='--chrs',
  short='-c',
  type='character',
  default=NA,
  help='Insert the chromosomes here to order the peaks by coordinate. Default: 24 canonical chromosomes of human genome')

args <- parse_args(parser)

script_dir <- args$script_dir
peaks_dir <- args$peaksdir
fprefix <- args$file_prefix
gtf_file <- args$gtf
outdir <- args$outdir
genes_file<- args$turef
peaks_file<- args$peaks
extn <- args$e
chrs <- args$chrs

output_file<-paste0(outdir,"/",fprefix,"_jtu.rds")

if (is.na(chrs) == TRUE) {
  chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
           'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
           'chr20','chr21','chr22','chrX','chrY')
} else {
  chrs <- chrs %>% strsplit(split=',') %>% unlist
}

source(paste0(script_dir,'scPASU_functions.R'))

# Create dir if doesn't exist

if(!dir.exists(outdir))
        {dir.create(outdir, recursive = TRUE)}

### Create genes table ###
if(is.na(genes_file) == FALSE)
  {
  cat(paste('Loading red_ens from',genes_file,'\n'))
  red_ens<-readRDS(genes_file)
  }else{
    g<-rtracklayer::import(gtf_file)
    gtf=as.data.frame(g)
    genes <- reformat_gtf(gtf,chrs=chrs)
    red_ens <- create_TU_from_gtf(genes,flank=5000,outdir=outdir,save=TRUE,
                                  final_tu_annot='gene_symbol',chrs=chrs)
}

## Peaks file ##
if(is.na(peaks_file) == FALSE)
  {
  cat(paste('Loading peaks file from',peaks_file,'\n'))
  merged<-fread(peaks_file)
  if(!('polya' %in% colnames(merged)))
    {
    merged$polya<-rep('polya',nrow(merged))}

  }else{
    
  files_m<-grep('.*minus.*bed',list.files(paste0(peaks_dir,'minus/4_polya_supported_peak_ref/'), full.names = TRUE),value=TRUE) %>% as.list
  file_lst_m<-lapply(files_m,fread,header=TRUE)
  merged_m<-do.call(bind_rows,file_lst_m)
  
  files_p<-grep('.*plus.*bed',list.files(paste0(peaks_dir,'plus/4_polya_supported_peak_ref/'), full.names = TRUE),value=TRUE) %>% as.list
  file_lst_p<-lapply(files_p,fread,header=TRUE)
  merged_p<-do.call(bind_rows,file_lst_p)
  
  merged<-rbind(merged_m,merged_p)
  
  # Rename columns since GRange recognizes start, end and strand only
  cat('Renaming colnames to create GRange obj \n')
  colnames(merged)<-gsub('\\bpeak_chr\\b','chr',colnames(merged))
  colnames(merged)<-gsub('\\bpeak_start\\b','start',colnames(merged))
  colnames(merged)<-gsub('\\bpeak_end\\b','end',colnames(merged))
  colnames(merged)<-gsub('\\bpeak_strand\\b','strand',colnames(merged))
  
  
  # Save the dataframe
  cat('Saving table \n')
  write.table(merged,paste0(outdir,fprefix,'_peak_universe.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

}

#convert to GRange

cat('Creating GRange obj \n')
peaks<-makeGRanges(merged,strand=T)

cat('Assign TU to peaks \n')

jtu <- joinTus_peaks(allpeaks=peaks,rg=red_ens)

cat('Save TU-assigned peaks, all peaks and TU reference \n')
save(jtu,peaks,red_ens,file=output_file)
load(output_file)

# Create Peak reference with all relevant columns merged after TU assignment #
cat('Add other relevant cols \n')

jtu$join<-as.data.frame(jtu$join)
# Remove multi-TU peaks
cat('Remove peaks that are assigned to multiple TUs, i.e. non-unique_peak TUs \n')
r1<-which(jtu$join$unique_peak==FALSE)

if(length(r1)!=0) {
 tu_tab<-jtu$join[-r1,] 
 }else{
 tu_tab<-jtu$join
 }

# Unique TU (TU contains only peaks assigned uniquely to TU) status is determined before removing peaks with unique_peak==FALSE so there can still be FALSE unique_tu peaks left

# Create peak per TU count
tu_peak<-jtu$join %>% dplyr::select(peak,tu) %>% group_by(tu) %>% tally()

# Peaks (saved from merging all peaks from MACS2 output)
peaks<-as.data.frame(jtu$polya_peaks)
found<-match(tu_tab$peak,peaks$peakID)
cat(paste(length(found),'peaks overlap uniquely with the reference TU. Only move forward with these peaks \n'))
matched_peaks<-peaks[found,]

# Sort all tables
# tu_tab<-tu_tab[mixedorder(tu_tab$peak),] Peak numbering, in this new approach, is not necessarily reflective of their relative position
# matched_peaks<-matched_peaks[mixedorder(matched_peaks$peakID),]

# Now lets transfer peak info to TU table
stopifnot(identical(tu_tab$peak,matched_peaks$peakID)) # TRUE

tu_tab$chr<-matched_peaks$seqnames
tu_tab$start<-matched_peaks$start
tu_tab$end<-matched_peaks$end
tu_tab$strand<-matched_peaks$strand
tu_tab$gene<-strsplit(tu_tab$tu_anno,split=':') %>% sapply(.,'[[',2)

# Also transfer other info for future use. This table can then serve as comprehensive features table

tu_tab$peak_width<-matched_peaks$peak_width
tu_tab$pr_chr<-matched_peaks$pr_chr
tu_tab$pr_start<-matched_peaks$pr_start
tu_tab$pr_end<-matched_peaks$pr_end
tu_tab$pr_strand<-matched_peaks$pr_strand
tu_tab$pr_width<-matched_peaks$pr_width
tu_tab$polya_count<-matched_peaks$polya_count
tu_tab$peakID<-matched_peaks$peakID

# Now add peakID to tu annotation column (multi peak TU will look like TU1:gene:P1, TU1:gene:P2 and TU1:gene:P3 while single peak TU will look like TU2:gene:P0)
cat('Sorting peaks by genomic coordinates before peak numbering \n')
tu_tab <- makeGRanges(tu_tab, strand = T)
seqlevels(tu_tab) <- chrs
tu_tab <- tu_tab[order(tu_tab)]
tu_tab <- as.data.frame(tu_tab)
tu_tab$width <- NULL
colnames(tu_tab)[1] <- 'chr'

# Get TU count
tu_count<-tu_tab %>% group_by(tu) %>% tally()

# Sort
tu_count<-tu_count[order(tu_count$n),]

#Single peak per TU
single<-which(tu_count$n==1)
s<-match(tu_count$tu[single],tu_tab$tu)
single_tu<-tu_tab[s,]
cat(paste(single_tu$tu %>% unique() %>% length(),'TUs have only one peaks. These peaks are hence annotated P0 \n'))
single_tu$final_annotation<-paste0(single_tu$tu_anno,':P0')

# Multi peak TU -
multi_tu<-tu_tab[-s,]

# Split by strand
multi_tu_p<-multi_tu[multi_tu$strand=='+',]
multi_tu_m<-multi_tu[multi_tu$strand=='-',]

# Now add final annotation col
multi_tu_p2<-create_final_annotation_col(multi_tu_p)
multi_tu_m2<-create_final_annotation_col(multi_tu_m,is_minus=TRUE)

# Merge plus and minus strand
multi_tu2<-rbind(multi_tu_p2,multi_tu_m2)

# Add P0 genes back to ref
merged_tu3<-rbind(single_tu,multi_tu2)

# Save
merged_tu3<-merged_tu3 %>% as.data.frame()

# Assess how many bp should be extended upstream
tu_ref <- red_ens$tu %>% makeDT
tu_spl <- split(merged_tu3, merged_tu3$tu)
cat('Calculate distances between peaks within each TU\n')
tu_btwn_peak_distance <- lapply(tu_spl, function(x) {
  num_peak <- nrow(x)
  if (num_peak > 1) {
    TU_id <- x$tu %>% unique
    TU_coord <- tu_ref[which(tu_ref$tu == TU_id),]
    strand <- x$strand %>% as.character %>% unique
    if (strand == '+') {
      start_increasing <- sort(x$start, decreasing = FALSE)
      if (identical(x$start,start_increasing) == FALSE) {
        cat(paste('Peak starts within TU',TU_id,', which is on + strand, are not sorted increasingly. Sorting... \n'))
        sorted_order <- match(start_increasing,x$start)
        x <- x[sorted_order,]
      }
      tu_start <- TU_coord$start
      btwn_dist <- c()
      #btwn_dist <- c(btwn_dist,x[1,'start'] - tu_start)
      for (i in 1:(nrow(x)-1)) {
        d <- x[i+1,'start'] - x[i,'end']
        btwn_dist <- c(btwn_dist,d)
      }
    } else {
      end_decreasing <- sort(x$end, decreasing = TRUE)
      if (identical(x$end, end_decreasing) == FALSE) {
        cat(paste('Peak ends within TU',TU_id,', which is on - strand, are not sorted decreasingly. Sorting... \n'))
        sorted_order <- match(end_decreasing, x$end)
        x <- x[sorted_order,]
      }
      tu_end <- TU_coord$end
      btwn_dist <- c()
      #btwn_dist <- c(btwn_dist,tu_end - x[1,'end'])
      for (i in 1:(nrow(x)-1)) {
        d <- x[i,'start'] - x[i+1,'end']
        btwn_dist <- c(btwn_dist,d)
      }
    }
    min_dist <- min(btwn_dist)
    max_dist <- max(btwn_dist)
    mean_dist <- mean(btwn_dist)
    tu_btwn_dist <- data.frame(TU_id, strand, min_dist, max_dist, mean_dist, num_peak)
    return(tu_btwn_dist)
  }
})
tu_btwn_peak_distance[sapply(tu_btwn_peak_distance, is.null)] <- NULL
tu_btwn_peak_distance <- do.call('rbind',tu_btwn_peak_distance)
tu_btwn_peak_distance <- tu_btwn_peak_distance[order(tu_btwn_peak_distance$min_dist),]
write.table(tu_btwn_peak_distance,paste0(outdir,fprefix,'_tu_btwn_peak_distance.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

cat('Calculate distances from peaks to their TU start site\n')
tu_tss_peak_distance <- lapply(tu_spl, function(x) {
  num_peak <- nrow(x)
  TU_id <- x$tu %>% unique
  TU_coord <- tu_ref[which(tu_ref$tu == TU_id),]
  strand <- x$strand %>% as.character %>% unique
  if (strand == '+') {
    start_increasing <- sort(x$start, decreasing = FALSE)
    if (identical(x$start,start_increasing) == FALSE) {
      cat(paste('Peak starts within TU',TU_id,', which is on + strand, are not sorted increasingly. Sorting... \n'))
      sorted_order <- match(start_increasing,x$start)
      x <- x[sorted_order,]
    }
    tu_start <- TU_coord$start
    tss_dist <- c()
    for (i in 1:(nrow(x))) {
      d <- x[i,'start'] - tu_start
      tss_dist <- c(tss_dist,d)
    }
  } else {
    end_decreasing <- sort(x$end, decreasing = TRUE)
    if (identical(x$end, end_decreasing) == FALSE) {
      cat(paste('Peak ends within TU',TU_id,', which is on - strand, are not sorted decreasingly. Sorting... \n'))
      sorted_order <- match(end_decreasing, x$end)
      x <- x[sorted_order,]
    }
    tu_end <- TU_coord$end
    tss_dist <- c()
    for (i in 1:(nrow(x))) {
      d <- tu_end - x[i,'end']
      tss_dist <- c(tss_dist,d)
    }
  }
  min_dist <- min(tss_dist)
  max_dist <- max(tss_dist)
  mean_dist <- mean(tss_dist)
  tss_dist <- data.frame(TU_id, strand, min_dist, max_dist, mean_dist, num_peak)
  return(tss_dist)
})
tu_tss_peak_distance <- do.call('rbind',tu_tss_peak_distance)
tu_tss_peak_distance <- tu_tss_peak_distance[order(tu_tss_peak_distance$min_dist,decreasing = FALSE),]
write.table(tu_tss_peak_distance,paste0(outdir,fprefix,'_tu_tss_peak_distance.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

# Extend peaks upstream
if (is.na(extn) == FALSE) {
  extn <- strsplit(extn, split = ',') %>% unlist
  direction <- extn[1]
  if ((direction %in% c('3','5')) == FALSE) {stop('Direction has to be either 3 or 5 \n')}
  ext_len <- extn[2] %>% as.numeric
  if (is.na(ext_len) == TRUE) {stop('Extension length has to be a number \n')}
  cat(paste('Extend peaks by',ext_len,'bp',direction,'prime \n'))
  merged_tu3 <- peak_ext(merged_tu3, ext_len=100, direction='5')
}

write.table(merged_tu3,paste0(outdir,fprefix,'_peak_universe_updated.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

# Modify flank regions of the genes.rds file according to the peaks

red_ens_flank_updated <- update_flank(merged_tu3,red_ens,chrs=chrs)

# Create SAF format peak ref too

# Select relevant columns
cols<-c('final_annotation','chr','start','end','strand')
select<-match(cols,colnames(merged_tu3))
saf_ref<-merged_tu3[,select]
colnames(saf_ref)<-c('GeneID','Chr','Start','End','Strand') 

cat('Creating SAF ref file \n')
write.table(saf_ref,paste0(outdir,fprefix,'_peak_universe_updated.saf'),sep='\t',quote=FALSE,row.names=FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
