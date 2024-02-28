library(magrittr)
library(data.table)
library(dplyr)
library(goldmine)
library(gtools)
library(GenomicRanges)
library(parallel)
library(argparser)

parser<-arg_parser(name="3c_merge_two_prongs.R",description="Merge peaks from two prongs")

parser<-add_argument(
  parser,
  arg='--indir',
  short = '-i',
  type="character",
  default='./',
  help="Enter the path to the output directory where the last peak references were created. Format: path/to/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='scPASU',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--chrs',
  short='-c',
  type='character',
  default=NA,
  help='Insert the chromosomes here to order the peaks by coordinate. Default: 24 canonical chromosomes of human genome')

parser<-add_argument(
  parser,
  arg='--ncores',
  short = '-n',
  default=4,
  type='numeric',
  help="Number of cores to use when processing in parallel. Default = 4")

args <- parse_args(parser)

indir <- args$indir
fprefix <- args$file_prefix
outdir <- args$outdir
chrs <- args$chrs
ncores <- args$ncores %>% as.numeric

if (is.na(chrs) == TRUE) {
  chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
           'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
           'chr20','chr21','chr22','chrX','chrY')
} else {
  chrs <- chrs %>% strsplit(split=',') %>% unlist
}

if(!dir.exists(outdir))
{dir.create(outdir, recursive = TRUE)}

# Load peaks from both prongs

all_filtered_reads_peaks <- fread(paste0(indir,'all_filtered_reads/',fprefix,'_all_filtered_reads_3rdtu_assigned_peak_universe_updated.txt'))
polyA_reads_peaks <- fread(paste0(indir,'polyA_reads/',fprefix,'_polyA_reads_3rdtu_assigned_peak_universe_updated.txt'))

all_filtered_reads_peaks$approach <- 'all_filtered_reads'
polyA_reads_peaks$approach <- 'polyA_reads'

# Count the number of unique TUs to each prong as well as intersect TU
all_filtered_reads_peaks_tu_freq <- all_filtered_reads_peaks %>% dplyr::select(peak,tu) %>% group_by(tu) %>% tally()
colnames(all_filtered_reads_peaks_tu_freq)[2] <- 'all_filtered_reads'
polyA_reads_peaks_tu_freq <- polyA_reads_peaks %>% dplyr::select(peak,tu) %>% group_by(tu) %>% tally()
colnames(polyA_reads_peaks_tu_freq)[2] <- 'polyA_reads'
union_tus <- union(all_filtered_reads_peaks_tu_freq$tu,polyA_reads_peaks_tu_freq$tu)
cat('There are total',length(union_tus),'TUs between two prongs \n')
intersect_tus <- intersect(all_filtered_reads_peaks_tu_freq$tu,polyA_reads_peaks_tu_freq$tu)
cat('There are',length(intersect_tus),'TUs shared between two prongs \n')

stopifnot(identical(colnames(all_filtered_reads_peaks),colnames(polyA_reads_peaks))) # TRUE

# Identify TUs unique to each prong. We don't have to merge these peaks for obvious reasons
all_filtered_reads_peaks_uniq_tu <- all_filtered_reads_peaks[-which(all_filtered_reads_peaks$tu %in% intersect_tus),]
cat('There are',length(unique(all_filtered_reads_peaks_uniq_tu$tu)),'TUs unique to all filtered reads prong, equal',nrow(all_filtered_reads_peaks_uniq_tu),'peaks\n')
polyA_reads_peaks_uniq_tu <- polyA_reads_peaks[-which(polyA_reads_peaks$tu %in% intersect_tus),]
cat('There are',length(unique(polyA_reads_peaks_uniq_tu$tu)),'TUs unique to polyA reads prong, equal',nrow(polyA_reads_peaks_uniq_tu),'peaks\n')

# Merge peaks from both prongs that share the same TUs
both_prongs_peaks <- rbind(all_filtered_reads_peaks,polyA_reads_peaks)
both_prongs_peaks <- both_prongs_peaks[which(both_prongs_peaks$tu %in% intersect_tus),]
both_prongs_peaks_spl <- split(both_prongs_peaks,both_prongs_peaks$tu)

merged_peaks <- mclapply(both_prongs_peaks_spl,function(x){
  tu_id <- x$tu %>% unique
  cols <- c(colnames(x),'cluster_id','merge')
  # Sanity check. These clusters should have been separated out
  if (length(unique(x$approach)) == 1) {
    cluster_id <- paste0(tu_id,'_C0')
    x$cluster_id <- cluster_id
    x$merge <- FALSE
    return(x)
  } else {
    peaks <- makeGRanges(x,strand=T)
    union_ranges <- reduce(peaks,with.revmap=FALSE)
    ovl <- findOverlaps(union_ranges,peaks) %>% as.data.frame()

    final_peaks <- matrix(nrow = 0, ncol = ncol(x)+2) %>% as.data.frame()
    colnames(final_peaks) <- cols
    
    for (q in unique(ovl$queryHits)) {
      cluster_id <- paste0(tu_id,'_C',q)
      sub_hits <- ovl[ovl$queryHits == q,'subjectHits'] %>% unique
      # Subset out peaks that lie within this cluster
      peaks_sub_hits <- peaks[sub_hits]
      peaks_sub_hits <- peaks_sub_hits[order(peaks_sub_hits)]
      leftmost_peak <- peaks_sub_hits[1] %>% as.data.frame()
      rightmost_peak <- peaks_sub_hits[length(peaks_sub_hits)] %>% as.data.frame()
      
      if (length(unique(peaks_sub_hits$approach)) == 1) {
        # For cluster comprising of peaks from just one prong, no need to merge or extend
        merged_peaks <- peaks_sub_hits %>% makeDT
        merged_peaks$merge <- FALSE
      } else {
        all_reads_idx <- which(peaks_sub_hits$approach == 'all_filtered_reads')
        polyA_reads_idx <- which(peaks_sub_hits$approach == 'polyA_reads')
        
        # Count peaks of each prong and pick the prong with the most number of peaks
        peak_count_by_prong <- c(all_filtered_reads = length(all_reads_idx), polyA_reads = length(polyA_reads_idx))
        
        chosen_prong <- which(peak_count_by_prong == max(peak_count_by_prong)) %>% names()
        if (length(chosen_prong) > 1) {
          # break tie by choosing the biggest peak width mean and if still does not break tie just go with all filtered read peaks
          peaks_sub_hits_all_filtered_reads <- peaks_sub_hits[all_reads_idx]
          peaks_sub_hits_polyA_reads <- peaks_sub_hits[polyA_reads_idx]
          mean_peak_width_by_prong <- c(all_filtered_reads = mean(peaks_sub_hits_all_filtered_reads$peak_width), 
                                        polyA_reads = mean(peaks_sub_hits_polyA_reads$peak_width))
          chosen_prong <- which(mean_peak_width_by_prong == max(mean_peak_width_by_prong)) %>% names()
          if (length(chosen_prong) > 1) {chosen_prong <- 'all_filtered_reads'}
        }
        
        merged_peaks <- peaks_sub_hits[which(peaks_sub_hits$approach == chosen_prong)]
        merged_peaks <- merged_peaks[order(merged_peaks)] %>% makeDT
        
        # Extend peaks at the edges of each cluster to cover the entire cluster range if they do not already
        # If extend peaks, also extend pr range if possible
        cluster_range <- union_ranges[q] %>% makeDT
        
        merged_peaks[1,'start'] <- cluster_range$start
        merged_peaks[nrow(merged_peaks),'end'] <- cluster_range$end
        
        if (leftmost_peak$pr_start < merged_peaks[1,'pr_start']) {merged_peaks[1,'pr_start'] = leftmost_peak$pr_start}
        if (leftmost_peak$pr_end > merged_peaks[1,'pr_end'] & leftmost_peak$pr_end < merged_peaks[1,'end']) {
          merged_peaks[1,'pr_end'] = leftmost_peak$pr_end}
        merged_peaks[1,'pr_width'] <- merged_peaks[1,'pr_end'] - merged_peaks[1,'pr_start']
        
        if (rightmost_peak$pr_end > merged_peaks[nrow(merged_peaks),'pr_end']) {merged_peaks[nrow(merged_peaks),'pr_end'] = rightmost_peak$pr_end}
        if (rightmost_peak$pr_start < merged_peaks[nrow(merged_peaks),'pr_start'] & rightmost_peak$pr_start > merged_peaks[nrow(merged_peaks),'start']) {
          merged_peaks[nrow(merged_peaks),'pr_start'] = rightmost_peak$pr_start}
        merged_peaks[nrow(merged_peaks),'pr_width'] <- merged_peaks[nrow(merged_peaks),'pr_end'] - merged_peaks[nrow(merged_peaks),'pr_start']
        
        merged_peaks$merge <- TRUE
      }

      merged_peaks <- as.data.frame(merged_peaks)
      colnames(merged_peaks)[1] <- 'chr'
      merged_peaks$cluster_id <- cluster_id
      merged_peaks$peak_width <- merged_peaks$end - merged_peaks$start
      merged_peaks$width <- NULL
      
      reorder_col_idx <- match(cols,colnames(merged_peaks))
      merged_peaks <- merged_peaks[,reorder_col_idx]
      final_peaks <- rbind(final_peaks,merged_peaks)
    }
    return(final_peaks)
    }
  }, mc.cores = ncores)

saveRDS(merged_peaks,paste0(outdir,fprefix,'_both_prongs_shared_tu.rds'))
merged_peaks_df <- do.call('rbind',merged_peaks) %>% select(-c('cluster_id'))
peaks_updated <- rbind(all_filtered_reads_peaks_uniq_tu,polyA_reads_peaks_uniq_tu)
peaks_updated$merge <- FALSE
stopifnot(identical(colnames(peaks_updated),colnames(merged_peaks_df)))
cat('Combine peaks from TUs unique to each prong with peaks from shared TUs\n')
peaks_updated <- rbind(peaks_updated,merged_peaks_df)
cat('Sort final peaks by chromosome and peak coordinate \n')
cols <- colnames(peaks_updated)
peaks_updated <- makeGRanges(peaks_updated,strand = T)
seqlevels(peaks_updated) <- chrs
peaks_updated <- peaks_updated[order(peaks_updated)]
peaks_updated <- peaks_updated %>% as.data.frame
colnames(peaks_updated)[1] <- 'chr'
peaks_updated <- peaks_updated %>% select(-c('width'))
peaks_updated <- peaks_updated[,match(cols,colnames(peaks_updated))]

write.table(peaks_updated,paste0(outdir,fprefix,'_both_prongs_merged_peak_universe_updated.txt'),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
