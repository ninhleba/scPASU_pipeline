library(dplyr)
library(tidyr)
library(goldmine)
library(igraph)
library(parallel)
library(argparser)

parser<-arg_parser(name="3e_identify_fragmented_peaks.R",description="Identify fragmented peaks because of aligning exome-reads against genome")

parser<-add_argument(
  parser,
  arg='--transcripts',
  short = '-t',
  type="character",
  default=NA,
  help="Enter the path to the transcripts causing peak fragmentation. Format: path/to/transcript_file")

parser<-add_argument(
  parser,
  arg='--realign_peaks',
  short = '-p',
  type="character",
  default=NA,
  help="Enter the path to the realign peaks. Format: path/to/realign_peak_file")

parser<-add_argument(
  parser,
  arg='--peak_count_dir',
  short = '-d',
  type="character",
  default='./',
  help="Enter directory path to where the read counts and spliced read counts for merge candidates and realign peak counts are. Format: path/to/dir/ Default= ./")

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
  arg='--spliced_read_pct_thres',
  short = '-s',
  default=40,
  type='numeric',
  help="Cutoff threshold for spliced read count percentage to filter merge candidates. Default = 40")

parser<-add_argument(
  parser,
  arg='--realign_peak_read_pct_thres',
  short = '-r',
  default=40,
  type='numeric',
  help="Cutoff threshold for realign peak read count percentage of total spliced read count to filter merge candidates. Default = 40")

parser<-add_argument(
  parser,
  arg='--peaks',
  short = '-b',
  type="character",
  default=NA,
  help="Enter the path to the peak file from the last filtering step. Format: path/to/peaks_file")

parser<-add_argument(
  parser,
  arg='--ncores',
  short = '-n',
  default=4,
  type='numeric',
  help="Number of cores to use when processing in parallel. Default = 4")

args <- parse_args(parser)

transcripts_file <- args$transcripts
realign_peaks_file <- args$realign_peaks
peak_count_dir <- args$peak_count_dir
fprefix <- args$file_prefix
outdir <- args$outdir
ncores <- args$ncores %>% as.numeric
spliced_read_pct_thres <- args$spliced_read_pct_thres
realign_peak_read_pct_thres <- args$realign_peak_read_pct_thres
peaks_file<- args$peaks

### Peak count loading and processing ###
plus_realign_peak_count_file <- list.files(peak_count_dir, pattern = paste0(fprefix,'_spliced_reads_realign_peaks_plus_peak_count.rds'), full.names = TRUE)
minus_realign_peak_count_file <- list.files(peak_count_dir, pattern = paste0(fprefix,'_spliced_reads_realign_peaks_minus_peak_count.rds'), full.names = TRUE)
plus_realign_peak_count <- readRDS(plus_realign_peak_count_file)
minus_realign_peak_count <- readRDS(minus_realign_peak_count_file)
colnames(plus_realign_peak_count$counts) <-'realign_peak_count'
colnames(minus_realign_peak_count$counts) <- 'realign_peak_count'
plus_realign_peak_count <- data.frame(realign_peak = row.names(plus_realign_peak_count$counts),plus_realign_peak_count$counts)
minus_realign_peak_count <- data.frame(realign_peak = row.names(minus_realign_peak_count$counts),minus_realign_peak_count$counts)
realign_peak_count <- rbind(plus_realign_peak_count,minus_realign_peak_count)

plus_peak_count_file <- list.files(peak_count_dir, pattern = paste0(fprefix,'_plus_subset_peak_count.rds'), full.names = TRUE)
minus_peak_count_file <- list.files(peak_count_dir, pattern = paste0(fprefix,'_minus_subset_peak_count.rds'), full.names = TRUE)
plus_peak_count <- readRDS(plus_peak_count_file)
minus_peak_count <- readRDS(minus_peak_count_file)
colnames(plus_peak_count$counts) <-'peak_count'
colnames(minus_peak_count$counts) <- 'peak_count'
plus_peak_count <- data.frame(peak = row.names(plus_peak_count$counts),plus_peak_count$counts)
minus_peak_count <- data.frame(peak = row.names(minus_peak_count$counts),minus_peak_count$counts)
peak_count <- rbind(plus_peak_count,minus_peak_count)

plus_spliced_only_peak_count_file <- list.files(peak_count_dir, pattern = paste0(fprefix,'_plus_subset_spliced_reads_peak_count.rds'), full.names = TRUE)
minus_spliced_only_peak_count_file <- list.files(peak_count_dir, pattern = paste0(fprefix,'_minus_subset_spliced_reads_peak_count.rds'), full.names = TRUE)
plus_spliced_only_peak_count <- readRDS(plus_spliced_only_peak_count_file)
minus_spliced_only_peak_count  <- readRDS(minus_spliced_only_peak_count_file)
colnames(plus_spliced_only_peak_count$counts) <-'spliced_only_peak_count'
colnames(minus_spliced_only_peak_count$counts) <- 'spliced_only_peak_count'
plus_spliced_only_peak_count <- data.frame(peak = row.names(plus_spliced_only_peak_count$counts),plus_spliced_only_peak_count$counts)
minus_spliced_only_peak_count <- data.frame(peak = row.names(minus_spliced_only_peak_count$counts),minus_spliced_only_peak_count$counts)
spliced_only_peak_count <- rbind(plus_spliced_only_peak_count,minus_spliced_only_peak_count)

###
cat('Splice in transcripts of interest and renumber the genomic coordinates contiguously \n')
transcripts <- read.delim(transcripts_file)
transcripts <- split(transcripts, transcripts$isoform.id)
transcripts <- mclapply(transcripts, function(x) {
  x$width <- (x$end - x$start)+1
  new_end <- cumsum(x$width)
  x$end <- new_end
  x$start <- (x$end - x$width)+1
  return(x)
}, mc.cores = ncores)

transcripts <- do.call('rbind',transcripts)
colnames(transcripts)[9] <- 'chr'
transcripts <- makeGRanges(transcripts, strand = T)

cat('Decide what split candidates to merge based on realign peaks \n')
peaks <- fread(realign_peaks_file)
peaks <- makeGRanges(peaks, strand = T)

ovl <- findOverlaps(transcripts,peaks) %>% as.data.frame()

merge <- transcripts[ovl$queryHits] %>% as.data.frame()
merge$spliced_in_peak <- peaks[ovl$subjectHits]$peakID

merge <- merge[!sapply(merge$peak, function(x) nchar(x) == 0),]

merge <- merge %>% separate_rows(peak, sep = ',') %>% mutate(peak = trimws(peak))

merge$realign_peak_count <- realign_peak_count$realign_peak_count[match(merge$spliced_in_peak,realign_peak_count$realign_peak)]
merge$peak_count <- peak_count$peak_count[match(merge$peak,peak_count$peak)]
merge$spliced_only_peak_count <- spliced_only_peak_count$spliced_only_peak_count[match(merge$peak,spliced_only_peak_count$peak)]
merge$spliced_read_pct <- merge$spliced_only_peak_count / merge$peak_count * 100

realign_peaks_annot <- merge %>% distinct(spliced_in_peak, peak, spliced_only_peak_count) %>% group_by(spliced_in_peak)%>% 
  mutate(total_spliced_read_count = sum(spliced_only_peak_count))
merge$total_spliced_read_count <- realign_peaks_annot$total_spliced_read_count[match(merge$spliced_in_peak,realign_peaks_annot$spliced_in_peak)]
merge$realign_peak_read_pct <- merge$realign_peak_count / merge$total_spliced_read_count * 100

spliced_read_pct <- data.frame(spliced_read_pct = merge$spliced_read_pct[!duplicated(merge$peak)])
h1 <- ggplot(spliced_read_pct, aes(x = spliced_read_pct)) + 
  geom_histogram(binwidth = 5, color="black", fill="lightblue") + scale_x_continuous(breaks=seq(0,100,5)) +
  geom_vline(aes(xintercept=spliced_read_pct_thres),color="blue", linetype="dashed", linewidth=1) + 
  theme_classic() +
  labs(x="Spliced read percentage", y = "Count")

ggsave(paste0(outdir,fprefix,'_spliced_read_pct_hist.png'),width = 7, height = 7, units = 'in',h1)

realign_peak_read_pct <- data.frame(realign_peak_read_pct = merge$realign_peak_read_pct[!duplicated(merge$spliced_in_peak)])
realign_peak_read_pct$realign_peak_read_pct <- pmin(realign_peak_read_pct$realign_peak_read_pct, 100)
h2 <- ggplot(realign_peak_read_pct, aes(x = realign_peak_read_pct)) + 
  geom_histogram(binwidth = 5, color="black", fill="pink") + scale_x_continuous(breaks=seq(0,100,5)) +
  geom_vline(aes(xintercept=realign_peak_read_pct_thres),color="red", linetype="dashed", linewidth=1) + 
  theme_classic() +
  labs(x="Realign peak read percentage of total spliced read", y = "Count")
  
ggsave(paste0(outdir,fprefix,'_realign_peak_read_pct_hist.png'),width = 7, height = 7, units = 'in',h2)

cat('Only consider peak candidates with spliced read percent and corresponding realign peak read percent more than',spliced_read_pct_thres,'% and',realign_peak_read_pct_thres,'% respectively \n')

merge <- split(merge, merge$spliced_in_peak)

merge_filtered <- lapply(merge, function(x) {
  x <- subset(x, spliced_read_pct > spliced_read_pct_thres & realign_peak_read_pct > realign_peak_read_pct_thres)
  uniq_peaks <- x$peak %>% unique
  if (length(uniq_peaks) > 1) {
    df <- data.frame(peak = uniq_peaks, spliced_in_peak = unique(x$spliced_in_peak), 
                     tu = unique(x$tu), gene.id = unique(x$gene.id), name = unique(x$name),
                     realign_peak_read_pct = unique(x$realign_peak_read_pct))
    df$spliced_read_pct <- x$spliced_read_pct[match(df$peak,x$peak)]
    return(df)
    }
})

merge_filtered <- merge_filtered[!(lapply(merge_filtered, function(x) is.null(x)) %>% unlist)]
merge_filtered <- do.call('rbind', merge_filtered)

merge_bytu <- split(merge_filtered, merge_filtered$tu)
merge_adj <- lapply(merge_bytu, function(x) {
  peaks <- unique(x$peak)
  num_peak <- length(peaks)
  adj_matrix <- matrix(0, nrow = num_peak, ncol = num_peak,
                       dimnames = list(peaks, peaks))
  # Adjacency matrix
  for (i in 1:nrow(adj_matrix)) {
    for (j in 1:ncol(adj_matrix)) {
      int <- intersect(x$spliced_in_peak[x$peak == rownames(adj_matrix)[i]],
                       x$spliced_in_peak[x$peak == colnames(adj_matrix)[j]])
      if (length(int) > 0) {
        adj_matrix[i, j] <- 1
      }
    }
  }
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  clusters <- clusters(graph)
  group_membership <- membership(clusters)
  groups <- lapply(unique(group_membership), function(x) peaks[group_membership == x])
  # Remove group of one peak. This should not be necessary since we already remove group of single peaks earlier
  single <- lapply(groups, function(x) length(x)) %>% unlist == 1
  groups <- groups[!single]
  return(groups)
})

merge_adj <- lapply(merge_adj, function(x){
  names(x) <- paste0('mP',1:length(x))
  return(x)
})

merge_final <- matrix(ncol = 3, nrow = 0) %>% as.data.frame()
colnames(merge_final) <- c('TU','merged_peak','peak')

for (tu in names(merge_adj)) {
  for (mpeak in names(merge_adj[[tu]])) {
    df <- data.frame(TU = tu, merged_peak = paste0(tu,':',mpeak), peak = merge_adj[[tu]][[mpeak]])
    merge_final <- rbind(merge_final, df)
  }
}

cat('Identify',length(merge_final$peak),'peaks on',length(unique(merge_final$TU)),'TUs to merge \n')
cat('Annotate split peaks to merge in the final reference file \n')
ref <- fread(peaks_file)
ref$merged_peak <- merge_final$merged_peak[match(ref$final_annotation,merge_final$peak)]

write.table(ref,paste0(outdir,fprefix,'_final_peak_universe_updated.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

# Create SAF format peak ref too

# Select relevant columns
saf_ref<-ref %>% select(final_annotation,chr,start,end,strand)
colnames(saf_ref)<-c('GeneID','Chr','Start','End','Strand') 

cat('Creating SAF ref file \n')
write.table(saf_ref,paste0(outdir,fprefix,'_final_peak_universe_updated.saf'),sep='\t',quote=FALSE,row.names=FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
