options(scipen = 999, digits = 20)

library(magrittr)
library(data.table)
library(dplyr)
library(goldmine)
library(gtools)
library(GenomicRanges)
library(parallel)
library(argparser)

parser<-arg_parser(name="2d_split_peaks.R",description="Split peaks before step 3")

parser<-add_argument(
  parser,
  arg='--prong',
  short = '-p',
  type="character",
  default='all_filtered_reads',
  help="Enter the name of the prong. Default=all_filtered_reads")

parser<-add_argument(
  parser,
  arg='--bg_dir',
  short = '-b',
  type="character",
  default='./',
  help="Enter the directory path of bedGraph files. Default=./")

parser<-add_argument(
  parser,
  arg='--summit_dir',
  short = '-s',
  type="character",
  default='./',
  help="Enter the path to the summits directory. Format: path/to/dir/ Default= ./")

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
  arg='--extn',
  short = '-e',
  type="character",
  default=NA,
  help="Enter direction and number of base pair to extend, separated by a comma. For example, 5,100")

parser<-add_argument(
  parser,
  arg='--gap_threshold',
  short = '-g',
  type="numeric",
  default=0,
  help="Value under this will be considered gaps. Default: 0")

parser<-add_argument(
  parser,
  arg='--chrs',
  short='-c',
  type='character',
  default=NA,
  help='Only split peaks on these chromosomes. Default: 24 canonical chromosomes of human genome')

parser<-add_argument(
  parser,
  arg='--ncores',
  short = '-n',
  default=4,
  type='numeric',
  help="Number of cores to use when processing in parallel. Default = 4")

args <- parse_args(parser)

prong <- args$prong
bg_dir <- args$bg_dir
summit_dir <- args$summit_dir
fprefix <- args$file_prefix
outdir <- args$outdir
gap_threshold <- args$gap_threshold %>% as.numeric
extn <- args$extn
chrs <- args$chrs
ncores <- args$ncores %>% as.numeric

if (is.na(chrs) == TRUE) {
  chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
           'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
           'chr20','chr21','chr22','chrX','chrY')
} else {
  chrs <- chrs %>% strsplit(split=',') %>% unlist
}

cat('Only split peaks on these chromosomes:',paste0(chrs,collapse=', '),'\n')

if(!dir.exists(outdir))
{dir.create(outdir, recursive = TRUE)}

### Process summit dir
cat('Load peaks with summits plus and minus strand separately \n')
summits_m <- fread(grep('.*minus.*peaks_summits.txt',list.files(summit_dir,full.names = TRUE),value=TRUE))
if (ncol(summits_m)!=10) {stop('Input summit file on minus strand should contain 10 columns from the original peaks.xls file \n')}
summits_p <- fread(grep('.*plus.*peaks_summits.txt',list.files(summit_dir,full.names = TRUE),value=TRUE))
if (ncol(summits_p)!=10) {stop('Input summit file on plus strand should contain 10 columns from the original peaks.xls file \n')}
summits_m$strand <- '-'
summits_p$strand <- '+'
summits <- rbind(summits_p,summits_m)
summits <- summits[,c(1:5,10,11)]
colnames(summits) <- c('chr','peak_start','peak_end','peak_width','summit','summitID','strand')
summits <- summits[which(summits$chr %in% chrs),]

if (is.na(extn) == FALSE) {
  extn <- strsplit(extn, split = ',') %>% unlist
  direction <- extn[1]
  ext_len <- extn[2] %>% as.numeric
  cat(paste('Extend peaks by',ext_len,'bp',direction,'prime \n'))
  if (direction == '3') {
    summits$peak_end <- ifelse(summits$strand=='+',summits$peak_end+ext_len,summits$peak_end)
    summits$peak_start <- ifelse(summits$strand=='-',summits$peak_start-ext_len,summits$peak_start)
  } else if (direction == '5') {
    summits$peak_start <- ifelse(summits$strand=='+',summits$peak_start-ext_len,summits$peak_start)
    summits$peak_end <- ifelse(summits$strand=='-',summits$peak_end+ext_len,summits$peak_end)
  } else {stop('Direction has to be either 3 or 5 \n')}
  summits$peak_width <- summits$peak_end - summits$peak_start
}

summits$peakID <- ifelse(str_detect(summits$summitID, "[a-z]$"), 
                         substring(summits$summitID,1,nchar(summits$summitID)-1), 
                         summits$summitID)
summits_freq <- summits %>% select(peakID) %>% group_by(peakID) %>% tally %>% as.data.frame
# Obtain valley ranges from summit file
cat('Obtain valley ranges, i.e. regions between two consecutive summits in a peak, from summit files\n')
summits_multimodal <- summits[which(summits$peakID %in% summits_freq[which(summits_freq$n > 1),'peakID']),]

summits_spl <- split(summits_multimodal,summits_multimodal$peakID)
valley <- mclapply(summits_spl, function(x) {
  s <- x[order(x$summit),]
  peakID <- unique(x$peakID)
  peak_start <- unique(x$peak_start)
  peak_end <- unique(x$peak_end)
  chr <- unique(x$chr)
  strand <- unique(x$strand)
  num_val <- nrow(s)-1
  v <- matrix(nrow = num_val, ncol = 6) %>% as.data.frame
  colnames(v) <- c('chr','start','end','valleyID','strand','peakID')
  v$chr <- chr
  v$start <- s$summit[1:(nrow(s)-1)]
  v$end <- s$summit[2:nrow(s)]
  v$valleyID <- paste0(peakID,'_',1:num_val)
  v$strand <- strand
  v$peakID <- peakID
  v$peak_start <- peak_start
  v$peak_end <- peak_end
  return(v)
}, mc.cores = ncores)

valley_df <- do.call('rbind',valley)

### Use bigwig/bedgraph to split peaks
cat('Load and process bedGraph files from both strands \n')
bg_m_file <- grep('.*minus.*bedGraph',list.files(bg_dir, full.names = TRUE),value=TRUE)
bg_p_file <- grep('.*plus.*bedGraph',list.files(bg_dir, full.names = TRUE),value=TRUE)
bg_m <- fread(bg_m_file)
bg_p <- fread(bg_p_file)

bg_m$strand <- '-'
bg_p$strand <- '+'
bg <- rbind(bg_p,bg_m)

colnames(bg) <- c('chr','start','end','score','strand')

cat('Obtain gap regions, i.e. contiguous regions where every nucleotide has a score equal or less than',gap_threshold,'from bedGraph files \n')
bg <- bg[which(bg$chr %in% chrs),]
bg <- makeGRanges(bg, strand = T)
seqlevels(bg) <- chrs
bg <- bg[order(bg)]

bg <- bg %>% as.data.frame()
bg <- split(bg, bg$seqnames)

bg_gap <- mclapply(bg, function(x) {
  chr <- unique(x$seqnames)
  strand <- unique(x$strand)
  
  x$pos <- (x$score > gap_threshold)
  
  x <- x %>% group_by(peakID = rleid(pos)) %>%
    reframe(chr = unique(chr), start = min(start), end = max(end), width = max(end) - min(start),
            strand = unique(strand), sum_score = sum(score), sum_pos = all(pos)) %>%
    filter(sum_pos == FALSE) %>% as.data.frame()
}, mc.cores = ncores)

bg_gap <- do.call('rbind',bg_gap)
row.names(bg_gap) <- 1:nrow(bg_gap)
bg_gap$peakID <- ifelse(bg_gap$strand == '+',
                        paste0(prong,'_gap_plus_',row.names(bg_gap)),
                        paste0(prong,'_gap_minus_',row.names(bg_gap)))
bg_gap <- bg_gap[,c('chr','start','end','peakID','sum_score','strand')]

# Find gaps in macs2 peaks
gaps <- makeGRanges(bg_gap, strand = T)
seqlevels(gaps) <- chrs
gaps <- gaps[order(gaps)]
valley_rg <- makeGRanges(valley_df, strand = T)
seqlevels(valley_rg) <- chrs
valley_rg <- valley_rg[order(valley_rg)]

cat('Obtain gap regions that lie totally within valley ranges \n')
ovl <- findOverlaps(gaps,valley_rg,type='within') %>% as.data.frame()
valley_rg <- as.data.frame(valley_rg)
gaps <- as.data.frame(gaps)

cat('Integrate peak and gap information to do peak splitting with these gaps \n')
ovl$peakID<- valley_rg[ovl$subjectHits,'peakID']
ovl$peak_chr <- valley_rg[ovl$subjectHits,'seqnames']
ovl$peak_start <- valley_rg[ovl$subjectHits,'peak_start']
ovl$peak_end <- valley_rg[ovl$subjectHits,'peak_end']
ovl$peak_strand <- valley_rg[ovl$subjectHits,'strand']

ovl$gapID <- gaps[ovl$queryHits,'peakID']
ovl$gap_chr <- gaps[ovl$queryHits,'seqnames']
ovl$gap_start <- gaps[ovl$queryHits,'start']
ovl$gap_end <- gaps[ovl$queryHits,'end']
ovl$gap_strand <- gaps[ovl$queryHits,'strand']

ovl_spl <- split(ovl, ovl$peakID)

split_peaks <- mclapply(ovl_spl, function(x) {
  peak_gaps <- x[order(x$gap_start),]
  peakID <- unique(x$peakID)
  chr <- unique(x$peak_chr)
  strand <- unique(x$peak_strand)
  starts <- c(unique(x$peak_start),x$gap_end)
  ends <- c(x$gap_start,unique(x$peak_end))
  stopifnot(length(starts)==length(ends))
  num_peak <- length(starts)
  split_peaks <- matrix(nrow = num_peak, ncol = 5) %>% as.data.frame
  colnames(split_peaks) <- c('chr','start','end','peakID','strand')
  split_peaks$chr <- chr
  split_peaks$start <- starts
  split_peaks$end <- ends
  split_peaks$peakID <- paste0(peakID,letters[1:num_peak])
  split_peaks$strand <- strand
  split_peaks$peak_width <- split_peaks$end - split_peaks$start
  return(split_peaks)
}, mc.cores = ncores)

split_peaks_df <- do.call('rbind',split_peaks)
split_peaks_df <- split_peaks_df[,c('chr','start','end','peakID','peak_width','strand')]
split_peaks_p <- split_peaks_df[which(split_peaks_df$strand == '+'),]
split_peaks_m <- split_peaks_df[which(split_peaks_df$strand == '-'),]

cat('Obtain peaks that are not split \n')
peaks_not_split <- summits[!(summits$peakID %in% unique(ovl$peakID)),]
peaks_not_split <- peaks_not_split[,c('chr','peak_start','peak_end','peakID','peak_width','strand')]
peaks_not_split <- peaks_not_split[!duplicated(peaks_not_split$peakID),]

colnames(peaks_not_split) <- c('chr','start','end','peakID','peak_width','strand')

cat('Merging split and not split peaks in an updated peak ref \n')
peaks_updated <- rbind(peaks_not_split,split_peaks_df)
cat('Split peak ref by strand \n')
peaks_updated_p <- peaks_updated[which(peaks_updated$strand == '+'),]
peaks_updated_m <- peaks_updated[which(peaks_updated$strand == '-'),]

write.table(bg_gap,
            paste0(outdir,fprefix,'_gaps.txt'),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
write.table(split_peaks_p,
            paste0(outdir,fprefix,'_plus_split.narrowPeak'),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
write.table(split_peaks_m,
            paste0(outdir,fprefix,'_minus_split.narrowPeak'),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
write.table(peaks_updated_p,
            paste0(outdir,fprefix,'_plus_peaks.narrowPeak'),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
write.table(peaks_updated_m,
            paste0(outdir,fprefix,'_minus_peaks.narrowPeak'),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
