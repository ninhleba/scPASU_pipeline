library(dplyr)
library(goldmine)
library(argparser)

parser<-arg_parser(name="3e_peak_calling_from_bedGragh.R",description="Call peak from alignment signals")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='scPASU',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--bg_dir',
  short = '-b',
  type="character",
  default='./',
  help="Enter the directory path of bedGraph files. Default=./")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--gap_threshold',
  short = '-g',
  type="numeric",
  default=0,
  help="Value under this will be considered gaps. Default: 0")

parser<-add_argument(
  parser,
  arg='--min_cov',
  short = '-m',
  default=10,
  type='numeric',
  help="Minimum coverage to retain peaks on each transcript. Default = 10")

args <- parse_args(parser)

### Arguments
fprefix <- args$file_prefix
bg_dir <- args$bg_dir
outdir <- args$outdir
gap_threshold <- args$gap_threshold %>% as.numeric
min_cov <- args$min_cov %>% as.numeric

cat('Load and process bedGraph files from both strands \n')
bg_m_file <- grep('.*minus.*bedGraph',list.files(bg_dir, full.names = TRUE),value=TRUE)
bg_p_file <- grep('.*plus.*bedGraph',list.files(bg_dir, full.names = TRUE),value=TRUE)
bg_m <- fread(bg_m_file)
bg_p <- fread(bg_p_file)

bg_m$strand <- '-'
bg_p$strand <- '+'
bg <- rbind(bg_p,bg_m)

colnames(bg) <- c('chr','start','end','score','strand')
transcripts <- bg$chr %>% unique

cat('Obtain peaks, i.e. contiguous regions with components larger than gap_threshold, from bedGraph files \n')
bg <- makeGRanges(bg, strand = T)
seqlevels(bg) <- transcripts
bg <- bg[order(bg)]

bg <- bg %>% as.data.frame()
bg <- split(bg, bg$seqnames)

bg_peak <- lapply(bg, function(x) {
  x$pos <- (x$score > gap_threshold)
  
  x <- x %>% group_by(peakID = rleid(pos)) %>%
    reframe(chr = unique(seqnames), start = min(start), end = max(end), width = max(end) - min(start),
            strand = unique(strand), sum_score = sum(score), sum_pos = all(pos)) %>% 
    filter(sum_pos == TRUE) %>% as.data.frame()
  
  if (nrow(x) == 1) {
    return(x)
  } else if (nrow(x) > 1) {
    perct <- x$sum_score / sum(x$sum_score) * 100
    x <- x[(perct > min_cov),]
    return(x)
    }
})

bg_peak <- do.call('rbind',bg_peak)
row.names(bg_peak) <- 1:nrow(bg_peak)
bg_peak$peakID <- ifelse(bg_peak$strand == '+',
                         paste0('spliced_reads_realign_plus_peak_',row.names(bg_peak)),
                         paste0('spliced_reads_realign_minus_peak_',row.names(bg_peak)))
bg_peak <- bg_peak[,c('chr','start','end','peakID','sum_score','strand')]

write.table(bg_peak, paste0(outdir,fprefix,'_spliced_reads_realign_peaks.txt'),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')