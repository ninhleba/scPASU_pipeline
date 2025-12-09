library(dplyr)
library(stringr)
library(goldmine)
library(parallel)
library(argparser)

parser<-arg_parser(name="phastCons30way_perbin_computing.R",description="Calculate phastCons30way score for each bin")

parser<-add_argument(
  parser,
  arg='--phastcon30way_dir',
  short = '-p',
  type="character",
  help="phastcon30way directory")

parser<-add_argument(
  parser,
  arg='--working_dir',
  short = '-d',
  type="character",
  default='./',
  help="Working directory, where the bed ranges for meme motif discovery inputs are stored. Default: ./")

parser<-add_argument(
  parser,
  arg='--meme_dir',
  short = '-o',
  type="character",
  help="Meme output directory.")

parser<-add_argument(
  parser,
  arg='--fprefix',
  short = '-f',
  type="character",
  default='WUI_APA_genes_UTR3',
  help="Enter the file prefix of meme outputs and bed ranges files. Default: WUI_APA_genes_UTR3")

parser<-add_argument(
  parser,
  arg='--bin_size',
  short = '-b',
  type="numeric",
  default=10,
  help="bin size. Default: 10")

parser<-add_argument(
  parser,
  arg='--num_motif',
  short = '-m',
  type="numeric",
  default=5,
  help="Number of motifs. Default: 5")

parser<-add_argument(
  parser,
  arg='--ncores',
  short = '-n',
  type="numeric",
  default=32,
  help="Number of cores for parallel computing. Default: 32")

args <- parse_args(parser)

phastcon30way_dir <- args$phastcon30way_dir
wd <- args$working_dir
meme_dir <- args$meme_dir
fprefix <- args$fprefix
bin_size = args$bin_size %>% as.numeric
num_motif = args$num_motif %>% as.numeric
ncores = args$ncores %>% as.numeric

utr3_ranges <- read.delim(paste0(wd,fprefix,'.bed'),
                          header = F)
colnames(utr3_ranges) <- c('chr','start','end','seq_name','score','strand')

for (memename in paste0('MEME-',1:num_motif)) {
  if (memename == 'MEME-1') {
    meme <- read.delim(paste0(meme_dir,fprefix,'_',memename,'.txt'), header = F)
    meme$meme_name <- memename
  } else {
    df <- read.delim(paste0(meme_dir,fprefix,'_',memename,'.txt'), header = F)
    df$meme_name <- memename
    meme <- rbind(meme,df)
    rm(df)
  }
}

colnames(meme) <- c('seq_name','site_strand','site_start','pval','flank1','site','flank2','meme_name')
meme$strand <- str_extract_all(meme$seq_name, "\\(([^)]+)\\)") %>% unlist
meme$strand <- gsub("\\(|\\)","",meme$strand)
meme$cis <- meme$site_strand==meme$strand
meme$seq_name <- gsub("\\(([^)]+)\\)","",meme$seq_name)
meme$chr <- utr3_ranges$chr[match(meme$seq_name,utr3_ranges$seq_name)]
meme$range_start <- utr3_ranges$start[match(meme$seq_name,utr3_ranges$seq_name)]
meme$range_end <- utr3_ranges$end[match(meme$seq_name,utr3_ranges$seq_name)]
meme$site_width <- nchar(meme$site)
meme$start <- ifelse(meme$strand=='+',meme$range_start+meme$site_start-1,meme$range_end-meme$site_start+1-meme$site_width+1)
meme$end <- meme$start+meme$site_width-1

utr3_ranges_wmemes <- subset(utr3_ranges, seq_name %in% unique(meme$seq_name))
utr3_gr <- utr3_ranges_wmemes %>% select(chr,start,end,seq_name)
utr3_gr$start %>% class()

utr3_gr <- apply(utr3_gr,1,function(x) {
  start <- seq(x[2], x[3], by = bin_size)
  chr <- x[1]
  df <- data.frame(chr = rep(chr,times=length(start)),
                   start = start)
  df$end <- df$start+bin_size-1
  df$end[nrow(df)] <- x[3]
  df$seq_name <- x[4]
  return(df)
})

utr3_gr <- do.call('rbind',utr3_gr)
utr3_gr <- utr3_gr %>% mutate(start = as.numeric(start), end = as.numeric(end))
utr3_gr <- makeGRanges(utr3_gr,strand = F)

motif_gr <- meme %>% select(chr,start,end,meme_name)
motif_gr <- makeGRanges(motif_gr,strand = F)

ovl <- findOverlaps(utr3_gr,motif_gr) %>% as.data.frame()
ovl$meme <- motif_gr$X[ovl$subjectHits]
ovl <- as.data.table(ovl)
ovl <- ovl[,.(meme = paste0(meme,collapse = ',')),by=queryHits]
utr3_gr <- as.data.frame(utr3_gr)
utr3_gr$meme <- 'no_motif'
utr3_gr$meme[ovl$queryHits] <- ovl$meme
colnames(utr3_gr)[1] <- 'chr'
utr3_gr$chr <- as.character(utr3_gr$chr)

utr3_spl <- split(utr3_gr,utr3_gr$chr)

utr3_phastcon30way_perbin <- lapply(utr3_spl,function(x){
  chr <- unique(x$chr)
  cat('Processing bins on',chr,'\n')
  
  p <- read.delim(paste0(phastcon30way_dir,chr,'.phastCons30way.wigFix'), sep = ' ', header = F)
  
  fixedStep_rowidx <- grep('fixedStep',p$V1)
  range_idx <- data.frame(start = fixedStep_rowidx, end = c(fixedStep_rowidx[-1]-1,nrow(p)))
  
  cat('Sanity checking phastcon30way file... \n')
  apply(range_idx,1,function(x) {
    df <- p[x[1]:x[2],]
    num_fs <- grep('\\bfixedStep\\b',df$V1) %>% length()
    vals_col2 <- df$V2 %>% unique()
    vals_col3 <- df$V3 %>% unique()
    vals_col4 <- df$V4 %>% unique()
    stopifnot(num_fs==1)
    stopifnot(length(vals_col2)==2)
    stopifnot(length(vals_col3)==2)
    stopifnot(length(vals_col4)==2)
    stopifnot("" %in% vals_col2)
    stopifnot("" %in% vals_col3)
    stopifnot("" %in% vals_col4)
    stopifnot(any(grepl(paste0('\\bchrom=',chr,'\\b'),vals_col2)))
    stopifnot(any(grepl('^start=',vals_col3)))
    stopifnot(any(grepl('\\bstep=1\\b',vals_col4)))
  })
  
  cat('Convert phastcon30way into genomic coordinates... \n')
  phastcon30way_ranges <- apply(range_idx,1,function(x) {
    start <- gsub('^start=','',p$V3[x[1]]) %>% as.numeric()
    subset_p <- p[x[1]:x[2],]
    num_positions <- nrow(subset_p)-1
    df <- data.frame(chr = chr, start = start:(start+num_positions-1))
    df$end <- df$start+1
    df$phastCons30 <- subset_p$V1[2:nrow(subset_p)] %>% as.numeric()
    return(df)
  })
  
  phastcon30way_ranges <- do.call('rbind',phastcon30way_ranges)
  # phastcon30way is 1-based so make it 0-based to be compatible with bedtool getfasta and ucsc genome browser
  phastcon30way_ranges$start <- phastcon30way_ranges$start - 1
  phastcon30way_ranges$end <- phastcon30way_ranges$end - 1
  
  cat('Assigning phastcon30way score per bin... \n')
  x <- x[order(as.numeric(x$start)),]
  row.names(x) <- 1:nrow(x)
  x <- split(x,row.names(x))
  x <- mclapply(x,function(x) {
    stopifnot(nrow(x)==1)
    start <- x$start
    end <- x$end
    scores <- phastcon30way_ranges$phastCons30[which(phastcon30way_ranges$start %in% start:end)]
    sum_score <- sum(scores)
    avg_score <- sum_score/length(scores)
    x$phastcon30way <- avg_score
    return(x)
  }, mc.cores = ncores)
  
  x <- do.call('rbind',x)
  return(x)
})

saveRDS(utr3_phastcon30way_perbin,paste0(meme_dir,fprefix,'_phastcon30way_perbin.rds'))

sessionInfo()
