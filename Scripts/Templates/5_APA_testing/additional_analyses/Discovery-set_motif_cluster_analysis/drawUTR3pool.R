library(parallel)
library(dplyr)
library(goldmine)
library(argparser)

parser<-arg_parser(name="drawUTR3pool_arg.R",description="Randomly sample sequences from UTR3")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Working directory. Default=./")

parser<-add_argument(
  parser,
  arg='--turef',
  short = '-t',
  type="character",
  default=NA,
  help="Enter the path to the TU reference file. Format: path/to/tu_ref_file")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='utr3_ranges_NULL',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--bedfile',
  short = '-b',
  type="character",
  default='filename.bed',
  help="bedfile name. Default filename.bed")

parser<-add_argument(
  parser,
  arg='--ntimes',
  short = '-m',
  type="numeric",
  default=10,
  help="Ntimes x Query size = Sample size.")

parser<-add_argument(
  parser,
  arg='--ncores',
  short = '-n',
  type="numeric",
  default=32,
  help="Number of cores for parallel computing.")

args <- parse_args(parser)

outdir <- args$outdir
turef <- args$turef
ncores <- args$ncores
bedfile <- args$bedfile
fprefix <- args$file_prefix
n <- args$ntimes

setwd(outdir)

turef <- readRDS(turef)

utr3_range <- read.delim(bedfile,header = F)
colnames(utr3_range) <- c('chr','start','end','annotation','score','strand')
utr3_range <- makeGRangesFromDataFrame(utr3_range,ignore.strand = FALSE)
#utr3_range <- utr3_range[1:10]

utr3_ref <- as.data.frame(turef$utr3)
colnames(utr3_ref)[1] <- 'chr'

drawUTR3Pool <- function(target.gr,n=10,utr3_ref=utr3_ref,ncores=4) {
  chrs <- unique(seqnames(target.gr)) %>% as.character
  cat('Draw sequences from UTR3 on',paste0(chrs,collapse = ', '),'\n')
  utr3_ref$chr <- as.character(utr3_ref$chr)
  
  lens <- width(target.gr)
  lens <- rep(lens, n)
  chrs <- sample(chrs, length(lens), replace = T)
  dt <- data.frame(len = lens, chr = chrs)
  
  dodrawUTR3 <- function(dt,utr3_ref,ncores=ncores) {
    if (nrow(dt) != 0) {
      utr_colidx <- grep('utr',colnames(utr3_ref))
      width_colidx <- grep('width',colnames(utr3_ref))
      
      dt <- split(dt,row.names(dt))
      
      random_utr3_seq <- mclapply(dt,function(x) {
        len <- as.numeric(x$len)
        chr <- x$chr
        utr3_to_sample <- utr3_ref[which(utr3_ref$chr == chr & utr3_ref$width >= len),]
        utr3_pool <- apply(utr3_to_sample,1,function(x) {rep(x[utr_colidx],x[width_colidx])}) %>% unlist %>% unname
        utr3 <- sample(utr3_pool,1)
        utr3 <- utr3_ref[which(utr3_ref$utr == utr3),c('chr','start','end','strand','utr')]
        start <- sample(utr3$start:(utr3$end - len + 1),1)
        x$start <- start
        x$strand <- utr3$strand
        x$utr3 <- utr3$utr
        x$utr3_chr <- utr3$chr
        return(x)
      }, mc.cores = ncores)
      
      random_utr3_seq <- do.call('rbind',random_utr3_seq)
      
      if(!identical(random_utr3_seq$chr,random_utr3_seq$utr3_chr)) {
        print(random_utr3_seq)
        stop('The chromosome of the UTR3 sampled do not match the query chromosome')
      }
      
      random_utr3_seq <- random_utr3_seq %>% mutate(start = as.numeric(start), len = as.numeric(len))
      draw.gr <- GRanges(random_utr3_seq$chr, IRanges(start = random_utr3_seq$start, 
                                                      width = random_utr3_seq$len),
                         strand = random_utr3_seq$strand,
                         utr3 = random_utr3_seq$utr3,
                         utr3_chr = random_utr3_seq$utr3_chr)
      return(draw.gr)
    }
  }
  
  draw.gr <- dodrawUTR3(dt,utr3_ref,ncores)
  
  isbad <- TRUE
  iter <- 1
  while (isbad) {
    dups <- data.table(as.data.frame(findOverlaps(draw.gr, 
                                                  draw.gr)))
    dups <- dups[queryHits != subjectHits, ]$queryHits
    
    toget <- unique(dups)
    if (length(toget) == 0) {
      isbad <- FALSE
    }
    else {
      message("Draw iter ", iter, " had ", length(dups), 
              " overlapping. Redrawing these...")
      iter <- iter + 1
      dt.toget <- as.data.frame(draw.gr[c(toget)]) %>% select(width,seqnames)
      colnames(dt.toget) <- c('len','chr')
      dt.toget$chr <- as.character(dt.toget$chr)
      draw.gr <- c(draw.gr[-c(toget)], dodrawUTR3(dt.toget,utr3_ref,ncores))
    }
  }
  
  if (!all((table(width(target.gr)) * n) == table(width(draw.gr)))) {
    stop("Final length tables did not match, the drawing did not work.")
  }
  return(draw.gr)
}

utr3_range_NULL <- drawUTR3Pool(target.gr = utr3_range, n = n, utr3_ref = utr3_ref, ncores = ncores)
saveRDS(utr3_range_NULL,paste0(fprefix,'.rds'))

utr3_range_NULL <- as.data.frame(utr3_range_NULL)

bed <- utr3_range_NULL %>% mutate(score = 1, annotation = paste0('Seq',1:nrow(utr3_range_NULL))) %>% select(seqnames,start,end,annotation,score,strand)
write.table(bed,paste0(fprefix,'.bed'),sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

sessionInfo()
