library(dplyr)
library(goldmine)
library(parallel)
library(argparser)

parser<-arg_parser(name="3d_peak_classification.R",description="Classify peaks and identify split peaks to merge")

parser<-add_argument(
  parser,
  arg='--peaks',
  short = '-b',
  type="character",
  default=NA,
  help="Enter the peak file from the last step. Format: path/to/peaks_file")

parser<-add_argument(
  parser,
  arg='--outdir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='scPASU',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--turef',
  short = '-t',
  type="character",
  default=NA,
  help="Enter the path to the TU reference file. Format: path/to/tu_ref_file")

parser<-add_argument(
  parser,
  arg='--classification_ref',
  short = '-r',
  type="character",
  default=NA,
  help="Enter the path to the classification reference file. Will skip making one if provided. Format: path/to/classification_ref_file")

parser<-add_argument(
  parser,
  arg='--identify_split',
  short = '-z',
  flag=TRUE,
  help="Flag set to identify candidates for split peaks to merge and transcripts likely causing those fragmentations.")

parser<-add_argument(
  parser,
  arg='--dist_pct',
  short = '-p',
  default=0.1,
  type='numeric',
  help="Proportion of TU length from its TSS to be considered close. Default = 0.1")

parser<-add_argument(
  parser,
  arg='--frag_length',
  short = '-l',
  default=500,
  type='numeric',
  help="Fragment length to help identify transcripts likely to cause peak fragmentation. Default = 500")

parser<-add_argument(
  parser,
  arg='--ncores',
  short = '-n',
  default=4,
  type='numeric',
  help="Number of cores to use when processing in parallel. Default = 4")

args <- parse_args(parser)

peakref_file <- args$peaks
outdir <- args$outdir
fprefix <- args$file_prefix
turef_file <- args$turef
classref_file <- args$classification_ref
identify_split <- args$identify_split
dist.pct <- args$dist_pct %>% as.numeric
frag.length <- args$frag_length %>% as.numeric
ncores <- args$ncores

if(!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

peakref <- read.delim(peakref_file)
peakref <- split(peakref,peakref$tu)
peakref <- lapply(peakref, function(x) makeGRanges(x, strand = T))

if (is.na(turef_file) == FALSE) {
  turef <- readRDS(turef_file)
  } else {
    if (is.na(classref_file) == TRUE) {
      stop('TU reference file is not provided but classification reference file is not provided either')
    }
  }

### Create classification ref file ###
if(is.na(classref_file) == FALSE)
{
  cat(paste('Loading classification ref from',classref_file,'\n'))
  classification_ref<-readRDS(classref_file)
  utr3 <- classification_ref$utr3
  tss_range <- classification_ref$tss_range
  exons <- classification_ref$exons
  introns <- classification_ref$introns
}else{
  cat('Create classification reference from TU reference \n')
  ## TU & Flank ##
  cat('Adjust TU coordinate based on its flank region \n')
  tu <- turef$tu %>% as.data.frame()
  colnames(tu)[1] <- 'chr'
  tu$width <- NULL
  
  flank <- turef$flank %>% as.data.frame()
  
  tu_flank <- tu
  tu_flank$flank_start <- flank[match(tu_flank$tu,flank$tu),'start']
  tu_flank$flank_end <- flank[match(tu_flank$tu,flank$tu),'end']
  
  tu_flank$end <- ifelse(tu_flank$strand == '+',tu_flank$flank_end,tu_flank$end)
  tu_flank$start <- ifelse(tu_flank$strand == '-',tu_flank$flank_start,tu_flank$start)  
  tu_flank$tu_width <- tu_flank$end - tu_flank$start
  tu_flank <- split(tu_flank,tu_flank$tu)
  
  ## Exons ##
  cat('Build exon reference from the transcripts slot (genes.tu) in TU reference \n')
  genes.tu <- turef$genes.tu %>% as.data.frame()
  colnames(genes.tu)[1] <- 'chr'
  genes.tu$width <- NULL
  genes.tu <- split(genes.tu, genes.tu$tu)
  exons <- mclapply(genes.tu,function(x) {
    chr <- unique(x$chr)
    strand <- unique(x$strand)
    tu.id <- unique(x$tu)
    gene.id <- unique(x$gene.id)
    name <- unique(x$name)
    cols <- c('chr','start','end','strand','tu','gene.id','name','isoform.id')
    exon_count_idx <- grep('exonCount',colnames(x))
    exonStarts_idx <- grep('exonStarts',colnames(x))
    exonEnds_idx <- grep('exonEnds',colnames(x))
    isoform.id_idx <- grep('isoform.id',colnames(x))
    e <- apply(x,1,function(y) {
      t <- matrix(ncol = length(cols), nrow = as.numeric(y[exon_count_idx])) %>% as.data.frame()
      colnames(t) <- cols
      t$chr <- chr
      t$start <- strsplit(y[exonStarts_idx], split = ',') %>% unlist %>% as.numeric
      t$end <- strsplit(y[exonEnds_idx], split = ',') %>% unlist %>% as.numeric
      t$isoform.id <- y[isoform.id_idx]
      if (strand == '+') {
        t <- t[order(t$start, decreasing = FALSE),]
      } else if (strand == '-')  {t <- t[order(t$start, decreasing = TRUE),]
      }
      t$exon.id <- paste0(t$isoform.id,':Exon',1:nrow(t))
      return(t)
    })
    e <- do.call('rbind',e)
    e$strand <- strand
    e$tu <- tu.id
    e$gene.id <- gene.id
    e$name <- name
    return(e)
  }, mc.cores = ncores)
  
  exons_uniq <- mclapply(exons, function(x) {
    uniq <- x[,1:5] %>% group_by(tu) %>% distinct() %>% as.data.frame()
    uniq <- makeGRanges(uniq,strand=T)
    return(uniq)
  }, mc.cores = ncores)
  
  exons <- lapply(exons, function(x) makeGRanges(x,strand = T))
  
  ## TSS ##
  cat('Obtain range within',dist.pct,'the TU length from its TSS \n')
  tss_range <- mclapply(tu_flank, function(x) {
    tss <- ifelse(x$strand=='+',x$start,x$end) %>% as.numeric()
    dist = round(dist.pct*x$tu_width)
    x$start <- tss - dist
    x$end <- tss + dist
    x <- x[,c('chr','start','end','strand','tu','gene.id','name','tu_width')]
    return(makeGRanges(x,strand = T))
  }, mc.cores = ncores)
  
  
  ## UTR3 ##
  cat('Extract UTR3 reference from the UTR3 slot (utr3) in TU reference \n')
  utr3 <- turef$utr3 %>% as.data.frame()
  colnames(utr3)[1] <- 'chr'
  utr3$width <- NULL
  utr3 <- split(utr3,utr3$tu)
  utr3 <- mclapply(utr3,function(x) makeGRanges(x, strand = T))
  
  ## Internal ## (Regions outside UTR3)
  #cat('Extract internal (regions outside UTR3 on each TU) from the internal slot (internal) in TU reference \n')
  #internal <- turef$internal %>% as.data.frame()
  #colnames(internal)[1] <- 'chr'
  #internal$width <- NULL
  #cat('Genes without coding transcripts do not have UTR3 regions, so their internal regions are the entire TU \n')
  #genes.wo.internal <- setdiff(tu$tu,internal$tu)
  #internal <- rbind(internal,tu[tu$tu %in% genes.wo.internal,colnames(internal)])
  #internal <- makeGRanges(internal,strand=T)
  
  ### Introns
  cat('Build intron reference by extracting the regions that are not exonic and effectively not UTR3 \n')
  introns <-  mclapply(tu$tu,function(x) IRanges::setdiff(turef$tu[turef$tu$tu==x],exons_uniq[[x]]), mc.cores = ncores)
  names(introns) <- tu$tu
  
  classification_ref <- list(utr3=utr3,tss_range=tss_range,exons=exons,introns=introns)
  saveRDS(classification_ref,paste0(outdir,fprefix,'_classification_ref.rds'))
}

## Classify
cat('Classify PR sites by these categories: UTR3, TSS proximity, Intronic, Exonic and Flank \n')
peakref_class <- mclapply(peakref, function(x) {
  tu.id <- x$tu %>% unique
  x$utr3 <- NA
  x$tss.prox <- NA
  x$exonic <- NA
  x$intronic <- NA
  x$flank <- NA
  x$isoform <- NA
  x$exon <- NA
  if (!is.null(utr3[[tu.id]])) {
    ovl <- findOverlaps(x,utr3[[tu.id]]) %>% as.data.frame()
    if (nrow(ovl) != 0) {x[unique(ovl$queryHits)]$utr3 <- 1}
  }
  if (!is.null(tss_range[[tu.id]])) {
    ovl <- findOverlaps(x,tss_range[[tu.id]]) %>% as.data.frame()
    if (nrow(ovl) != 0) {x[unique(ovl$queryHits)]$tss.prox <- 2}
  }
  if (!is.null(exons[[tu.id]])) {
    e <- exons[[tu.id]]
    ovl <- findOverlaps(x,e) %>% as.data.frame()
    if (nrow(ovl) != 0) {
      x[unique(ovl$queryHits)]$exonic <- 3
      for (peak in unique(ovl$queryHits)) {
        x[peak]$isoform <- paste0(e[subset(ovl,queryHits==peak)$subjectHits]$isoform.id,collapse=',')
        x[peak]$exon <- paste0(e[subset(ovl,queryHits==peak)$subjectHits]$exon.id,collapse=',') #One exon can appear in multiple transcripts, so can be duplicates here
      }
    }
  }
  if (!is.null(introns[[tu.id]])) {
    ovl <- findOverlaps(x,introns[[tu.id]]) %>% as.data.frame()
    if (nrow(ovl) != 0) {x[unique(ovl$queryHits)]$intronic <- 4}
  }
  x$flank[which(nchar(x$flank_tus) != 0)] <- 5
  return(x)
}, mc.cores = ncores)

peakref_class_df <- mclapply(peakref_class, function(x) {
  x <- x[order(x)]
  x <- as.data.frame(x)
  colnames(x)[1] <- 'chr'
  x$width <- NULL
  
  x <- x %>% mutate(top_classification = pmin(utr3,tss.prox,exonic,intronic,flank,na.rm = TRUE))
  x <- x %>%
    mutate(final_classification = case_when(top_classification == 1 ~ "UTR3",
                                            top_classification == 2 ~ "TSS-proximal",
                                            top_classification == 3 ~ "Exonic",
                                            top_classification == 4 ~ "Intronic",
                                            top_classification == 5 ~ "Flank",
                                            TRUE ~ "Unclassified"
                                         ))
  x$top_classification <- NULL
  return(x)
}, mc.cores = ncores)

peakref_class_df <- do.call('rbind',peakref_class_df)

cat("Classification outcomes: \n")
print(table(peakref_class_df$final_classification))

write.table(peakref_class_df,paste0(outdir,fprefix,'_peak_universe_classified.txt'),sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

if (identify_split) {
  cat('Identify candidates for split peaks to merge and transcripts likely causing those fragmentations \n')
  spliced_candidate <- mclapply(peakref_class, function(x) {
    if (length(x) > 1) {
      all_isoforms <- x$isoform
      all_isoforms <- all_isoforms[!is.na(all_isoforms)]
      all_isoforms <- paste0(all_isoforms,collapse=',')
      all_isoforms <- all_isoforms %>% strsplit(split = ',') %>% unlist
      multipeak_isoforms <- all_isoforms[duplicated(all_isoforms)] %>% unique #Because isoforms cannot duplicate within each row, duplication has to come from different peaks
      if (length(multipeak_isoforms) > 0) {
        tu.id <- x$tu %>% unique
        mi <- exons[[tu.id]][exons[[tu.id]]$isoform.id %in% multipeak_isoforms]
        for (i in 1:length(mi)) {
          exon_pattern <- paste0('\\b',mi$exon.id[i],'\\b')
          peak_ovl_exon <- x$final_annotation[grep(exon_pattern,x$exon)]
          mi$peak[i] <- paste0(peak_ovl_exon,collapse=',')
        }
        return(mi)
      }
    }
  }, mc.cores = ncores)
  
  spliced_candidate_null.rm <- spliced_candidate[!(lapply(spliced_candidate, function(x) is.null(x)) %>% unlist)]
  spliced_candidate_null.rm <- lapply(spliced_candidate_null.rm, function(x) as.data.frame(x))
  spliced_candidate_null.rm <- do.call('rbind',spliced_candidate_null.rm)
  
  cat('Remove transcripts whose peak-overlapped exons sum up to a length bigger than fragment length',frag.length,'bp\n')
  cat('Also remove transcripts whose peak-overlapped exons overlap with only one peak\n')
  spliced_candidate_null.rm_filtered <- split(spliced_candidate_null.rm, spliced_candidate_null.rm$isoform.id)
  
  spliced_candidate_null.rm_filtered <- lapply(spliced_candidate_null.rm_filtered, function(x) {
    exons_ovl_idx <- sapply(x$peak, function(x) nchar(x)!=0)
    exons_ovl_total_len <- x$width[exons_ovl_idx] %>% sum()
    ovl_peaks <- x$peak[(nchar(x$peak)!=0)]
    ovl_peaks <- paste0(ovl_peaks,collapse=',') %>% strsplit(split = ',') %>% unlist
    ovl_peaks <- unique(ovl_peaks)
    if (exons_ovl_total_len <= frag.length & length(ovl_peaks) > 1 & length(exons_ovl_idx) > 1) {
      return(x)
    }
  })
  spliced_candidate_null.rm_filtered <- spliced_candidate_null.rm_filtered[!(lapply(spliced_candidate_null.rm_filtered, 
                                                                                    function(x) is.null(x)) %>% unlist)]
  spliced_candidate_null.rm_filtered <- do.call('rbind',spliced_candidate_null.rm_filtered)
  write.table(spliced_candidate_null.rm_filtered,paste0(outdir,fprefix,'_transcripts_causing_peak_fragmentation.txt'),
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')