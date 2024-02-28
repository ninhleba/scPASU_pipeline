library(Biostrings)
library(goldmine)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(argparser)
library(parallel)

parser<-arg_parser(name="nuc_freq_plot.R",description="Nucleotide frequency plotting")

parser<-add_argument(
  parser,
  arg='--cores',
  short = '-c',
  default=30,
  type='numeric',
  help="Number of cores to use when processing samples in parallel. Default = 30")

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
  arg='--ref',
  short = '-r',
  type="character",
  help="Enter the path to the peak ref. Format: path/to/peak_ref_file")

parser<-add_argument(
  parser,
  arg='--kmer',
  short = '-k',
  type="numeric",
  default=201,
  help="Length of region to plot nucleotide frequency. Ranges will be extended by half this length on both side.")

parser<-add_argument(
  parser,
  arg='--probs',
  short = '-p',
  type="character",
  default='0.05,0.5,0.95',
  help="Lower, Plotted value, Higher. Default: 0.05, 0.5, 0.95")

parser<-add_argument(
  parser,
  arg='--process_ref',
  short = '-b',
  flag=TRUE,
  help="Set flag if need to process standard ref to make GRanges object. The goal is to have chr, start, end and strand columns.")

parser<-add_argument(
  parser,
  arg='--min.gapwidth',
  short = '-d',
  type="numeric",
  default=NA,
  help="Min gap width to reduce ranges.")

parser<-add_argument(
  parser,
  arg='--kmer_to_plot',
  short = '-n',
  type="character",
  default="all",
  help="What kmer to plot. Values are all (default), first, or last.")

parser<-add_argument(
  parser,
  arg='--ref_genome',
  short = '-g',
  type="character",
  default="hg38",
  help="Your reference genome. Values are hg38 and mm10.")

args <- parse_args(parser)

cores<-args$cores %>% as.numeric
fprefix<-args$file_prefix
outdir<-args$outdir
ref<-args$ref
k<-args$kmer %>% as.numeric
probs <- args$probs %>% strsplit(split = ',') %>% unlist %>% as.numeric
extn<-(k-1)/2
process_ref<-args$process_ref
min.gapwidth <- args$min.gapwidth %>% as.numeric
kmer_to_plot <- args$kmer_to_plot
ref_genome <- args$ref_genome

if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

cat('Load peak rerefence from',ref,'\n')
if (process_ref) {
  cat('Rename columns...\n')
  file<-fread(ref,header=TRUE)
  file2<-file
  colnames(file2)<-gsub('\\bchr\\b','peak_chr',colnames(file2))
  colnames(file2)<-gsub('\\bstart\\b','peak_start',colnames(file2))
  colnames(file2)<-gsub('\\bend\\b','peak_end',colnames(file2))
  colnames(file2)<-gsub('\\bpr_chr\\b','chr',colnames(file2))
  colnames(file2)<-gsub('\\bpr_start\\b','start',colnames(file2))
  colnames(file2)<-gsub('\\bpr_end\\b','end',colnames(file2))
} else {
  file2 <- fread(ref,header=TRUE)
}

gr<-makeGRangesFromDataFrame(file2,ignore.strand = FALSE)

if (is.na(min.gapwidth)==FALSE) {
  cat('Merging ranges... Ranges separated by a gap of at least',min.gapwidth,'positions are not merged.\n')
  gr<-reduce(gr,min.gapwidth=min.gapwidth)
}

fl <- resize(gr,extn+width(gr),fix="end")
fl <- resize(fl,extn+width(fl),fix="start")
if (ref_genome == 'hg38') {
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,names=fl)
} else if (ref_genome == 'mm10') {
  seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,names=fl)
} else {stop('Acceptable values are hg38 and mm10')}

seq_list <- as.list(as.character(seq))

cat('Extract',k,'mers \n')

if (kmer_to_plot == 'all') {
  kmer_list <- lapply(seq_list,function(x) {
    substring(as.character(x), 1:(nchar(x) - k + 1), k:nchar(x))
  })
} else if (kmer_to_plot == 'first') {
  kmer_list <- lapply(seq_list,function(x) {
    substring(as.character(x), 1, k)
  })
} else if (kmer_to_plot == 'last') {
  kmer_list <- lapply(seq_list,function(x) {
    substring(as.character(x), (nchar(x) - k + 1), nchar(x))
  })
} else {stop('Acceptable values are all (default), first, or last \n')}

calculate_nucleotide_proportions <- function(kmers,k=201) {
  # Initialize a matrix to store nucleotide counts
  nuc <- c("A", "T", "C", "G", "N")

  nucleotide_counts <- matrix(0, nrow = k, ncol = length(nuc), 
                              dimnames = list(NULL, nuc))
  
  # Loop through each position in the k-mers
  for (i in 1:k) {
    # Extract nucleotides at the current position
    nucleotides <- sapply(kmers, function(x) substr(x, i, i))
    
    # Count the occurrences of each nucleotide
    nuc_freq <- table(nucleotides)
    missing_nuc <- nuc[-which(nuc %in% names(nuc_freq))]
    if (length(missing_nuc) > 0) {
      nuc_freq <- c(nuc_freq,rep(0,times=length(missing_nuc)))
      names(nuc_freq)[(length(nuc)-length(missing_nuc)+1):length(nuc)] <- missing_nuc
    }
    nuc_idx <- match(nuc,names(nuc_freq))
    nuc_freq <- nuc_freq[nuc_idx]
    
    nucleotide_counts[i, ] <- nuc_freq
  }
  
  # Convert counts to proportions
  nucleotide_proportions <- prop.table(nucleotide_counts, margin = 1) %>% t() %>% as.data.frame()
  nucleotide_proportions$nuc <- row.names(nucleotide_proportions)
  row.names(nucleotide_proportions) <- 1:length(nuc)
  
  return(nucleotide_proportions)
}

cat('Compute nucleotide frequency on each PR... \n')
proportion_list <- mclapply(kmer_list, function(x) {
  prop <- calculate_nucleotide_proportions(x,k=k)
  return(prop)
  }, mc.cores = cores
)

prop_df <- do.call('rbind',proportion_list)
prop_spl <- split(prop_df,prop_df$nuc)
saveRDS(prop_spl,paste0(outdir,fprefix,'_prop_by_nuc.rds'))

if (kmer_to_plot == 'all') {
  plot.data <- lapply(prop_spl,function(x) {
    stats <- apply(x[,1:k],2,function(x) {
      quantile(x,probs = c(probs))}) %>% t %>% as.data.frame
    colnames(stats) <- c('lower','plot','higher')
    return(stats)
  }
  )
  
  plot.data <- do.call('cbind',plot.data)
  plot.data$position <- 1:k
  
  g <- ggplot(plot.data, aes(x = position)) +
    geom_line(aes(y = A.plot, color = "A"), linewidth = 1) +
    geom_ribbon(aes(ymin = A.lower, ymax = A.higher), alpha = 0.1, fill = "#E41A1C") +
    geom_line(aes(y = C.plot, color = "C"), linewidth = 1) +
    geom_ribbon(aes(ymin = C.lower, ymax = C.higher), alpha = 0.1, fill = "#4DAF4A") +
    geom_line(aes(y = G.plot, color = "G"), linewidth = 1) +
    geom_ribbon(aes(ymin = G.lower, ymax = G.higher), alpha = 0.1, fill = "#377EB8") +
    geom_line(aes(y = T.plot, color = "T"), linewidth = 1) +
    geom_ribbon(aes(ymin = T.lower, ymax = T.higher), alpha = 0.1, fill = "#984EA3") +
    geom_line(aes(y = N.plot, color = "N"), linewidth = 1) +
    geom_ribbon(aes(ymin = N.lower, ymax = N.higher), alpha = 0.1, fill = "#FF7F00") +
    labs(x = "Position", y = "Nucleotide Frequency", title = fprefix) +
    scale_color_manual(values = c("A" = "#E41A1C", "C" = "#4DAF4A", "G" = "#377EB8", 
                                  "T" = "#984EA3", 'N' = "#FF7F00"),
                       name = "Nucleotide") +
    scale_x_continuous(breaks = c(1,51,101,151,201),
                       labels = c(-100,-50,0,50,100)) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic()
} else {
  plot.data <- lapply(prop_spl,function(x) {
    stats <- apply(x[,1:k],2,function(x) {
      sum(x)/length(x)}) %>% as.data.frame
    colnames(stats) <- c('plot')
    return(stats)
  })
  
  cn <- names(plot.data)
  plot.data <- do.call('cbind',plot.data)
  colnames(plot.data) <- paste0(cn,'.plot')
  plot.data$position <- 1:k
  
  g <- ggplot(plot.data, aes(x = position)) +
    geom_line(aes(y = A.plot, color = "A"), linewidth = 1) +
    geom_line(aes(y = C.plot, color = "C"), linewidth = 1) +
    geom_line(aes(y = G.plot, color = "G"), linewidth = 1) +
    geom_line(aes(y = T.plot, color = "T"), linewidth = 1) +
    geom_line(aes(y = N.plot, color = "N"), linewidth = 1) +
    labs(x = "Position", y = "Nucleotide Frequency", title = fprefix) +
    scale_color_manual(values = c("A" = "#E41A1C", "C" = "#4DAF4A", "G" = "#377EB8", 
                                  "T" = "#984EA3", 'N' = "#FF7F00"),
                       name = "Nucleotide") +
    scale_x_continuous(breaks = c(1,51,101,151,201),
                       labels = c(-100,-50,0,50,100)) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic()
}


ggsave(paste0(outdir,fprefix,'_nuc_freq.png'), width = 10, height = 8, units = 'in', g, bg = 'white')

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
