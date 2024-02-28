library(data.table)
library(dplyr)
library(stringr)
library(argparser)

parser<-arg_parser(name="TES_proximity_check.R",description="Check whether or not peaks without hexamers are close to at least the TES of annotated transcripts")

parser<-add_argument(
  parser,
  arg='--distance',
  short = '-b',
  default=100,
  type='numeric',
  help="Distance from TES. Default: 100 bp")

parser<-add_argument(
  parser,
  arg='--turef',
  short = '-t',
  type="character",
  default=NA,
  help="Enter the path to the TU reference file. We're going to use the transcript slot (genes.tu) to check proximity of peaks to annotated TES. Format: path/to/tu_ref_file")

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
  arg='--ref_file',
  short = '-r',
  type="character",
  default='filtered_peak_ref.txt',
  help="Enter the path to the peaks filtered by PR width in step 3b1. Format: path/to/ref_file. Default: filtered_peak_ref.txt")

args <- parse_args(parser)

distance=args$distance %>% as.numeric
genes_file=args$turef
fprefix=args$file_prefix
outdir=args$outdir
peaks_file=args$ref_file

w_hexamer_file=paste0(outdir,fprefix,'_peaks_with_hexamer.txt')
wo_hexamer_file=paste0(outdir,fprefix,'_peaks_wo_hexamer.txt')

cat(paste('Loading red_ens from',genes_file,'\n'),'We are going to use the transcript slot (genes.tu) to check proximity of peaks to annotated TES \n')
red_ens <- readRDS(genes_file)
genes <- red_ens$genes.tu %>% as.data.frame
genes$width <- NULL
colnames(genes)[1] <- 'chr'

genes <- split(genes,genes$tu_anno)

peaks <- fread(peaks_file)

w_hexamer <- fread(w_hexamer_file,header = FALSE)
colnames(w_hexamer) <- 'final_annotation'
w_hexamer$final_annotation <- gsub('^>','',w_hexamer$final_annotation)
w_hexamer$final_annotation <- str_sub(w_hexamer$final_annotation, start = 1, end = nchar(w_hexamer$final_annotation)-3)
cat(paste('Obtain full table of peaks that have hexamers, which are',nrow(w_hexamer),'in number \n'))
peaks_w_hexamer <- peaks[which(peaks$final_annotation %in% w_hexamer$final_annotation),]
write.table(peaks_w_hexamer, w_hexamer_file, sep='\t', row.names = FALSE, quote = FALSE)

wo_hexamer <- fread(wo_hexamer_file,header = FALSE)
colnames(wo_hexamer) <- 'final_annotation'
wo_hexamer$final_annotation <- gsub('^->','',wo_hexamer$final_annotation)
wo_hexamer$final_annotation <- str_sub(wo_hexamer$final_annotation, start = 1, end = nchar(wo_hexamer$final_annotation)-3)
peaks_wo_hexamer <- peaks[which(peaks$final_annotation %in% wo_hexamer$final_annotation),]

cols <- colnames(peaks)
tu_anno_col <- grep('tu_anno',cols)
strand_col <- grep('strand',cols)[1]
start_col <- grep('start',cols)[1]
end_col <- grep('end',cols)[1]
peak_col <- grep('final_annotation',cols)

cat(paste('Find peaks whose starts or ends are within the',distance,'bp distance of the TES of at least one annotated transcript \n'))
peaks_near_TES <- apply(peaks_wo_hexamer,1,function(x, dist = distance) {
  peak <- x[peak_col] %>% as.character
  tu_anno <- x[tu_anno_col] %>% as.character
  transcripts <- genes[[tu_anno]]
  strand <- x[strand_col] %>% as.character
  start <- x[start_col] %>% as.numeric
  end <- x[end_col] %>% as.numeric
  if (strand == '+') {
    transcripts <- transcripts[order(transcripts$end, decreasing = TRUE),]
    for (i in 1:nrow(transcripts)) {
      TES <- transcripts[i,'end'] %>% as.numeric
      start_to_TES <- abs(TES - start)
      end_to_TES <- abs(TES - end)
      if (start_to_TES <= dist | end_to_TES <= dist) {
        return(peak)
        break
      }
    }
  } else if (strand == '-') {
    transcripts <- transcripts[order(transcripts$start, decreasing = FALSE),]
    for (i in 1:nrow(transcripts)) {
      TES <- transcripts[i,'start'] %>% as.numeric
      start_to_TES <- abs(TES - start)
      end_to_TES <- abs(TES - end)
      if (start_to_TES <= dist | end_to_TES <= dist) {
        return(peak)
        break
      }
    }
  } else {stop('Strand has to be either + or -')}
})

peaks_near_TES <- unlist(peaks_near_TES)
peaks_near_TES_full_tbl <- peaks_wo_hexamer[which(peaks_wo_hexamer$final_annotation %in% peaks_near_TES),]
cat(paste('Out of',nrow(peaks_wo_hexamer),'peaks without hexamer,',nrow(peaks_near_TES_full_tbl),'are within',distance,'bp distance of the TES of at least one annotated transcript\n'))
write.table(peaks_near_TES_full_tbl, paste0(outdir,fprefix,'_peaks_wo_hexamer_near_TES.txt'), sep='\t', row.names = FALSE, quote = FALSE)

peaks_PAS_and_TES_filtered <- rbind(peaks_w_hexamer,peaks_near_TES_full_tbl)
cat(paste(nrow(peaks_PAS_and_TES_filtered),'peaks survive PAS and TES filtering \n'))
write.table(peaks_PAS_and_TES_filtered, paste0(outdir,fprefix,'_peaks_PAS_and_TES_filtered.txt'), sep='\t', row.names = FALSE, quote = FALSE)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')
