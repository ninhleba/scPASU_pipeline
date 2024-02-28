## Adapted from https://github.com/hwanglab/apa_atingLab2019 ##

library(Rsamtools)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape)
library(matrixStats)
library(edgeR)
library(DEXSeq)
library(GenomicRanges)
library(goldmine)
library(rtracklayer)
library(dplyr)
library(gtools)

### This function creates TUs ###
reduceGenes <- function(genes,chrs,flank=5000,ncore=1)
{
  
  genes <- makeGRanges(genes,strand=T)
  genes <- genes[seqnames(genes) %in% chrs]
  genes$cdsStart <- genes$cdsStart+1
  genes$cdsEnd <- genes$cdsEnd+1
  stopifnot(genes$strand %in% c("+","-"))
  
  # Some genes may have disjoint txUnits, let's treat them as separate genes to make things easier
  cat("Clustering TUs \n")
  df <- as.data.frame(table(genes$gene.id))
  genes.singles <- genes[genes$gene.id %in% df[df$Freq==1,]$Var1]
  genes.tored <- genes[genes$gene.id %in% df[df$Freq>1,]$Var1]
  spl <- split(genes.tored,genes.tored$gene.id)
  spl.red <- lapply(spl,function(x) reduce(x,with.revmap=FALSE))#collect union boundary of the same gene
  for(i in 1:length(spl.red))
  {
    spl.red[[i]]$gene.id <- names(spl)[i]
  }
  
  genes.red <- do.call(c,unname(spl.red))
  genes.red$name <- genes[match(genes.red$gene.id,genes$gene.id)]$name
  cat("Single isoform annotated:",length(unique(genes.singles$gene.id)),'\n')
  cat("Multiple isoform annotated:",length(unique(genes.tored$gene.id)),'\n')
  dt <- makeDT(genes.red)
  genes.red.counts <- dt[,list(nPos=length(strand),nStrands=length(unique(strand))),by="gene.id"]
  cat("Reduced to disjoint clusters on same strand only:",nrow(genes.red.counts[(nPos>1)&(nStrands==1),]),'\n')
  cat("Reduced to disjoint clusters on different strands:",nrow(genes.red.counts[(nStrands>1),]),'\n')
  # Assign isoforms to TUs
  tu <- genes.singles
  values(tu) <- NULL # removes all metadata
  tu$gene.id <- genes.singles$gene.id
  tu$name <- genes.singles$name
  tu <- c(tu,genes.red)
  seqlevels(tu) <- chrs
  tu <- tu[order(tu)]
  tu$tu <- paste0("TU",1:length(tu))
  tu$tu_anno<-paste0(tu$tu,':',tu$name)
  
  # From now on, must aggregate by TU rather than gene.id, as there can be duplication in gene.id 
  cat("Assigning transcripts to TUs \n")
  genes.singles$tu <- tu[match(genes.singles$gene.id,tu$gene.id),]$tu
  genes.singles$tu_anno <- tu[match(genes.singles$gene.id,tu$gene.id),]$tu_anno
  
  fo <- data.table(as.data.frame(findOverlaps(genes.tored,tu)))
  fo$gene_gid <- genes.tored[fo$queryHits]$gene.id
  fo$tu_gid <- tu[fo$subjectHits]$gene.id
  fo$assign <- fo$gene_gid==fo$tu_gid
  fo <- fo[assign==TRUE,] #For the one w/ FALSE, a different gene is overlaped. The action with assign==TRUE still retains the one w/ overlapped with the other gene
  stopifnot(!duplicated(fo$queryHits))
  stopifnot(nrow(fo)==length(genes.tored))
  dt <- makeDT(genes.tored)
  dt[fo$queryHits,tu:= tu[fo$subjectHits]$tu] 
  dt[,tu_anno:= paste0(dt$tu,':',dt$name)] 
  genes.tu <- c(genes.singles,makeGRanges(dt,strand=TRUE))
  
  stopifnot(length(unique(genes.tu$tu))==length(tu))
  stopifnot(length(genes.tu)==length(genes))
  stopifnot(length(unique(genes$gene.id))==length(unique(genes.tu$gene.id)))
  seqlevels(genes.tu) <- chrs
  genes.tu <- genes.tu[order(genes.tu)]
  
  # Assign TUs as coding
  cat("Detecting coding TUs \n")
  genes.tu$coding <- genes.tu$cdsStart!=genes.tu$cdsEnd
  dt <- makeDT(genes.tu)
  dt <- dt[,list(coding=any(coding)),by="tu"]
  stopifnot(nrow(dt)==length(tu))
  tu$coding <- "new"
  tu[match(dt$tu,tu$tu)]$coding <- dt$coding
  
  # Get 3' end flanks (i.e., 5k bp segment right after TSE)
  cat("Generating TU flanks \n")
  tu.flank <- flank(tu,width=flank,start=FALSE,both=FALSE,ignore.strand=FALSE)
  
  # Get 3' UTRs
  cat("Clustering 3' UTRs for each coding TU \n")
  coding <- genes.tu[genes.tu$coding==TRUE]
  
  # If on (+), then want range between cdsEnd and txEnd
  # If on (-), then want range between txStart and cdsStart
  coding.p <- coding[strand(coding)=="+"]
  utr3.p <- GRanges(seqnames(coding.p),IRanges(coding.p$cdsEnd,end(coding.p)),strand=strand(coding.p),tu=coding.p$tu,gene.id=coding.p$gene.id,name=coding.p$name)
  coding.m <- coding[strand(coding)=="-"]
  utr3.m <- GRanges(seqnames(coding.m),IRanges(start(coding.m),coding.m$cdsStart),strand=strand(coding.m),tu=coding.m$tu,gene.id=coding.m$gene.id,name=coding.m$name)
  stopifnot((length(coding.p)+length(coding.m))==length(coding))
  
  redme <- function(myutrs)
  {
    df <- as.data.frame(table(myutrs$tu))
    myutrs.singles <- myutrs[myutrs$tu %in% df[df$Freq==1,]$Var1]
    myutrs.tored <- myutrs[myutrs$tu %in% df[df$Freq>1,]$Var1]
    spl <- split(myutrs.tored,myutrs.tored$tu)
    spl.red <- lapply(spl,reduce)
    for(i in 1:length(spl.red))
    {
      spl.red[[i]]$tu <- names(spl)[i]
    }
    myutrs.red <- do.call(c,unname(spl.red))
    myutrs.red$gene.id <- myutrs[match(myutrs.red$tu,myutrs$tu)]$gene.id
    myutrs.red$name <- myutrs[match(myutrs.red$tu,myutrs$tu)]$name
    c(myutrs.singles,myutrs.red)
  }
  utr3.red.m <- redme(utr3.m)
  utr3.red.p <- redme(utr3.p)
  utr3 <- c(utr3.red.p,utr3.red.m)
  utr3 <- utr3[order(utr3)]
  utr3 <- utr3[width(utr3)>0]
  utr3$utr <- paste0("UTR",1:length(utr3))
  
  # Internal
  cat("Computing coding gene internal regions via setdiff \n")
  spl <- split(utr3,utr3$tu)
  coding <- tu[tu$coding==TRUE]
  
  # Remove coding genes that don't have 3' UTRs
  coding <- coding[coding$tu %in% names(spl)]
  diff <- lapply(coding$tu,function(x) setdiff(coding[coding$tu==x],spl[[x]]))
  for(i in 1:length(diff))
  {
    diff[[i]]$tu <- coding$tu[i]
  }
  
  internal <- do.call(c,diff)
  internal$gene.id <- tu[match(internal$tu,tu$tu)]$gene.id
  internal$name <- tu[match(internal$tu,tu$tu)]$name
  
  # No-UTR gene bodies
  coding <- tu[tu$coding==TRUE]
  noutr <- coding[!(coding$tu %in% names(spl))]
  
  # Between annotated 3' UTRs
  cat("Computing between multi-UTR \n")
  dt <- makeDT(utr3)
  dt <- dt[,list(chr=chr[1],start=min(start),end=max(end),strand=strand[1],gene.id=gene.id[1],name=name[1],nUtr=length(utr)),by="tu"]
  multis <- dt[nUtr>1,]
  # Now have the maximal range for all the multis
  # If we setdiff out the real 3', then we have just the between ranges
  union <- makeGRanges(multis,strand=T)
  spl <- split(utr3,utr3$tu)
  diff <- lapply(union$tu,function(x) setdiff(union[union$tu==x],spl[[x]]))
  for(i in 1:length(diff))
  {
    diff[[i]]$tu <- union$tu[i]
  }
  
  btw <- do.call(c,diff)
  btw$gene.id <- tu[match(btw$tu,tu$tu)]$gene.id
  btw$name <- tu[match(btw$tu,tu$tu)]$name
  
  # Now remove the betweens from the internals
  internal <- internal[!(internal %in% btw)]
  
  ret <- list(tu=tu,genes.tu=genes.tu,flank=tu.flank,utr3=utr3,between=btw,internal=internal,noutr3=noutr)
  
  #writeBEDFromGRanges(ret$tu,file="output/tu_tu.bed",name="tu")
  #writeBEDFromGRanges(ret$flank,file="output/tu_flank.bed",name="tu")
  #writeBEDFromGRanges(ret$utr3,file="output/tu_utr3.bed",name="tu")
  #writeBEDFromGRanges(ret$between,file="output/tu_between.bed",name="tu")
  #writeBEDFromGRanges(ret$internal,file="output/tu_internal.bed",name="tu")
  #writeBEDFromGRanges(ret$noutr3,file="output/tu_noutr3.bed",name="tu")
  return(ret)
}
### This function assigns peaks to TUs ###

joinTus_peaks<- function(allpeaks,rg)
{
  # Retain only polya supported peaks
  cat('Retain only polyA supported peaks. This step should be redundant because the input peaks are only polyA peaks from the earlier step \n')
  keep<-which(allpeaks$polya=='polya')
  peaks<-allpeaks[keep,]
  
  cat('Annotate peaks based on which TU or flanking region they lie in \n')
  fo_tu <- data.table(as.data.frame(findOverlaps(peaks,rg$tu)))
  fo_flank <- data.table(as.data.frame(findOverlaps(peaks,rg$flank)))
  
  fo_tu$set <- "tu"
  fo_flank$set <- "flank"
  
  fo_tu$tu <- rg$tu[fo_tu$subjectHits]$tu
  fo_tu$coding <- rg$tu[fo_tu$subjectHits]$coding
  cat('Assigning TU anno \n')
  fo_tu$tu_anno <- rg$tu[fo_tu$subjectHits]$tu_anno
  
  fo_flank$tu <- rg$flank[fo_flank$subjectHits]$tu
  fo_flank$coding <- rg$flank[fo_flank$subjectHits]$coding
  cat('Assigning TU anno from flank \n')
  fo_flank$tu_anno <- rg$flank[fo_flank$subjectHits]$tu_anno
  
  # Remove cases where same PR links to same TU and flank of that TU
  cat('For peaks that lie at the junction of a TU and its 3 flanking region, count them as part of the TU body \n')
  fo_flank$key <- paste0(fo_flank$queryHits,"+",fo_flank$tu)
  fo_tu$key <- paste0(fo_tu$queryHits,"+",fo_tu$tu)
  fo_flank <- fo_flank[!(fo_flank$key %in% fo_tu$key),]
  
  ov <- rbind(fo_tu,fo_flank)
  
  # Join table that links each PR to a TU/flank
  cat('Create data table \n')
  join <- data.table(peak=peaks[ov$queryHits]$peakID,tu=ov$tu,tu_anno=ov$tu_anno, type=ov$set,coding=ov$coding)
  
  # Reduced assignment table that assigns the uniques
  cat('Check if a peak is assigned to only one TU, i.e. unique_peak and if a TU contains only unique peaks, i.e. unique_tu \n')
  join <- join[,list(tu=tu,tu_anno=tu_anno,type=type,coding=coding,unique_peak=length(tu)==1,over_tus=toString(tu[type=="tu"]),flank_tus=toString(tu[type=="flank"])),by="peak"]
  join <- join[,list(peak=peak,type=type,coding=coding,unique_peak=unique_peak,unique_tu=all(unique_peak),over_tus=over_tus,flank_tus=flank_tus,tu_anno=tu_anno),by="tu"]
  
  ret <-list(allpeaks=allpeaks,polya_peaks=peaks,join=join)
  return(ret)
}

create_final_annotation_col<-function(mtu, is_minus=FALSE)
{
  # sort
  mtu<-mtu[mixedorder(mtu$tu),]
  
  # Sort minus strand differently for correct annotation
  if(is_minus==TRUE)
  {
    mtu<-dplyr::arrange(mtu,chr,desc(start))
  }
  
  tu<-unique(mtu$tu)
  cnt<-mtu %>% group_by(tu) %>% tally()
  cnt<-cnt[order(cnt$n),]
  
  mtu$final_annotation<-rep('none',nrow(mtu))
  
  for(i in 1:length(tu))
  {
    
    indx<-match(tu[i],cnt$tu)
    
    n<-cnt$n[indx]
    p<-paste0('P',(1:n))
    indx2<-which(mtu$tu %in% tu[i])
    new_anno<-paste0(mtu$tu[indx2],':',mtu$gene[indx2],':',p)
    mtu$final_annotation[indx2]<-new_anno
    
  }
  return(mtu)
  
}


get_histogram<-function(dat,binw,col='royal blue',fill='salmon', x_lim='none',y_lim='none',title,x_lab)
{
  if(length(x_lim)==1)
  {
    x_lim<-c(0,max(dat))
  }
  rem<-dat[dat>x_lim[2]] %>% length()
  rem_perc<-((rem/length(dat)) *100) %>% round(.,digits = 1)
  
  if(length(x_lim)==1)
  {h<-ggplot()+
    geom_histogram(aes(dat),binwidth=binw, colour=col,fill= fill)+
    xlim(x_lim)
  }else{
    h<-ggplot()+
      geom_histogram(aes(dat),binwidth=binw, colour=col,fill= fill)+
      xlim(x_lim) + ylim(y_lim)
  }
  
  min<-min(dat, na.rm = TRUE)
  max<-max(dat, na.rm = TRUE)
  
  h<-h+labs(title = title,
            subtitle = paste0('[binwidth:',binw,']', ' [min:',min,'; max:', max,']'),
            caption=paste0('[',rem_perc,'% datapoints removed by the given cutoff]'))+
    xlab(x_lab)+theme_bw()+
    theme(plot.caption=element_text(color='navy blue',
                                    size=8),
          text = element_text(size=30))
  return(h)
}


# These two functions convert cellranger reference gtf into a TU reference

reformat_gtf <- function(gtf, chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                                       'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                                       'chr20','chr21','chr22','chrX','chrY'), 
                         isoform.rm = NULL,
                         biotype.rt = c("lncRNA","protein_coding","IG_V_gene","IG_C_gene",
                                        "IG_J_gene","TR_C_gene","TR_J_gene","TR_V_gene",
                                        "TR_D_gene","IG_D_gene")) {
  cat('Retain isoforms on these chromosomes only:',paste(chrs, collapse = ', '),'\n')
  gtf <- gtf[gtf$seqnames %in% chrs,]
  
  if (is.null(biotype.rt) == FALSE) {
    cat('Retain isoforms of these biotypes only:',paste(biotype.rt, collapse = ', '),'\n')
    gtf_nogenes <- gtf[(gtf$transcript_type %in% biotype.rt),]
    gtf_nogenes$transcript_id <- paste(gtf_nogenes$transcript_id,gtf_nogenes$transcript_version,sep = '.')
  }
  
  #Filtering isoforms
  if (is.null(isoform.rm) == FALSE) {
    cat('Remove these isoforms: ',paste(isoform.rm, collapse = ', '),'\n')
    gtf_nogenes <- gtf_nogenes[!(gtf_nogenes$transcript_id %in% isoform.rm),]
  }
  
  gtf_nogenes$seqnames <- gtf_nogenes$seqnames %>% as.character
  gtf_nogenes$strand <- gtf_nogenes$strand %>% as.character
  
  cat('Convert gtf file from long to wide format... \n')
  spl <- split(gtf_nogenes,gtf_nogenes$transcript_id)
  genes <- lapply(spl, function(gtf_transcript_id) {
    genes_row <- matrix(nrow = 1, ncol = 20) %>% as.data.frame()
    colnames(genes_row) <- c("chr","start","end","strand","name","gene.id","isoform.id","cdsStart","cdsEnd","num_startcodon","num_stopcodon","startcodon_start","startcodon_end","stopcodon_start","stopcodon_end","utrStart","utrEnd","exonCount","exonStarts","exonEnds")
    
    genes_row[,c("chr","start","end","strand","name","gene.id","isoform.id")] <- gtf_transcript_id[which(gtf_transcript_id$type == 'transcript'),c('seqnames','start','end','strand','gene_name','gene_id','transcript_id')]
    
    if ('CDS' %in% gtf_transcript_id$type) {
      genes_row[,c("cdsStart")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'CDS'),]$start,collapse=',')
      genes_row[,c("cdsEnd")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'CDS'),]$end,collapse=',')
      stopifnot(strsplit(genes_row[,"cdsStart"],split = ',') %>% unlist %>% length == which(gtf_transcript_id$type == 'CDS') %>% length)
    }
    if ('start_codon' %in% gtf_transcript_id$type) {
      genes_row[,c("startcodon_start")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'start_codon'),]$start,collapse=',')
      genes_row[,c("startcodon_end")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'start_codon'),]$end,collapse=',')
      genes_row[,c("num_startcodon")] <- gtf_transcript_id[which(gtf_transcript_id$type == 'start_codon'),]$start %>% length
    }
    if ('stop_codon' %in% gtf_transcript_id$type) {
      genes_row[,c("stopcodon_start")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'stop_codon'),]$start,collapse=',')
      genes_row[,c("stopcodon_end")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'stop_codon'),]$end,collapse=',')
      genes_row[,c("num_stopcodon")] <- gtf_transcript_id[which(gtf_transcript_id$type == 'stop_codon'),]$start %>% length
    }
    if ('UTR' %in% gtf_transcript_id$type) {
      genes_row[,c("utrStart")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'UTR'),]$start,collapse=',')
      genes_row[,c("utrEnd")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'UTR'),]$end,collapse=',')
      stopifnot(strsplit(genes_row[,"utrStart"],split = ',') %>% unlist %>% length == which(gtf_transcript_id$type == 'UTR') %>% length)
    }
    stopifnot(gtf_transcript_id$exon_number %>% as.numeric %>% max(na.rm=TRUE) == which(gtf_transcript_id$type == 'exon') %>% length)
    genes_row[,c("exonCount")] <- gtf_transcript_id$exon_number %>% as.numeric %>% max(na.rm=TRUE)
    genes_row[,c("exonStarts")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'exon'),]$start,collapse=',')
    genes_row[,c("exonEnds")] <- paste(gtf_transcript_id[which(gtf_transcript_id$type == 'exon'),]$end,collapse=',')
    stopifnot(strsplit(genes_row[,"exonStarts"],split = ',') %>% unlist %>% length == genes_row[,"exonCount"])
    
    return(genes_row)
  })
  genes_df <- do.call('rbind',genes)
  
  #### Prepare UTR3
  genes_coding <- genes_df[!(is.na(genes_df$cdsStart)),]
  
  ### Plus strand 
  ## Have UTR annotation
  cat('Processing UTR3 regions on plus strand... \n')
  genes_coding_p_utr <- genes_coding[which(genes_coding$strand=='+' & is.na(genes_coding$utrStart) == FALSE),]
  genes_coding_p_utr_spl <- split(genes_coding_p_utr, genes_coding_p_utr$isoform.id)
  utr3_coding_p <- lapply(genes_coding_p_utr_spl, function(row) {
    utr3 <- matrix(nrow = 1, ncol = 16) %>% as.data.frame()
    colnames(utr3) <- c("chr","start","end","strand","name","gene.id","isoform.id","cdsStart","cdsEnd","utrStart","utrEnd","method","end_w_transcript","exonCount","exonStarts","exonEnds")
    utr3[,c("chr","start","end","strand","name","gene.id","isoform.id","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")] <- row[,c("chr","start","end","strand","name","gene.id","isoform.id","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")]
    utr3[,c("start","end","exonCount")] <- utr3[,c("start","end","exonCount")] %>% as.numeric
    
    utrstarts <- row$utrStart %>% strsplit(split = ',') %>% unlist %>% as.numeric
    utrends <- row$utrEnd %>% strsplit(split = ',') %>% unlist %>% as.numeric
    stopifnot(length(utrstarts)==length(utrends))
    
    cdsends <- utr3$cdsEnd %>% strsplit(split = ',') %>% unlist %>% as.numeric
    max_cdsend <- max(cdsends)
    
    if (!is.na(row$num_stopcodon)) {
      stopcodonstarts <- row$stopcodon_start %>% strsplit(split = ',') %>% unlist %>% as.numeric
      stopifnot(any(stopcodonstarts < max_cdsend) == FALSE) #Sanity check if all the stop codon starts are after the last cds end
      if (is.unsorted(stopcodonstarts) == FALSE) {stopcodonstarts <- sort(stopcodonstarts)} # Make sure stop codon starts are in increasing order
      for (i in 1:length(stopcodonstarts)) {
        if (match(stopcodonstarts[i],utrstarts) %>% is.na == FALSE) {
          match_index <- match(stopcodonstarts[i],utrstarts)
          utr3_candidate_start <- utrstarts[match_index]
          match_index <- c(match_index, which(utrstarts > utr3_candidate_start))
          utr3_candidate_start <- utrstarts[match_index]
          
          utr3_candidate_end <- utrends[match_index]
          utr3[,c("utrStart","utrEnd")] <- c(paste(utr3_candidate_start,collapse = ','),paste(utr3_candidate_end,collapse=','))
          utr3$end_w_transcript <- max(utr3_candidate_end) == utr3$end
          utr3$method <- 'stop_codon'
          break 
        }
      }
    } else { #no stop codon 
      match_index <- which(utrstarts >= max_cdsend) 
      if (length(match_index) != 0) {
        utr3_candidate_start <- utrstarts[match_index]
        utr3_candidate_end <- utrends[match_index]
        utr3[,c("utrStart","utrEnd")] <- c(paste(utr3_candidate_start,collapse = ','),paste(utr3_candidate_end,collapse=','))
        utr3$end_w_transcript <- max(utr3_candidate_end) == utr3$end
        utr3$method <- 'coding_sequence'
      } else {
        utr3_candidate_start <- max_cdsend + 1
        if (max(utrends) > utr3_candidate_start) {
          utr3_candidate_end <- max(utrends)
          method <- 'substract_utr_by_cds'
        } else {
          utr3_candidate_end <- utr3$end
          method <- 'substract_transcript_by_cds'
        }
        utr3[,c("utrStart","utrEnd")] <- c(paste(utr3_candidate_start,collapse = ','),paste(utr3_candidate_end,collapse=','))
        utr3$end_w_transcript <- utr3_candidate_end == utr3$end
        utr3$method <- method
      }
    } 
    return(utr3)
  })
  
  utr3_coding_p_df <- do.call('rbind',utr3_coding_p)
  cat('Isoforms with UTR3 regions computed using stop_codon method:',length(which(utr3_coding_p_df$method == 'stop_codon')),'\n')
  cat('Isoforms with UTR3 regions computed using coding_sequence method:',length(which(utr3_coding_p_df$method == 'coding_sequence')),'\n')
  cat('Isoforms with UTR3 regions computed using substract_transcript_by_cds method:',length(which(utr3_coding_p_df$method == 'substract_transcript_by_cds')),'\n')
  
  ### Minus strand
  ## Have UTR annotation
  cat('Processing UTR3 regions on minus strand...')
  genes_coding_m_utr <- genes_coding[which(genes_coding$strand=='-' & is.na(genes_coding$utrStart) == FALSE),]
  genes_coding_m_utr_spl <- split(genes_coding_m_utr, genes_coding_m_utr$isoform.id)
  utr3_coding_m <- lapply(genes_coding_m_utr_spl, function(row) {
    utr3 <- matrix(nrow = 1, ncol = 16) %>% as.data.frame
    colnames(utr3) <- c("chr","start","end","strand","name","gene.id","isoform.id","cdsStart","cdsEnd","utrStart","utrEnd","method","start_w_transcript","exonCount","exonStarts","exonEnds")
    utr3[,c("chr","start","end","strand","name","gene.id","isoform.id","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")] <- row[,c("chr","start","end","strand","name","gene.id","isoform.id","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")]
    utr3[,c("start","end","exonCount")] <- utr3[,c("start","end","exonCount")] %>% as.numeric
    
    utrstarts <- row$utrStart %>% strsplit(split = ',') %>% unlist %>% as.numeric
    #  utrstarts_de <- sort(utrstarts, decreasing = TRUE)
    #  stopifnot(identical(utrstarts,utrstarts_de)) #Check if UTR ends are in decreasing order
    
    utrends <- row$utrEnd %>% strsplit(split = ',') %>% unlist %>% as.numeric
    #  utrends_de <- sort(utrends, decreasing = TRUE)
    #  stopifnot(identical(utrends,utrends_de)) #Check if UTR ends are in decreasing order
    
    stopifnot(length(utrstarts)==length(utrends))
    
    cdsstarts <- utr3$cdsStart %>% strsplit(split = ',') %>% unlist %>% as.numeric
    min_cdsstarts <- min(cdsstarts)
    
    if (!is.na(row$num_stopcodon)) {
      stopcodonends <- row$stopcodon_end %>% strsplit(split = ',') %>% unlist %>% as.numeric
      stopifnot(any(stopcodonends > min_cdsstarts) == FALSE) #Sanity check if all the stop codon ends are before the first cds start
      stopcodonends_de <- sort(stopcodonends, decreasing = TRUE)
      if (identical(stopcodonends,stopcodonends_de) == FALSE) {stopcodonends <- stopcodonends_de} #Make sure stop codon ends are in decreasing order
      
      for (i in 1:length(stopcodonends)) {
        if (match(stopcodonends[i],utrends) %>% is.na == FALSE) {
          match_index <- match(stopcodonends[i],utrends)
          utr3_candidate_end <- utrends[match_index]
          match_index <- c(match_index,which(utrends < utr3_candidate_end))
          utr3_candidate_end <- utrends[match_index]
          
          utr3_candidate_start <- utrstarts[match_index]
          utr3[,c("utrEnd","utrStart")] <- c(paste(utr3_candidate_end,collapse = ','),paste(utr3_candidate_start,collapse=','))
          utr3$start_w_transcript <- min(utr3_candidate_start) == utr3$start
          utr3$method <- 'stop_codon'
          break 
        }
      }
    } else { #no stop codon 
      match_index <- which(utrends <= min_cdsstarts)
      if (length(match_index) != 0) {
        utr3_candidate_end <- utrends[match_index]
        utr3_candidate_start <- utrstarts[match_index]
        utr3[,c("utrEnd","utrStart")] <- c(paste(utr3_candidate_end,collapse = ','),paste(utr3_candidate_start,collapse=','))
        utr3$start_w_transcript <- min(utr3_candidate_start) == utr3$start
        utr3$method <- 'coding_sequence'
      } else {
        utr3_candidate_end <- min_cdsstarts - 1
        if (min(utrstarts) < utr3_candidate_end) {
          utr3_candidate_start <- min(utrstarts)
          method <- 'substract_utr_by_cds'
        } else {
          utr3_candidate_start <- utr3$start
          method <- 'substract_transcript_by_cds'
        }
        utr3[,c("utrEnd","utrStart")] <- c(paste(utr3_candidate_end,collapse = ','),paste(utr3_candidate_start,collapse=','))
        utr3$start_w_transcript <- utr3_candidate_start == utr3$start
        utr3$method <- method
      } 
    }
    return(utr3)
  })
  
  utr3_coding_m_df <- do.call('rbind',utr3_coding_m)
  cat('Isoforms with UTR3 regions computed using stop_codon method:',length(which(utr3_coding_m_df$method == 'stop_codon')),'\n')
  cat('Isoforms with UTR3 regions computed using coding_sequence method:',length(which(utr3_coding_m_df$method == 'coding_sequence')),'\n')
  cat('Isoforms with UTR3 regions computed using substract_transcript_by_cds method:',length(which(utr3_coding_m_df$method == 'substract_transcript_by_cds')),'\n')
  
  utr3_coding_df <- rbind(utr3_coding_p_df[,c('isoform.id','utrStart','utrEnd')],
                          utr3_coding_m_df[,c('isoform.id','utrStart','utrEnd')])
  
  genes_df[,c('utrStart','utrEnd')] <- utr3_coding_df[match(genes_df$isoform.id,utr3_coding_df$isoform.id),c('utrStart','utrEnd')]
  genes_df <- genes_df %>% dplyr::select(-c("num_startcodon","num_stopcodon","startcodon_start","startcodon_end","stopcodon_start","stopcodon_end"))
  
  return(genes_df)
}

separate_consecutive <- function(vec) {
  numeric_vec <- gsub('TU','',vec) %>% as.numeric
  # Find the differences between consecutive elements
  diffs <- c(0, diff(numeric_vec))
  
  # Identify consecutive elements and non-consecutive elements
  consecutive <- split(vec, cumsum(diffs != 1))
  non_consecutive <- split(vec, cumsum(diffs == 1))
  
  return(consecutive)
}

merge_mTU <- function(tus_sorted,tu,index=NA) {
  gene.id <- tu$gene.id[which(tu$tu %in% tus_sorted)] %>% unique
  gene.symbol <- tu$name[which(tu$tu %in% tus_sorted)] %>% unique
  rep_tu <- tus_sorted[1]
  rep_gene.id <- paste0(gene.id,if(is.na(index)==FALSE){paste0('.',index)})
  rep_gene.symbol <- paste0(gene.symbol,if(is.na(index)==FALSE){paste0('.',index)})
  rep_strand <- tu$strand[which(tu$tu == rep_tu)]
  rep_chr <- tu$chr[which(tu$tu == rep_tu)]
  last_tu <- tus_sorted[length(tus_sorted)]
  new_tu <- c(rep_chr,tu$start[which(tu$tu == rep_tu)],tu$end[which(tu$tu == last_tu)],
              rep_strand,rep_gene.id,rep_gene.symbol,rep_tu,paste0(rep_tu,':',rep_gene.symbol))
  tu[which(tu$tu == rep_tu),] <- new_tu
  tu <- tu[!(tu$tu %in% tus_sorted[-1]),]
  
  return(tu)
}

rename_mTU <- function(mtus,tu) {
  tus <- unique(mtus)
  tus <- tus[mixedorder(tus)]
  gene.id <- tu$gene.id[which(tu$tu %in% tus)] %>% unique
  if (length(gene.id)!=1) {stop(paste0(tus,collapse=', '),' do not share the same gene ID')}
  consecutive_tus <- separate_consecutive(tus)
  if (length(consecutive_tus) > 1) {
    for (i in 1:length(consecutive_tus)) {
      if (length(consecutive_tus[[i]]) > 1) {
        tu <- merge_mTU(consecutive_tus[[i]],tu,index=i)
      } else {
        tu[which(tu$tu == consecutive_tus[[i]]),c('gene.id','name','tu_anno')] <- apply(tu[which(tu$tu == consecutive_tus[[i]]),c('gene.id','name','tu_anno')],1,
                                                                                        function(x) paste0(x,'.',i)) 
      }
    }
  } else {
    tu <- merge_mTU(consecutive_tus[[1]],tu,index=NA)
  }
  return(tu)
}

create_TU_from_gtf <- function(genes,flank=5000,outdir='./',use.genesymbol=FALSE,save=TRUE,
                               final_tu_annot='gene.id',tu.create.only=FALSE,
                               chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                                        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                                        'chr20','chr21','chr22','chrX','chrY')){
  cat('Use',if(use.genesymbol){paste0('gene symbol')}else{paste0('Ensembl ID')},'to collapse transcripts by \n')
  if (use.genesymbol) {
    genes$gene.id <- genes$name
  }
  # Rename genes on different chromosomes and strands with the same symbols
  strand_tally <- genes %>% group_by(gene.id) %>% summarise(nStrand = length(unique(strand)))
  genes_mStrand <- strand_tally$gene.id[which(strand_tally$nStrand > 1)]
  if (length(genes_mStrand) != 0) {
    cat('There are',length(genes_mStrand),'genes on different strands with the same symbol. Rename them \n')
    genes[which(genes$gene.id %in% genes_mStrand),c('strand','name','gene.id')] <- apply(genes[which(genes$gene.id %in% genes_mStrand),c('strand','name','gene.id')],1,function(x){
      if (x[1] == '+') {new_name <- paste0(x[-1],'.1')} else if (x[1] == '-') {new_name <- paste0(x[-1],'.2')}
      return(c(x[1],new_name))
    }) %>% t()
  }
  chr_tally <- genes %>% group_by(gene.id) %>% summarise(nChr = length(unique(chr)))
  genes_mChr <- chr_tally$gene.id[which(chr_tally$nChr > 1)]
  if (length(genes_mChr) != 0) {
    cat('There are',length(genes_mChr),'genes on multiple chromosomes with the same symbol. Rename them \n')
    for (g in genes_mChr) {
      mChr <- genes$chr[which(genes$gene.id == g)] %>% unique
      for (i in 1:length(mChr)) {
        genes[which(genes$gene.id == g & genes$chr == mChr[i]),c('name','gene.id')] <- apply(genes[which(genes$gene.id == g & genes$chr == mChr[i]),c('name','gene.id')],
                                                                                             1, function(x) paste0(x,'.',i))
      }
    }
  }
  
  genes <- makeGRanges(genes,strand=T)
  stopifnot(genes$strand %in% c("+","-"))
  
  # Some genes may have disjoint txUnits, let's treat them as separate genes to make things easier
  cat("Clustering TUs \n")
  df <- as.data.frame(table(genes$gene.id))
  genes.singles <- genes[genes$gene.id %in% df[df$Freq==1,]$Var1] 
  genes.tored <- genes[genes$gene.id %in% df[df$Freq>1,]$Var1] 
  spl <- split(genes.tored,genes.tored$gene.id)
  spl.red <- lapply(spl,function(x) reduce(x,with.revmap=FALSE)) #collect union boundary of the same gene, i.e. the discrete range(s) that together cover all the isoforms' starts and ends. Ninh: Drop everything but seqnames. ranges and strand
  for(i in 1:length(spl.red))
  {
    spl.red[[i]]$gene.id <- names(spl)[i]
  }
  
  genes.red <- do.call(c,unname(spl.red)) # Merge all the elements in a list. This object contains genes with multiple isoforms that have been reduced into union boundaries
  genes.red$name <- genes[match(genes.red$gene.id,genes$gene.id)]$name # Transfer gene symbols
  cat("Number of genes with single isoform annotated:",length(unique(genes.singles$gene.id)),'\n')
  cat("Number of genes with multiple isoform annotated:",length(unique(genes.tored$gene.id)),'\n')
  dt <- makeDT(genes.red) #Convert a GRanges object into a dataframe
  genes.red.counts <- dt[,list(nPos=length(strand),nStrands=length(unique(strand))),by="gene.id"]
  
  # Assign isoforms to TUs
  tu <- genes.singles
  values(tu) <- NULL # removes all metadata
  tu$gene.id <- genes.singles$gene.id
  tu$name <- genes.singles$name
  tu <- c(tu,genes.red)
  seqlevels(tu) <- chrs 
  tu <- tu[order(tu)] #order orders by strand, seqnames and ranges
  tu$tu <- paste0("TU",1:length(tu)) # So each TU corresponds to one isoform. If a gene has multiple isoforms, such gene can have multiple TU, each corresponding to one discrete range
  if (final_tu_annot == 'gene_symbol') {
    cat('Use gene symbol for final TU annotation \n')
    tu$tu_anno<-paste0(tu$tu,':',tu$name) 
  } else if (final_tu_annot == 'gene.id') {
    cat('Use whatever ID in gene.id column for final TU annotation. In this run, it is',if(use.genesymbol){paste0('gene symbol')}else{paste0('Ensembl ID')},'\n')
    tu$tu_anno<-paste0(tu$tu,':',tu$gene.id)
  } else {stop('Acceptable values are either gene_symbol or gene.id, with the latter meaning use whatever ID in gene.id column for final TU annotation.')}
  
  
  # From now on, must aggregate by TU rather than gene.id
  cat("Assigning transcripts to TUs \n")
  
  genes.singles$tu <- tu[match(genes.singles$gene.id,tu$gene.id),]$tu
  genes.singles$tu_anno <- tu[match(genes.singles$gene.id,tu$gene.id),]$tu_anno
  
  fo <- data.table(as.data.frame(findOverlaps(genes.tored,tu)))#to see how multiple-isoform genes are overlapped with the tu-assigned set; as a default, strand information is also included in the criteria
  fo$gene_gid <- genes.tored[fo$queryHits]$gene.id
  fo$tu_gid <- tu[fo$subjectHits]$gene.id
  fo$assign <- fo$gene_gid==fo$tu_gid
  fo_F <- fo[assign==FALSE,]
  write.table(fo_F, paste0(outdir,'TUs_overlap_with_othergenes.csv'),sep='\t',row.names=FALSE,quote=FALSE)
  fo <- fo[assign==TRUE,] #For the one w/ FALSE, a different gene is overlaped. The action with assign==TRUE still retains the one w/ overlapped with the other gene. An isoform despite belonging to one gene (its TU) can overlap with TU from a different gene. The gene/TU it belongs to still shows up as TRUE so we don't lose any isoform
  stopifnot(!duplicated(fo$queryHits)) #to make sure queryHits, i.e. queries from genes.tored, i.e, all isoforms from multiple-isoform genes are not duplicated in findOverlaps
  stopifnot(nrow(fo)==length(genes.tored)) #to make sure all multiple-isoform genes are present in findOverlaps
  dt <- makeDT(genes.tored)
  dt[fo$queryHits,tu:= tu[fo$subjectHits]$tu] #bring tu id.  Ninh: Assign TU that isoforms from multiple-isoform genes overlap with
  #dt[,tu_anno:= paste0(dt$tu,':',dt$gene.id)] #Add tu anno
  dt$tu_anno <- tu[match(dt$tu,tu$tu)]$tu_anno
  genes.tu <- c(genes.singles,makeGRanges(dt,strand=TRUE))
  
  stopifnot(length(unique(genes.tu$tu))==length(tu))
  stopifnot(length(genes.tu)==length(genes))
  stopifnot(length(unique(genes$gene.id))==length(unique(genes.tu$gene.id)))
  seqlevels(genes.tu) <- chrs
  genes.tu <- genes.tu[order(genes.tu)] #Osbervation: TU is numbered based on a sorting by strand, chromosome, and coordinates
  
  cat("Number of genes with isoforms reduced to disjoint clusters/TUs on same strand only:",nrow(genes.red.counts[(nPos>1)&(nStrands==1),]),'\n') #The number of genes with more than 1 union boundaries all on the same strand
  # Reduced to disjoint clusters on same strand only: 
  cat("Number of genes with isoforms reduced to disjoint clusters/TUs on different strands:",nrow(genes.red.counts[(nStrands>1),]),'\n') #The number of genes with more than 1 union boundaries on the different strands
  # Reduced to disjoint clusters on different strands: 0 # Should be empty as we fix it earlier
  if(nrow(genes.red.counts[(nPos>1)&(nStrands==1),]) != 0) {
    # Examine genes with more than one TU
    genes_with_multiple_TU <- genes.tu[which(genes.tu$gene.id %in% genes.red.counts[(nPos>1)&(nStrands==1),]$gene.id)]
    genes_with_multiple_TU <- makeDT(genes_with_multiple_TU)
    write.table(genes_with_multiple_TU, paste0(outdir,'genes_with_multiple_TU.csv'),sep='\t',quote=FALSE)
    cat("For genes with multiple TUs, treat each set of consecutive TUs as one separate genes. Check genes_with_multiple_TU.csv. \n")
    genes_mTU_spl <- split(genes_with_multiple_TU, genes_with_multiple_TU$gene.id)
    
    tu <- as.data.frame(tu)
    colnames(tu)[1] <- 'chr'
    tu$chr <- as.character(tu$chr)
    tu$strand <- as.character(tu$strand)
    tu$width <- NULL
    
    for (i in 1:length(genes_mTU_spl)) {tu <- rename_mTU(genes_mTU_spl[[i]]$tu,tu)}
    tu$start <- as.numeric(tu$start)
    tu$end <- as.numeric(tu$end)
    
    cat("Fix annotations in isoform TU assignment as well. \n")
    genes.tu <- as.data.frame(genes.tu)
    colnames(genes.tu)[1] <- 'chr'
    genes.tu$chr <- as.character(genes.tu$chr)
    genes.tu$strand <- as.character(genes.tu$strand)
    genes.tu$width <- NULL
    for (i in 1:length(genes_mTU_spl)) {
      mtus <- unique(genes_mTU_spl[[i]]$tu)
      mtus <- mtus[mixedorder(mtus)]
      consecutive_tus <- separate_consecutive(mtus)
      for (i in 1:length(consecutive_tus)) {
        rep_tu <- intersect(consecutive_tus[[i]],tu$tu)
        rep_tu_anno <- as.list(tu[which(tu$tu == rep_tu),c('gene.id','name','tu','tu_anno')])
        genes.tu[which(genes.tu$tu %in% consecutive_tus[[i]]),c('gene.id','name','tu','tu_anno')] <- rep_tu_anno
      }
    }
    
    nTU_pergene <- tu %>% group_by(gene.id) %>% summarise(nTU=n())
    cat("Number of genes with more than 1 TU:",nrow(nTU_pergene[which(nTU_pergene$nTU>1),]),'\n')
    stopifnot(length(setdiff(genes.tu$tu,tu$tu)) == 0)
    
    tu <- makeGRanges(tu,strand=T)
    genes.tu <- makeGRanges(genes.tu,strand=T)
  }
  if (tu.create.only) {
    cat('Do not annotate TU \n')
    ret <- list(tu = tu, genes.tu = genes.tu)
    if(save) {saveRDS(ret,paste0(outdir,'genes.rds'))}
    return(ret)
  }
  
  # Assign TUs as coding
  cat("Detecting coding TUs. A TU will be counted as coding if it contain at least one coding isoform.\n")
  genes.tu$coding <- !(genes.tu$cdsStart %>% is.na) 
  dt <- makeDT(genes.tu)
  dt <- dt[,list(coding=any(coding)),by="tu"] #Check if any coding isoform falls under each tu. If a tu has one coding isoform, count it as TRUE for coding
  stopifnot(nrow(dt)==length(tu))
  tu$coding <- "new"
  tu[match(dt$tu,tu$tu)]$coding <- dt$coding
  
  # Get 3' end flanks (i.e., 5k bp segment right after TSE. Leave it be, will be modified with the peaks)
  cat("Generating TU flanks of ",flank," bp on the 3' end \n")
  tu.flank <- flank(tu,width=flank,start=FALSE,both=FALSE,ignore.strand=FALSE) # after on forwards strand and before on reverse strand
  
  # Get 3' UTRs
  cat("Clustering 3' UTRs for each coding TU \n")
  coding <- genes.tu[genes.tu$coding==TRUE] %>% as.data.frame() # Get the coding isoform
  wutr3 <- coding[which(is.na(coding$utrStart) == FALSE),] # Coding isoforms with UTR
  
  spl <- split(wutr3,wutr3$tu)
  spl.utr3 <- lapply(spl, function(x) {
    utr3_start <- paste(x$utrStart, collapse = ',') %>% strsplit(split =',') %>% unlist
    utr3_end <- paste(x$utrEnd, collapse = ',') %>% strsplit(split =',') %>% unlist
    stopifnot(length(utr3_start) == length(utr3_end))
    utr3.ranges <- matrix(nrow = length(utr3_start), ncol = 7) %>% as.data.frame
    colnames(utr3.ranges) <- c('seqnames','start','end','strand','gene.id','name','tu')
    utr3.ranges$seqnames <- unique(x$seqnames %>% as.character)
    utr3.ranges$strand <- unique(x$strand %>% as.character)
    utr3.ranges$gene.id <- unique(x$gene.id %>% as.character)
    utr3.ranges$name <- unique(x$name %>% as.character)
    utr3.ranges$tu <- unique(x$tu %>% as.character)
    utr3.ranges$start <- utr3_start
    utr3.ranges$end <- utr3_end
    return(utr3.ranges)
  })
  spl.utr3 <- do.call('rbind', spl.utr3)
  utr3.all <- GRanges(spl.utr3$seqnames,IRanges(spl.utr3$start %>% as.numeric, spl.utr3$end %>% as.numeric),
                      strand=spl.utr3$strand,tu=spl.utr3$tu,gene.id=spl.utr3$gene.id,name=spl.utr3$name)
  
  redme <- function(myutrs)
  {
    df <- as.data.frame(table(myutrs$tu))
    myutrs.singles <- myutrs[myutrs$tu %in% df[df$Freq==1,]$Var1]
    myutrs.tored <- myutrs[myutrs$tu %in% df[df$Freq>1,]$Var1]
    spl <- split(myutrs.tored,myutrs.tored$tu)
    spl.red <- lapply(spl,reduce)
    for(i in 1:length(spl.red))
    {
      spl.red[[i]]$tu <- names(spl)[i]
    }
    myutrs.red <- do.call(c,unname(spl.red))
    myutrs.red$gene.id <- myutrs[match(myutrs.red$tu,myutrs$tu)]$gene.id
    myutrs.red$name <- myutrs[match(myutrs.red$tu,myutrs$tu)]$name
    c(myutrs.singles,myutrs.red)
  }
  
  utr3 <- redme(utr3.all)
  utr3 <- utr3[order(utr3)]
  utr3 <- utr3[width(utr3)>0] #Filter out UTR regions that start right after end of isoforms, i.e. no UTR. In this case, all UTRs retained
  utr3$utr <- paste0("UTR",1:length(utr3))
  
  # Internal
  cat("Computing internal regions of coding isoforms via setdiff. Internal regions are within a TU but outside its UTR3s \n")
  spl <- split(utr3,utr3$tu) #Split UTR3 by TU
  coding <- tu[tu$coding==TRUE]
  
  # Retrieve coding TUs that don't have 3' UTRs
  coding <- coding[coding$tu %in% names(spl)]
  diff <- lapply(coding$tu,function(x) IRanges::setdiff(coding[coding$tu==x],spl[[x]])) #Obtain the regions within each TU that are outside its UTRs 
  for(i in 1:length(diff))
  {
    diff[[i]]$tu <- coding$tu[i]
  }
  
  internal <- do.call(c,diff)
  internal$gene.id <- tu[match(internal$tu,tu$tu)]$gene.id
  internal$name <- tu[match(internal$tu,tu$tu)]$name
  
  # No-UTR gene bodies
  cat("Collects coding TUs with no UTR in noutr \n")
  coding <- tu[tu$coding==TRUE]
  noutr <- coding[!(coding$tu %in% names(spl))] #because TUs with no UTR have been filtered out when generating splt
  
  # Between annotated 3' UTRs
  cat("Computing between multi-UTR. Between here refers to regions lying in-between UTRs \n")
  dt <- makeDT(utr3)
  dt <- dt[,list(chr=chr[1],start=min(start),end=max(end),strand=strand[1],gene.id=gene.id[1],name=name[1],nUtr=length(utr)),by="tu"] #The range that encompasses all UTRs within a TU
  multis <- dt[nUtr>1,]
  # Now have the maximal range for all the multis
  # If we setdiff out the real 3', then we have just the between ranges
  union <- makeGRanges(multis,strand=T)
  spl <- split(utr3,utr3$tu)
  diff <- lapply(union$tu,function(x) IRanges::setdiff(union[union$tu==x],spl[[x]])) #Regions within the max UTR that did not overlap with the component UTRs
  for(i in 1:length(diff))
  {
    diff[[i]]$tu <- union$tu[i]
  }
  
  btw <- do.call(c,diff)
  btw$gene.id <- tu[match(btw$tu,tu$tu)]$gene.id
  btw$name <- tu[match(btw$tu,tu$tu)]$name
  
  # Now remove the betweens from the internals
  cat("Now remove the betweens from internals \n")
  internal <- internal[!(internal %in% btw)] # The regions that outsides UTRs but still lie within the max UTR region don't count as internal
  
  ret <- list(tu=tu,genes.tu=genes.tu,flank=tu.flank,utr3=utr3,between=btw,internal=internal,noutr3=noutr) # Summarize the steps that lead to tu, genes.tu, tu.flank, utr3, btw, internal and noutr
  
  if(save) {saveRDS(ret,paste0(outdir,'genes.rds'))}
  
  return(ret)
}

update_flank <- function(peaks_annot,genes, save=TRUE, 
                         chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                                  'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                                  'chr20','chr21','chr22','chrX','chrY')) {
  cat('Modify the flank regions so that it ends at the most downstream peak. If no peaks are found within the flank regions and the no peak overfloods the TU border, the regions are deleted.\n')
  tu <- as.data.frame(genes$tu)
  peaks_annot$tu_start <- tu$start[match(peaks_annot$tu,tu$tu)]
  peaks_annot$tu_end <- tu$end[match(peaks_annot$tu,tu$tu)]
  peaks_annot$tu_need_update[(peaks_annot$strand == '+' & (peaks_annot$end > peaks_annot$tu_end))] <- TRUE
  peaks_annot$tu_need_update[(peaks_annot$strand == '-' & (peaks_annot$start < peaks_annot$tu_start))] <- TRUE
  peaks_annot$tu_need_update[is.na(peaks_annot$tu_need_update)] <- FALSE
  
  flank_peaks <- peaks_annot[which(peaks_annot$tu_need_update),]
  flank_peaks <- split(flank_peaks, flank_peaks$tu)
  flank_peaks_lastpeakend <- lapply(flank_peaks, function(x) {
    if (x$strand %>% unique == '-') {
      lastpeakend <- min(x$start)
    } else {
      lastpeakend <- max(x$end)
    }
    return(lastpeakend)
  })
  tu_flanks <- makeDT(genes$flank)
  tu_coords <- makeDT(genes$tu)
  tu_flanks_updated <- lapply(split(tu_flanks,tu_flanks$tu), function (x) {
    tu_id <- x$tu
    if (tu_id %in% names(flank_peaks_lastpeakend)) {
      updated_3_prime_end <- flank_peaks_lastpeakend[[tu_id]] %>% as.numeric
      if (x$strand == '-') {x$start <- updated_3_prime_end} else {
        x$end <- updated_3_prime_end
      }
    } else {
      if (x$strand == '-') {
        tu_start <- tu_coords[which(tu_coords$tu == tu_id),'start']
        if (x$end != tu_start) {
          x$end <- x$end + 1
          x$start <- x$end
        }
      } else {
        tu_end <- tu_coords[which(tu_coords$tu == tu_id),'end']
        if (x$start != tu_end) {
          x$start <- x$start - 1
          x$end <- x$start
        }
      }
    }
    x$width <- x$end - x$start
    return(x)
  })
  tu_flanks_updated <- do.call('rbind',tu_flanks_updated)
  tu_flanks_updated <- tu_flanks_updated[mixedorder(tu_flanks_updated$tu),]
  tu_flanks_updated <- makeGRanges(tu_flanks_updated,strand=TRUE)
  seqlevels(tu_flanks_updated) <- chrs
  genes$flank <- tu_flanks_updated
  
  if(save) {saveRDS(genes,paste0(outdir,'genes_flankupdated.rds'))}
  
  return(genes)
}

peak_ext <- function(peaks, ext_len=100, direction='5') {
  plus <- peaks[which(peaks$strand == '+'),]
  minus <- peaks[which(peaks$strand == '-'),]
  if (direction=='5') {
    plus$start <- plus$start - ext_len
    minus$end <- minus$end + ext_len
  } else if (direction=='3') {
    plus$start <- plus$start + ext_len
    minus$end <- minus$end - ext_len
  } else (stop('Direction arg needs to be either 3 or 5'))
  peaks <- rbind(plus,minus)
  peaks$peak_width <- peaks$end - peaks$start
  return(peaks)
}