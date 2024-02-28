### New Functions ###

# Flag features expressed by < percent of cells. Currently using counts as a proxy for that filtering
stratify_matrix <- function(merged_counts,meta,vars='seurat_clusters',cutoff_pct=10,min_cell_per_group=10) {
  
  vars_idx <- which(colnames(meta) %in% vars)
  for (i in 1:length(vars_idx)) {if (i == 1) {group <- meta[,vars_idx[i]]} else {group <- paste0(group,'_',meta[,vars_idx[i]])}}
  meta$group <- group
  meta_spl <- split(meta,meta$group)
  
  min_cell_exp_rows<-function(peak_mat,cutoff=cutoff_pct)
  {
    n<-ncol(peak_mat)
    rsum<-rowSums(peak_mat)
    min_cell<-((cutoff/100)*n) %>% round()
    r<-rownames(peak_mat)[which(rsum>=min_cell)]
    return(r)
  }
  
  cat('Filter out peaks/genes whose total counts across cells <=',cutoff_pct,'% number of cells \n')
  
  merged_counts_spl <- lapply(meta_spl, function(x) {
    group_name <- unique(x$group)
    if (nrow(x) < min_cell_per_group) {
      mtx_filtered <- NA
      cat('Remove',group_name,'because it has fewer than',min_cell_per_group,'cells \n')
    } else {
      mtx <- merged_counts[,row.names(x)]
      peaks_to_keep <- min_cell_exp_rows(mtx,cutoff=cutoff_pct)
      mtx_filtered <- mtx[peaks_to_keep,]
    }
    return(mtx_filtered)
  })
  
  return(merged_counts_spl) 
}

APA_feature_filter <- function(peak_mat,min_cell_expr_pct=10,expr.thres=1) {
  pct.cell.expr.features <- apply(peak_mat,1,function(x) {
    pct.cell.expr <- length(which(x >= expr.thres))/length(x) * 100
    return(pct.cell.expr)
  })
  features.flag <- pct.cell.expr.features < min_cell_expr_pct
  df <- data.frame(pct.expr = pct.cell.expr.features, flag = features.flag)
  df$TU <- sapply(str_split(row.names(df),pattern = ':', n = 3), `[`, 1)
  tu_all_flag <- df %>% group_by(TU) %>% summarise(low.expr = all(flag)) %>% as.data.frame()
  tu.rm <- tu_all_flag[which(tu_all_flag$low.expr == TRUE),'TU']
  tu.left <- tu_all_flag[which(tu_all_flag$low.expr == FALSE),'TU']
  peak.rm <- df[which(df$TU %in% tu.rm),] %>% row.names()
  peak.left <- df[-which(df$TU %in% tu.rm),] %>% row.names()
  cat('Remove',length(tu.rm),'TUs, consisting of',length(peak.rm),'peaks, for having all of their peaks expressed by fewer than',min_cell_expr_pct,'% of cells \n')
  peak_mat_flt <- peak_mat[peak.left,]
  cat('There are',length(tu.left),'TUs, consisting of',nrow(peak_mat_flt),'peaks left \n')
  low.expr.peak.left <- df[which(row.names(df) %in% peak.left & df$flag == TRUE),] %>% row.names()
  cat('Out of them,',length(low.expr.peak.left),'peaks are expressed by fewer than',min_cell_expr_pct,'% of cells \n')
  peak_mat_flt[low.expr.peak.left,] <- NA
  return(peak_mat_flt)
}

create_replicates<-function(clus_name,peak_mat,nrep=2,p=0.7,method='random',num.vars=2000,
                            elbow=TRUE,numpc=20,pcs=seq(1,20),save.clus=TRUE,outdir='./',fprefix,seed=123){
  bc<-colnames(peak_mat)
  ncell<-length(bc)
  
  rep<-list()
  
  if (method == 'random') {
    for(i in 1:nrep)
    {
      set.seed(i)
      subset_bc<-sample(bc,(p*ncell),replace = FALSE)
      subset_mat<-peak_mat[,subset_bc]
      rep[[i]]<-rowSums(subset_mat) %>% as.data.frame()
      names(rep)[[i]]<-paste0(clus_name,'_',i)
    }
  } else {stop('Acceptable value is: random (Sampling cell barcodes without replacement)')}
  
  return(rep)
}

create_test_inputs <- function(test='APA',all_groups,ident1,ident2,min_cell_expr_pct=10,expr.thres=1,pseudo.count=1,APA.feature.filter=TRUE,
                               replicate='random',meta,nrep=3,p=0.7,min_cell_per_rep=10,num.vars=2000,elbow=TRUE,
                               numpc=20,pcs=seq(1,20),outdir='./',save.clus=TRUE) {
  cbind_intersect_peaks <- function(groups) {
    peaks <- lapply(groups, function(x) {return(row.names(x))})
    common_peaks <- Reduce(intersect,peaks)
    groups_intersect <- lapply(groups, function(x) {return(x[common_peaks,])})
    
    names(groups_intersect) <- NULL
    groups_intersect <- do.call('cbind',groups_intersect)
    
    return(groups_intersect)
  }
  
  #####
  ident1_group <- all_groups[which(names(all_groups) %in% ident1)]
  ident1_group_name <- paste(ident1,collapse = '_')
  
  if (paste(ident2,collapse='_')=='rest') {
    ident2_group <- all_groups[-which(names(all_groups) %in% ident1)]
    ident2_group_name <- 'rest'
  } else {
    ident2_group <- all_groups[which(names(all_groups) %in% ident2)]
    ident2_group_name <- paste(ident2,collapse = '_')
  }
  
  if (replicate %in% c('random','sample')) {
    if (test == 'APA') {
      ident1_group <- do.call('cbind',ident1_group)
      ident2_group <- do.call('cbind',ident2_group)
      if (APA.feature.filter) {
        ident1_group <- APA_feature_filter(ident1_group,min_cell_expr_pct=min_cell_expr_pct,expr.thres=expr.thres)
        ident2_group <- APA_feature_filter(ident2_group,min_cell_expr_pct=min_cell_expr_pct,expr.thres=expr.thres)
      }
    } else if (test == 'DEG') {
      ident1_group <- cbind_intersect_peaks(ident1_group)
      ident2_group <- cbind_intersect_peaks(ident2_group)
    } else {stop('Acceptable values are: APA (Preparing inputs for APA test), DEG (Preparing inputs for DEG test)')}
  } else {stop('Acceptable values are: sample (Treat each sample as a replicate), random (Sampling cell barcodes without replacement using create_replicates)')}
  
  
  if (replicate == 'random') {
    cat('Create pseudoreplicates by random sampling...\n')
    
    ident1_group_reps <- create_replicates(ident1_group_name,ident1_group,nrep=nrep,p=p)
    ident1_group_reps <- do.call('cbind',ident1_group_reps)
    colnames(ident1_group_reps) <- paste0(ident1_group_name,'_rep',1:nrep)
    
    ident2_group_reps <- create_replicates(ident2_group_name,ident2_group,nrep=nrep,p=p)
    ident2_group_reps <- do.call('cbind',ident2_group_reps)
    colnames(ident2_group_reps) <- paste0(ident2_group_name,'_rep',1:nrep)
  } else if (replicate == 'sample') {
    cat('Split cluster by original sample...\n')
    
    cells.of.interest <- c(colnames(ident1_group),colnames(ident2_group))
    meta <- meta[which(row.names(meta) %in% cells.of.interest),]
    meta_spl <- split(meta,meta$sample)
    
    cat('Remove replicates in each group with fewer than',min_cell_per_rep,'cells \n')
    
    ident1_group_reps <- lapply(meta_spl,function(x) {
      cells.subset <- intersect(colnames(ident1_group),row.names(x))
      if (length(cells.subset) >= min_cell_per_rep) {
        sample_rep <- ident1_group[,cells.subset]
        sample_rep <- rowSums(sample_rep)
        return(sample_rep)
      }
    })
    ident1_group_reps <- do.call('cbind',ident1_group_reps)
    colnames(ident1_group_reps) <- paste0(ident1_group_name,'_rep',1:ncol(ident1_group_reps))
    
    ident2_group_reps <- lapply(meta_spl,function(x) {
      cells.subset <- intersect(colnames(ident2_group),row.names(x))
      if (length(cells.subset) >= min_cell_per_rep) {
        sample_rep <- ident2_group[,cells.subset]
        sample_rep <- rowSums(sample_rep)
        return(sample_rep)
      }
    })
    ident2_group_reps <- do.call('cbind',ident2_group_reps)
    colnames(ident2_group_reps) <- paste0(ident2_group_name,'_rep',1:ncol(ident2_group_reps))
  }
  
  common_peaks <- intersect(row.names(ident1_group_reps),row.names(ident2_group_reps))
  
  ident1_group_reps <- ident1_group_reps[common_peaks,]
  ident2_group_reps <- ident2_group_reps[common_peaks,]
  
  test_input <- cbind(ident1_group_reps,ident2_group_reps)
  
  rows_with_na <- apply(is.na(test_input),1,any)
  test_input[rows_with_na,] <- replace(test_input[rows_with_na,],is.na(test_input[rows_with_na,]),pseudo.count)
  
  return(test_input)
}

APA_DEXseq_test <- function(ident1,ident2,peak_mat,vars=NA,adjust.var=NULL,min_peak=2,ncpu=4,dispersion.plot.save=TRUE,peak_ref=peak_ref,outdir='./') {
  cat('DEXseq testing:',ident1,'_v_',ident2,'\n')
  sd <- data.frame(sample = colnames(peak_mat),
                   group = gsub('_rep\\d+$','',colnames(peak_mat)))
  sd$sample <- factor(sd$sample, levels = colnames(peak_mat))
  sd$group <- factor(sd$group, levels = c(ident2,ident1))
  sd <- arrange(sd,sample)
  
  if (is.na(vars)==FALSE) {
    sd <- cbind(sd,str_split_fixed(sd$group, '_', length(vars)))
    colnames(sd)[3:(3+length(vars)-1)] <- vars
  }
  
  row.names(sd) <- sd$sample
  sd$sample <- NULL
  
  cat('Filter peak reference...\n')
  
  # Filter peaks based on read count within samples to be tested #. No need cause already filter with min_cell_exp_rows
  # Retain only those peaks in this new peak ref that made the read cutoff (default 10) as well as assigned to unique TU (Transcription Unit) #
  # This might be misleading because unique_tu status was determined before unique_peak == FALSE was filtered. So remove this filtering
  # Filter for peaks in the comparison
  peak_ref_use <- peak_ref[which(peak_ref$peak %in% row.names(peak_mat)),]
  
  # Sanity checks #
  stopifnot(!str_detect(peak_ref_use$over_tus,","))
  stopifnot(!str_detect(peak_ref_use$flank_tus,","))
  stopifnot(all(peak_ref_use$unique_peak)) 
  stopifnot(length(unique(peak_ref_use$peak))==nrow(peak_ref_use))
  stopifnot(class(peak_ref_use$peak)=="character")
  stopifnot(peak_ref_use$peak %in% rownames(peak_mat))
  
  # Create data table for this new peak ref #
  peak_ref_use  <- data.table(peak_ref_use, key="tu")
  
  # Add relevant columns by TU #
  peak_ref_use <- peak_ref_use[,list(peak=peak,
                                     type=type,
                                     coding=coding,
                                     unique_peak=unique_peak,
                                     unique_tu=unique_tu,
                                     over_tus=over_tus,
                                     flank_tus=flank_tus,
                                     num_peak=length(peak),
                                     gene=gene,
                                     final_classification=final_classification,
                                     merged_peak=merged_peak),
                               by="tu"]
  
  # Remove TU with 1 peak remaining after applying all filters #
  cat('Remove from peak reference TUs with fewer than',min_peak,'peaks found in this comparison \n')
  peak_ref_use <- peak_ref_use[num_peak >= min_peak,] 
  # Create genomic range from the remaining peaks
  peak<-peak_ref[match(peak_ref_use$peak,peak_ref$peak),c('chr','start','end','strand','peak')]
  rownames(peak) <- 1:nrow(peak) 
  peak_gr<-makeGRangesFromDataFrame(peak,keep.extra.columns = TRUE)
  
  exon_num <- peak_ref_use[,list(paste0("E",1:length(peak)),peak=peak),by="tu"]
  stopifnot(exon_num$tu==peak_ref_use$tu)
  stopifnot(exon_num$peak==peak_ref_use$peak)
  peak_ref_use$exon_num <- exon_num$V1
  peak_ref_use$key <- paste(peak_ref_use$tu,peak_ref_use$exon_num,sep=":")
  
  cat("Building DEXSeqDataSet object \n")
  if(is.null(adjust.var))
  {
    design <- "~ sample + exon + group:exon"
  }else
  {
    design <- paste0("~ sample + exon + adjust:exon + group:exon")
    sd$adjust <- factor(with(sd,get(adjust.var)))
  }
  cat("Design Formula: ",design,'\n')
  
  # Now subset peak in counts and mm to include only the peaks from this new peak ref #
  cat('Subset inputs according to filtered peak reference \n')
  counts<-peak_mat[peak_ref_use$peak,]
  
  stopifnot(rownames(counts)==peak_ref_use$peak)
  stopifnot(peak_gr$peak==peak_ref_use$peak)
  
  cat('Create DEXSeq object \n')
  dxd <- DEXSeqDataSet(countData=counts, sampleData=sd, design=as.formula(design),
                       featureID=peak_ref_use$exon_num, groupID=peak_ref_use$tu, featureRanges=peak_gr)
  
  # Run tests
  cat("Estimating size factors \n")
  dxd <- DEXSeq::estimateSizeFactors(dxd)
  cat("Estimating dispersions \n")
  dxd <- DEXSeq::estimateDispersions(dxd)
  
  if (dispersion.plot.save) {
    png(paste0(outdir,ident1,'_v_',ident2,'_DispEsts.png'),width = 6, height = 6, units = 'in', res = 1200)
    plotDispEsts(dxd)
    dev.off()
  }
  
  if(is.null(adjust.var))
  {
    cat("Tesing w/o adjustment \n")
    dxd<-testForDEU(dxd)
  }else
  {
    cat("Tesing w/ adjustment for variable: ", adjust.var,'\n')
    dxd <- testForDEU(dxd,fullModel=as.formula("~ sample + exon + adjust:exon + group:exon"),reducedModel=as.formula("~ sample + exon + adjust:exon"))
  }
  
  dxd <- estimateExonFoldChanges(dxd,fitExpToVar="group",denominator=ident2,BPPARAM=MulticoreParam(workers=ncpu))
  
  # Pull APA results table
  cat("Extracting DEXSeqResults \n")
  dxr <- DEXSeqResults(dxd)
  
  # Fetch normalized counts
  counts.norm <- DEXSeq::counts(dxd, normalized=TRUE)
  counts.norm <- counts.norm[,1:nrow(sd)]
  colnames(counts.norm) <- row.names(sd)
  stopifnot(rownames(counts.norm)==peak_ref_use$key)
  rownames(counts.norm) <- peak_ref_use$peak
  
  # Start building nice results list
  dt <- data.table(as.data.frame(dxr))
  out <- with(dt,data.table(tu=groupID,peak=peak_ref_use$peak,p=pvalue,padj=padj,dexl2fc_ident1_over_ident2=get(paste0("log2fold_",ident1,"_",ident2)),
                            chr=genomicData.seqnames,start=genomicData.start,end=genomicData.end,strand=genomicData.strand,gene_name=peak_ref_use$gene,
                            final_classification=peak_ref_use$final_classification,merged_peak=peak_ref_use$merged_peak))
  
  stopifnot(out$peak==rownames(counts.norm))
  
  # Make usage matrices
  c2 <- data.table(counts)
  c2$peak <- out$peak
  c2$tu <- out$tu
  cm <- melt(c2,id.vars=c("peak","tu"))
  stopifnot(!is.na(cm$value))
  
  # Within each TU and each sample, get each percent as percent of total
  cm<-as.data.table(cm)
  m2 <- cm[,list(peak=peak,raw_count=value,use_frac=value/sum(value),tu_zero=sum(value)==0),by=c("tu","variable")]
  stopifnot(nrow(m2[is.na(use_frac),])==sum(m2$tu_zero))
  
  m2$use_frac <- ifelse(is.na(m2$use_frac),0,m2$use_frac)
  stopifnot(!is.na(m2$use_frac))
  
  m2$group <- sd[match(m2$variable,row.names(sd)),'group']
  stopifnot(!is.na(m2$group))
  
  ## Groupwise usage means
  # Fraction
  m3 <- m2[,list(tu=tu[1],mean_ident1=mean(use_frac[group==ident1]),mean_ident2=mean(use_frac[group==ident2])),by=c("peak")]
  m3$ident1_minus_ident2 <- m3$mean_ident1-m3$mean_ident2
  m3$l2_ident1_over_ident2 <- log2((m3$mean_ident1+0.0001)/(m3$mean_ident2+0.0001))
  stopifnot(out$tu==m3$tu)
  stopifnot(out$peak==m3$peak)
  setnames(m3,paste0(colnames(m3),"_frac"))
  
  cas <- cast(m2,formula="peak+tu~variable",value="use_frac")
  caskey <- paste(cas$tu,cas$peak,sep=":")
  outkey <- paste(out$tu,out$peak,sep=":")
  m <- match(outkey,caskey)
  stopifnot(!is.na(m))
  cas <- cas[m,]
  myuse.frac <- cas[,c(-1,-2)]
  row.names(myuse.frac) <- cas$peak
  colnames(myuse.frac) <- paste0(colnames(myuse.frac),"_frac")
  
  # Normalized counts
  mynorm <- counts.norm
  colnames(mynorm) <- paste0(colnames(mynorm),"_norm")
  
  nmeans_ident1 <- rowMeans(counts.norm[,row.names(sd)[which(sd$group==ident1)]])
  nmeans_ident2 <- rowMeans(counts.norm[,row.names(sd)[which(sd$group==ident2)]])
  
  dt <- data.table(nmeans_ident1=nmeans_ident1,
                   nmeans_ident2=nmeans_ident2,
                   ident1_minus_ident2=nmeans_ident1-nmeans_ident2,
                   l2_ident1_over_ident2=log2((nmeans_ident1+0.0001)/(nmeans_ident2+0.0001)))
  setnames(dt,paste0(colnames(dt),"_norm"))
  
  res <- cbind(out,mynorm,dt,myuse.frac,m3[,c("mean_ident1_frac","mean_ident2_frac","ident1_minus_ident2_frac","l2_ident1_over_ident2_frac"),with=F])
  colnames(res) <- gsub('ident1',ident1,colnames(res))
  colnames(res) <- gsub('ident2',ident2,colnames(res))
  
  # Return
  ret<-list(dxd=dxd,dxr=dxr,res=res,counts.raw=counts,counts.norm=counts.norm,use.frac=myuse.frac,ident1=ident1,ident2=ident2,sampleData=sd)
  
  return(ret)
}

APA_TTest <- function(myta,ident1,ident2){
  cat("T Testing of fraction used:\n")
  dat <- myta$use.frac
  ident1col_idx <- grep(ident1,colnames(DEXseq_res$use.frac))
  ident2col_idx <- grep(ident2,colnames(DEXseq_res$use.frac))
  # Found cases where both groups have zero variance and t-test failed on them, so have to flag this and pass 0 p-value
  t.p <- lapply(1:nrow(dat),function(x){
    
    avec <- dat[x,][,ident1col_idx]
    bvec <- dat[x,][,ident2col_idx]
    
    # Flag with 0 p-value cases where fractions used within each group equal each other and so have no variance and led to t test failing 
    if((sum(avec==max(avec))==length(ident1col_idx))&(sum(bvec==max(bvec))==length(ident2col_idx)))
    {return(0)} else {
      td <- t.test(x=avec,y=bvec)
      return(td$p.value)}
  })
  t.p <- do.call(c,t.p)
  t.padj <- p.adjust(t.p,method="fdr")
  myta$res$ttest_p <- t.p
  myta$res$ttest_padj <- t.padj
  return(myta)
}

APA_EdgeR_test <- function(myta,ident1,ident2,edger_prefix=paste0(ident1,"_v_",ident2),outdir='./'){
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  
  cat("EdgeR Testing:",edger_prefix,"\n")
  
  res <- myta$res
  ann <- data.frame(tu=res$tu,pr=res$peak)
  group <- as.character(myta$sampleData$group)
  #adj <- factor(myta$stab.sub$batch)
  cnt <- myta$counts.raw
  y.all <- DGEList(counts=cnt,genes=ann,group=group)
  # Lib size normalization
  y <- calcNormFactors(y.all)
  # Create design matrix
  design <- model.matrix(~0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  makeContrasts_cmd <- paste0('makeContrasts(',ident1,'-',ident2,',levels=design)')
  contr <- eval(parse(text=makeContrasts_cmd))
  # Estimate NB dispersion
  y <- estimateDisp(y, design, robust=TRUE)
  # Estimate QL dispersion
  fit <- glmQLFit(y, design, robust=TRUE)
  
  # Generate plot  
  pdf(file=paste0(outdir,edger_prefix,"_BCV_QLDisp_MDS.pdf"))
  plotBCV(y)
  plotQLDisp(fit)
  plotMDS(y)
  dev.off()
  
  # Make contrast & differential exon usage test
  sp <- diffSpliceDGE(fit, contrast = contr, geneid="tu", exonid="pr")
  
  myta$edger <- list(NB_disp = y, QL_disp = fit, diff_exon_res = sp)
  
  # Not all PR in sp
  # Differential exon test result
  idx<-match(names(sp$exon.p.value),myta$res$peak)
  myta$res$edger_exon_p<-rep(NA,nrow(myta$res))
  myta$res$edger_exon_p[idx] <- sp$exon.p.value
  myta$res$edger_exon_padj[idx] <- p.adjust(sp$exon.p.value,method="fdr")
  
  # Gene level test result from F test
  ftest <- sp$gene.p.value
  ftest.adj <- ftest
  ftest.adj[,1] <- p.adjust(ftest[,1],method="fdr")
  
  # Gene level test result from Simes test
  simes <- sp$gene.Simes.p.value
  names(simes) <- res[match(names(simes),res$peak),]$tu
  simes.adj <- p.adjust(simes,method="fdr")
  
  myta$res$edger_ftest_p <- ftest[match(res$tu,rownames(ftest))]
  myta$res$edger_ftest_padj <- ftest.adj[match(res$tu,rownames(ftest.adj))]
  
  myta$res$edger_simes_p <- simes[match(res$tu,names(simes))]
  myta$res$edger_simes_padj <- simes.adj[match(res$tu,names(simes.adj))]
  
  return(myta)
}

callsig <- function(myta,ident1,ident2,delta_frac_thres=0.1, padj_thres=0.0001,
                    l2fc_frac_thres=log2(1.5),mean_frac_thres=0.05)
{
  ident1_minus_ident2_frac <- paste0(ident1,'_minus_',ident2,'_frac')
  l2_ident1_over_ident2_frac <- paste0('l2_',ident1,'_over_',ident2,'_frac')
  mean_ident1_frac <- paste0('mean_',ident1,'_frac')
  mean_ident2_frac <- paste0('mean_',ident2,'_frac')
  
  # P-value and fraction delta
  cat('Check what peaks have absolute delta fraction used at least',delta_frac_thres,'AND DEXseq adjusted p value less than',padj_thres,'\n')
  myta$res$delta_sig <- ((abs(myta$res[[ident1_minus_ident2_frac]])>=delta_frac_thres)&(myta$res$padj<padj_thres))
  
  # P-value and fold change
  cat('Check what peaks have DEXseq adjusted p value less than',padj_thres,'AND absolute log2 fold change of fraction used at least',l2fc_frac_thres,
      '\nAND the mean fraction used of one or more groups at least',mean_frac_thres,'\n')
  myta$res$perfc_sig <- ((myta$res$padj<padj_thres)&(abs(myta$res[[l2_ident1_over_ident2_frac]])>=l2fc_frac_thres)&
                           ((myta$res[[mean_ident1_frac]]>=mean_frac_thres)|(myta$res[[mean_ident2_frac]]>=mean_frac_thres)))
  
  # P-value and fold change and delta fraction (intersect set we want to use)
  cat('Intersect of those two sets \n')
  myta$res$int_sig <- myta$res$delta_sig & myta$res$perfc_sig 
  
  cat('Peaks significate under DEXseq, T test, or edgeR (exon-level test) \n')
  myta$res$dexseq_sig <- myta$res$padj < padj_thres
  myta$res$ttest_sig <- myta$res$ttest_padj < padj_thres
  #myta$res$edger_sig <- myta$res$edger_exon_padj < padj_thres
  
  cat('Significant under both dexseq and t tests \n')
  #myta$res$all_tests_sig <- (myta$res$dexseq_sig & myta$res$ttest_sig & myta$res$edger_sig)
  myta$res$all_tests_sig <- (myta$res$dexseq_sig & myta$res$ttest_sig)
  
  return(myta)
}

### DEG
DEG_DESeq2_pseudobulk <- function(inputs,comp,test.used='LRT',
                                  padj_thres=0.0001,l2fc_frac_thres=log2(1.5),outdir) {
  sd <- data.frame(sample_id = colnames(inputs), group_id = gsub('_rep\\d+$','',colnames(inputs)))
  rownames(sd) <- sd$sample_id
  sd$group_id <- factor(sd$group_id, levels=comp)
  
  stopifnot(all(rownames(sd) == colnames(inputs)))  
  
  dds <- DESeqDataSetFromMatrix(inputs, 
                                colData = sd, 
                                design = ~ group_id)
  
  rld <- rlog(dds, blind=TRUE)
  DESeq2::plotPCA(rld, intgroup = "group_id")
  ggsave(paste0(outdir, levels(sd$group_id)[2], "_vs_", levels(sd$group_id)[1], "_specific_PCAplot.png"), width = 12, height = 8, units = 'in')
  
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  png(paste0(outdir, levels(sd$group_id)[2], "_vs_", levels(sd$group_id)[1], "_specific_heatmap.png"), width = 12, height = 8, units = 'in', res = 1200)
  pheatmap(rld_cor, annotation = sd[, c("group_id"), drop=F])
  dev.off()
  
  #Test
  
  median_abs_residual <- c()
  reduced.formula = as.formula(~ 1)
  
  fittypes = c('parametric','local','mean')
  if (test.used == 'LRT') {fittypes <- c(fittypes,'glmGamPoi')}
  for  (fittype in fittypes) {
    if (test.used == 'LRT') {
      dds <- DESeq(dds, test = test.used, fitType = fittype, sfType = 'ratio',reduced = reduced.formula)
    } else if (test.used == 'Wald') {dds <- DESeq(dds, test = test.used, fitType = fittype, sfType = 'ratio')}
    mar <- abs(mcols(dds)$dispGeneEst - mcols(dds)$dispFit) %>% median(na.rm=TRUE)
    median_abs_residual <- c(median_abs_residual,mar)
    png(paste0(outdir, levels(sd$group_id)[2], "_vs_", levels(sd$group_id)[1], "_dispersion_plot_fittype_",fittype,".png"), width = 6, height = 6, units = 'in', res = 1200)
    plotDispEsts(dds)
    dev.off()
  }
  names(median_abs_residual) <- fittypes
  chosen_fittype <- which(median_abs_residual == min(median_abs_residual)) %>% names()
  #chosen_fittype <- 'glmGamPoi'
  cat('The fit type with the least median absolute residual is',chosen_fittype,'\n')
  cat('Differential testing using',chosen_fittype,'as fit type and',test.used,'as test \n')
  
  if (test.used == 'LRT') {
    dds <- DESeq(dds, test = test.used, fitType = chosen_fittype, sfType = 'ratio',reduced = reduced.formula)
  } else if (test.used == 'Wald') {dds <- DESeq(dds, test = test.used, fitType = chosen_fittype, sfType = 'ratio')}
  
  levels(sd$group_id)[2]
  levels(sd$group_id)[1]
  
  resultsNames(dds)
  if (test.used == 'LRT') {
    res <- results(dds, alpha = 0.05)
  } else if (test.used == 'Wald') {
    contrast <- c("group_id", levels(sd$group_id)[2], levels(sd$group_id)[1])
    res <- results(dds, contrast = contrast, alpha = 0.05)
  }
  
  # betaPrior is set to zero so lfcShrink now 
  if (test.used == 'Wald') {
    coef <- paste("group_id",levels(sd$group_id)[2], "vs", levels(sd$group_id)[1], sep = '_')
    cat('betaZero is set to zero so lfcShrink with coef ',coef,'\n')
    res <- lfcShrink(dds, coef = coef, res=res)}
  
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  cat('Filter for genes with adjusted p values less than',padj_thres,' and absolute log2 fold change at least',l2fc_frac_thres,'\n')
  sig_genes <- res_tbl %>% filter(padj < padj_thres & abs(log2FoldChange) >= l2fc_frac_thres)
  
  # Save sig results
  write.table(sig_genes,paste0(outdir,levels(sd$group_id)[2], "_vs_", levels(sd$group_id)[1], "_",test.used,"_sig_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
  # Save all genes
  write.table(res_tbl,paste0(outdir,levels(sd$group_id)[2], "_vs_", levels(sd$group_id)[1], "_",test.used,"_all_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
  
  return(sig_genes)
}