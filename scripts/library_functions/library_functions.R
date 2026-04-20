# QC function
qc_exons <- function(Y, ref, 
                          cov_threshold = c(10, 500), 
                          length_threshold = c(50, 2000), 
                          mappability_threshold = 0.9, 
                          gc_threshold = c(0.2, 0.8)) {
  
  # Extract exon-level covariates
  exon_lengths <- width(ref)
  mapp <- ref$mapp
  gc <- ref$gc
  
  # Median coverage across all samples (no sample filtering)
  median_cov <- apply(Y, 1, median, na.rm = TRUE)  
  # Exon filters
  keep_cov <- median_cov >= cov_threshold[1] & median_cov <= cov_threshold[2]
  keep_len <- exon_lengths >= length_threshold[1] & exon_lengths <= length_threshold[2]
  keep_mapp <- mapp >= mappability_threshold
  keep_gc <- gc >= gc_threshold[1] & gc <= gc_threshold[2]
  
  # Final filter
  keep_all <- keep_cov & keep_len & keep_mapp & keep_gc
  
  # Log filter stats
  message(sprintf("Coverage filter: excluded %d exons", sum(!keep_cov)))
  message(sprintf("Length filter: excluded %d exons", sum(!keep_len)))
  message(sprintf("Mappability filter: excluded %d exons", sum(!keep_mapp)))
  message(sprintf("GC content filter: excluded %d exons", sum(!keep_gc)))
  message(sprintf("Total exons retained: %d / %d", sum(keep_all), length(keep_all)))
  
  # QC matrix
  qc_matrix <- data.frame(
    chr = as.character(seqnames(ref)),
    start = start(ref),
    end = end(ref),
    pass_all = keep_all,
    median_cov = median_cov,
    pass_cov = keep_cov,
    length_kb = exon_lengths / 1000,
    pass_len = keep_len,
    mapp = mapp,
    pass_mapp = keep_mapp,
    gc = round(gc, 2),
    pass_gc = keep_gc
  )
  
  # Return filtered objects
  list(
    Y_qc = Y[keep_all, , drop = FALSE],
    ref_qc = ref[keep_all],
    qc_matrix = qc_matrix
  )
}




# Normalization main function

normalization_bin_ref <- function(Y, bin_size, gc, map, ref_id){
  
  RCTL <- Y
  RC_norm=matrix(data=NA,nrow = dim(RCTL)[1],ncol = dim(RCTL)[2])
  rownames(RC_norm)=rownames(RCTL)
  colnames(RC_norm)=colnames(RCTL)
  for(sub in 1:dim(RCTL)[2]){
    ### Bin size normalization only if InTarget###
    step <- 50
    RCLNormListIn <- CorrectSize(RCTL[,sub],L=bin_size,step)
    RCLNormIn <- RCLNormListIn$RCNorm
    
    ### Mappability normalization ###
    step <- 0.02
    RCMAPNormListIn <- CorrectMAP(RCLNormIn,MAPContent=map,step)
    RCMAPNormIn <- RCMAPNormListIn$RCNorm
    
    ### GC-content Normalization ###
    step <- 0.02
    RCGCNormListIn <- CorrectGC(RCMAPNormIn,GCContent=gc,step)
    RCGCNormIn <- RCGCNormListIn$RCNorm
    
    RC_norm[,sub]=RCGCNormIn
  }
  # dont add the rgamma error term this time
  # RC_norm <- RC_norm + matrix(rgamma(nrow(RC_norm)*ncol(RC_norm),0.1,10),nrow(RC_norm),ncol(RC_norm))
  
  # Calculate L2R based on RC_norm and ref_id
  RC_norm1 <- RC_norm[, !(colnames(RC_norm) %in% ref_id)]
  ref_norm <- RC_norm[, ref_id]
  
  lrr=matrix(data=NA,nrow = dim(RC_norm1)[1],ncol = dim(RC_norm1)[2])
  rownames(lrr)=rownames(RC_norm1)
  colnames(lrr)=colnames(RC_norm1)
  for(i in 1:dim(RC_norm1)[1]){
    mean=mean(ref_norm[i,])
    rc=RC_norm1[i,]
    # lrr[i,]=log2((rc+1)/(mean+1))
    # change from plus 1 to a smaller value
    lrr[i,]=log2((rc+0.01)/(mean+0.01))
    
  }
  log2Rdata=lrr
  
  return(list(log2Rdata = log2Rdata, RC_norm =RC_norm))
}

CorrectGC<-function(RC,GCContent,step){
  stepseq<-seq(min(GCContent),max(GCContent),by=step)
  #stepseq<-seq(0,100,by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianGC<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(GCContent>=stepseq[i] & GCContent<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(GCContent>stepseq[i] & GCContent<=stepseq[i+1])
    }
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianGC[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianGC
  RCNormList$StepGC<-stepseq[1:(length(stepseq)-1)]
  RCNormList$RCNorm<-RCNormMedian
  RCNormList
}


CorrectSize<-function(RC,L,step){
  stepseq<-seq(min(L),max(L),by=step)
  #stepseq<-seq(0,max(L),by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianL<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(L>=stepseq[i] & L<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(L>stepseq[i] & L<=stepseq[i+1])
    }
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianL[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianL
  RCNormList$RCNorm<-RCNormMedian
  RCNormList
}


CorrectMAP<-function(RC,MAPContent,step){
  stepseq<-seq(min(MAPContent),max(MAPContent),by=step)
  #stepseq<-seq(0,100,by=step)
  MasterMedian<-median(RC,na.rm=T)
  MedianMAP<-rep(0,length(stepseq)-1)
  RCNormMedian<-RC
  for (i in 1:(length(stepseq)-1)){
    if (i==1){
      ind<-which(MAPContent>=stepseq[i] & MAPContent<=stepseq[i+1])
    }
    if (i!=1){
      ind<-which(MAPContent>stepseq[i] & MAPContent<=stepseq[i+1])
    }
    
    if (length(ind)>0){
      m<-median(RC[ind],na.rm=T)
      if (m>0){
        MedianMAP[i]<-m
        RCNormMedian[ind]<-RC[ind]*MasterMedian/m
      }
    }
  }
  RCNormList<-list()
  RCNormList$Median<-MedianMAP
  RCNormList$StepMAP<-stepseq[1:(length(stepseq)-1)]
  RCNormList$RCNorm<-RCNormMedian
  RCNormList
}