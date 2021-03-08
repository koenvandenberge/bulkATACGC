### GCQN, first implementation
FQnorm <- function(counts, type="mean"){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  if(type=="mean"){
    refdist <- apply(counts.sort,1,mean)
  } else if(type=="median"){
    refdist <- apply(counts.sort,1,median)
  }
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}


gcqn_mean_scone <- function(ei){
  gcBinNormCounts <- matrix(NA, nrow=nrow(ei), ncol=ncol(ei), dimnames=list(rownames(ei),colnames(ei)))
  for(ii in 1:nlevels(gcGroupsScone)){
    id <- which(gcGroupsScone==levels(gcGroupsScone)[ii])
    countBin <- ei[id,]
    normCountBin <- FQnorm(countBin, type='mean')
    # normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}

gcqn_median_20_scone <- function(ei){
  gcBinNormCounts <- matrix(NA, nrow=nrow(ei), ncol=ncol(ei), dimnames=list(rownames(ei),colnames(ei)))
  for(ii in 1:nlevels(gcGroupsScone)){
    id <- which(gcGroupsScone==levels(gcGroupsScone)[ii])
    countBin <- ei[id,]
    normCountBin <- FQnorm(countBin, type='median')
    # normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}

gcqn_median_10_scone <- function(ei){
  gcBinNormCounts <- matrix(NA, nrow=nrow(ei), ncol=ncol(ei), dimnames=list(rownames(ei),colnames(ei)))
  for(ii in 1:nlevels(gcGroupsScone10)){
    id <- which(gcGroupsScone10==levels(gcGroupsScone10)[ii])
    countBin <- ei[id,]
    normCountBin <- FQnorm(countBin, type='median')
    # normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}

gcqn_median_50_scone <- function(ei){
  gcBinNormCounts <- matrix(NA, nrow=nrow(ei), ncol=ncol(ei), dimnames=list(rownames(ei),colnames(ei)))
  for(ii in 1:nlevels(gcGroupsScone50)){
    id <- which(gcGroupsScone50==levels(gcGroupsScone50)[ii])
    countBin <- ei[id,]
    normCountBin <- FQnorm(countBin, type='median')
    # normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}

gcqn_median_100_scone <- function(ei){
  gcBinNormCounts <- matrix(NA, nrow=nrow(ei), ncol=ncol(ei), dimnames=list(rownames(ei),colnames(ei)))
  for(ii in 1:nlevels(gcGroupsScone100)){
    id <- which(gcGroupsScone100==levels(gcGroupsScone100)[ii])
    countBin <- ei[id,]
    normCountBin <- FQnorm(countBin, type='median')
    # normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}


gcqn_qsmooth_scone <- function(ei){
  gcBinNormCounts <- matrix(NA, nrow=nrow(ei), ncol=ncol(ei), dimnames=list(rownames(ei),colnames(ei)))
  for(ii in 1:nlevels(gcGroupsScone)){
    id <- which(gcGroupsScone==levels(gcGroupsScone)[ii])
    countBin <- ei[id,]
    qs <- qsmooth(countBin, group_factor=bio)
    normCountBin <- qs@qsmoothData
    # normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}

qsmooth_scone <- function(ei){
    qs <- qsmooth(ei, group_factor=bio)
    normCounts <- qs@qsmoothData
    normCounts[normCounts<0] <- 0
  return(normCounts)
}

cqn_scone <- function(ei){
  library(cqn)
  cqnModel <- cqn(ei, x=gcContentScone, lengths=1000, lengthMethod="fixed", 
                  sizeFactors = colSums(ei))
  countsCqn <- 2^(cqnModel$y + cqnModel$offset)
  return(countsCqn)
}

cqn_scone_length <- function(ei){
  library(cqn)
  cqnModel <- cqn(ei, x=gcContentScone, lengths=width(gr), 
                  sizeFactors = colSums(ei))
  countsCqn <- 2^(cqnModel$y + cqnModel$offset)
  return(countsCqn)
}


edaseq <- function(ei){
  library(EDASeq)
  wit1 <- withinLaneNormalization(ei, gcContentScone, which="full")
  bet1 <- betweenLaneNormalization(wit1, which="full")
  return(bet1)
}





sconeEvaluation <- function(counts, bio, qc=NULL, batch=NULL, negcon, gcContent, gcGroups, 
                            k_ruv=3, k_qc=0, adjust_bio="no", eval_kclust=4,
                            diffLength=FALSE, ...){

  exprScone <- as.matrix(counts)
  rownames(exprScone) <- rownames(counts)
  
  # scone can't handle peaks having identical counts so remove these

  if(is.null(qc)){
    if(is.null(batch)){
      my_scone <- SconeExperiment(exprScone,
                                  bio = bio,
                                  negcon_ruv = rownames(exprScone) %in% negcon)
    } else {
      my_scone <- SconeExperiment(exprScone,
                                  bio = bio,
                                  negcon_ruv = rownames(exprScone) %in% negcon,
                                  batch = batch)
    }
    
  } else if(!is.null(qc)){
    if(is.null(batch)){
      my_scone <- SconeExperiment(exprScone,
                                  bio = bio,
                                  negcon_ruv = rownames(exprScone) %in% negcon,
                                  qc = qc)
    } else {
      my_scone <- SconeExperiment(exprScone,
                                  bio = bio,
                                  negcon_ruv = rownames(exprScone) %in% negcon,
                                  qc = qc,
                                  batch = batch)
    }
  }
    
 
  
  # normalization methods to try
  if(!diffLength){
    scaling=list(none=identity, # Identity - do nothing
                 gcqn_mean = gcqn_mean_scone, # User-defined function
                 gcqn_median_20 = gcqn_median_20_scone,
                 cqn = cqn_scone,
                 gcqn_smooth = gcqn_qsmooth_scone,
                 qsmooth = qsmooth_scone,
                 sum = SUM_FN, # SCONE library wrappers...
                 tmm = TMM_FN,
                 uq = UQ_FN,
                 fq = FQT_FN,
                 deseq = DESEQ_FN,
                 edaseq = edaseq,
                 gcqn_median_10 = gcqn_median_10_scone,
                 gcqn_median_50 = gcqn_median_50_scone,
                 gcqn_median_100 = gcqn_median_100_scone)
  } else if(diffLength){
    scaling=list(none=identity, # Identity - do nothing
                 gcqn_mean = gcqn_mean_scone, # User-defined function
                 gcqn_median_20 = gcqn_median_20_scone,
                 cqn = cqn_scone,
                 cqn_length = cqn_scone_length,
                 gcqn_smooth = gcqn_qsmooth_scone,
                 qsmooth = qsmooth_scone,
                 sum = SUM_FN, # SCONE library wrappers...
                 tmm = TMM_FN,
                 uq = UQ_FN,
                 fq = FQT_FN,
                 deseq = DESEQ_FN,
                 edaseq = edaseq,
                 gcqn_median_10 = gcqn_median_10_scone,
                 gcqn_median_50 = gcqn_median_50_scone,
                 gcqn_median_100 = gcqn_median_100_scone)
  }
  
  
  # run
  my_scone <- scone(my_scone,
                    scaling=scaling,
                    k_ruv = k_ruv, 
                    k_qc=k_qc,
                    adjust_bio=adjust_bio,
                    eval_kclust=eval_kclust,
                    stratified_pam = TRUE,
                    ...)
  return(my_scone)
}

biplot_colorIk <- function (x, y, rank = TRUE, ties_method = c("max", "min", "first", 
                                             "last", "random"), choices = 1:2, expand = 1, ...) 
{
  if (rank) {
    ties_method <- match.arg(ties_method)
    y = rank(y, ties.method = ties_method)
  }
  else {
    if (any(abs(y - round(y)) > .Machine$double.eps^0.5)) {
      stop("ranks must be integer")
    }
    else {
      y = as.integer(y)
    }
    if (any(y <= 0)) {
      stop("ranks must be positive")
    }
    if (any(y > length(y))) {
      stop("ranks must be less than or equal to total number of elements")
    }
  }
  lam <- x$sdev[choices]
  n <- NROW(x$x)
  lam <- lam * sqrt(n)
  xx <- t(t(x$x[, choices])/lam)
  yy <- t(t(x$rotation[, choices]) * lam)
  ratio <- max(range(yy)/range(xx))/expand
  cols <- rev(colorRampPalette(c("black", "navyblue", "mediumblue", 
                                 "dodgerblue3", "aquamarine4", "green4", "yellowgreen", 
                                 "yellow"))(length(y)))[y]
  # plot(xx, pch = 19, col = cols, ...)
  plot(xx, ...)
  labs <- rownames(yy)
  text(yy/ratio, labels = labs, col = 2)
  arrows(0, 0, yy[, 1] * 0.8/ratio, yy[, 2] * 0.8/ratio, col = 2, 
         length = 0.1)
  invisible(xx)
}


# to quantify differential GC content effects we need to quantify it across 
# samples within a bin and then assess how much this varies across bins
# note that, due to the correlation between count and GC content, the variance
# of the IQR is expected to change across GC bins but this does not necessarily
# reflect GC bias. Therefore, evaluating the medians seems to be more appropriate.

rleGC_med <- function(counts, gcContent, binSize=4e3){
  logCounts <- log1p(counts)
  meds <- rowMedians(as.matrix(logCounts))
  rle <- sweep(as.matrix(logCounts),1,meds,FUN="-") # log scale
  gcBins <- Hmisc::cut2(gcContent, g=round(nrow(logCounts)/binSize))
  rle_med <- vector(length=ncol(logCounts))
  for(kk in 1:nlevels(gcBins)){
    rleBin <- rle[gcBins == levels(gcBins)[kk],]
    binMeds <- matrixStats::colMedians(rleBin)
    rle_med[kk] <- mean(binMeds^2)
  }
  rleMedsGC <- var(rle_med)
}

rleGC_iqr <- function(counts, gcContent, binSize=4e3){
  logCounts <- log1p(counts)
  meds <- rowMedians(as.matrix(logCounts))
  rle <- sweep(as.matrix(logCounts),1,meds,FUN="-") # log scale
  gcBins <- Hmisc::cut2(gcContent, g=round(nrow(logCounts)/binSize))
  rle_iqr <- vector(length=ncol(logCounts))
  for(kk in 1:nlevels(gcBins)){
    rleBin <- rle[gcBins == levels(gcBins)[kk],]
    binIQRs <- matrixStats::colIQRs(rleBin)
    rle_iqr[kk] <- var(binIQRs)
  }
  rleMedsGC <- var(rle_iqr)
}

rleGCh5 <- function(h5File, gcContent, type="med", ...){
  h5List <- h5ls(h5File)
  rleGCVals <- sapply(2:(nrow(h5List)-1), function(normId){
    normCounts <- h5read(h5File, h5List$name[normId])
    if(type == "med"){
      rleValue <- rleGC_med(counts=normCounts, gcContent=gcContent, ...)
    } else if(type == "iqr"){
      rleValue <- rleGC_iqr(counts=normCounts, gcContent=gcContent, ...)
    }
    return(rleValue)
  })
  names(rleGCVals) <- h5List$name[2:(nrow(h5List)-1)]
  return(rleGCVals)
}









