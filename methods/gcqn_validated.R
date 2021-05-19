logit <- function(x) log(x/(1-x))


### GCQN, first implementation
FQnorm <- function(counts, type="mean"){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  if(type=="mean"){
    # refdist <- apply(counts.sort,1,mean)
    refdist <- base::rowMeans(counts.sort)
  } else if(type=="median"){
    #refdist <- apply(counts.sort,1,median)
    refdist <- matrixStats::rowMedians(counts.sort)
  }
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}


gcqn <- function(counts, gcGroups, summary='mean', round=TRUE){
  gcBinNormCounts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts), dimnames=list(rownames(counts),colnames(counts)))
  for(ii in 1:nlevels(gcGroups)){
    id <- which(gcGroups==levels(gcGroups)[ii])
    if(length(id) == 1){
      normCountBin <- counts[id,]
      if(round) normCountBin <- round(normCountBin)
      gcBinNormCounts[id,] <- normCountBin
      next
    }
    countBin <- counts[id,,drop=FALSE]
    if(summary=="mean"){
      normCountBin <- FQnorm(countBin, type='mean')
    } else if(summary=="median"){
      normCountBin <- FQnorm(countBin, type='median')
    }
    if(round) normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}


gcqn_qsmooth <- function(counts, gcGroups, bio){
  gcBinNormCounts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts), dimnames=list(rownames(counts),colnames(counts)))
  for(ii in 1:nlevels(gcGroups)){
    id <- which(gcGroups==levels(gcGroups)[ii])
    countBin <- counts[id,]
    qs <- qsmooth(countBin, group_factor=bio)
    normCountBin <- qs@qsmoothData
    normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}



