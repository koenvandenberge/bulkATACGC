## startified split of samples according to bio variable
splitSamples <- function(bio, seed=15){
  grp <- rep(NA, length(bio)) # group assignments
  grpVals <- c("A", "B")
  for(bb in 1:nlevels(bio)){
    id <- which(bio == levels(bio)[bb])
    if(length(id) > 1){
      halfLength <- length(id)/2
      g1 <- sample(grpVals, 1)
      if(halfLength%%1 == 0){
        grp[id[1:halfLength]] <- g1
        grp[id[(halfLength+1):length(id)]] <- grpVals[!grpVals %in% g1]
      } else {
        grp[id[1:floor(halfLength)]] <- g1
        grp[id[ceiling(halfLength):length(id)]] <- grpVals[!grpVals %in% g1]
      }
    }
  }
  return(grp)
}

## simulate data
simulateData <- function(counts, 
                         grp, 
                         nDE = round(nrow(counts)/10),
                         nTotal = ncol(counts),
                         meanlog=0.8, 
                         sdlog=0.1, 
                         plot=FALSE,
                         seed=15){
  
  if(plot) par(mfrow=c(1,2))
  set.seed(seed)
  n1 <- 1:(nTotal/2)
  n2 <- (nTotal/2+1):nTotal
  fractions <- sweep(counts, 2, colSums(counts), "/")
  simCounts <- matrix(NA, nrow=nrow(counts), ncol=nTotal,
                           dimnames=list(rownames(counts), paste0("sample", 1:nTotal)))
  remainingSamplesA <- which(grp == "A")
  remainingSamplesB <- which(grp == "B")
  
  deId <- sample(x=1:nrow(counts), size=nDE)
  delta <- rep(0, nrow(counts))
  sign <- rep(c(1,-1), each=nDE/2)
  if(length(sign) != nDE) sign <- c(sign,1)
  delta[deId] <- sign*log(rlnorm(n=nDE, meanlog=meanlog, sdlog=sdlog))
  if(plot) hist(delta[deId],  breaks=40)
  
  # group 1
  for(ss in n1){
    sampleId <- sample(x=remainingSamplesA, size=1)
    remainingSamplesA <- remainingSamplesA[!remainingSamplesA %in% sampleId]
    curLS <- runif(n=1, min(colSums(counts)), max(colSums(counts)))
    curFracs <- fractions[,sampleId]
    curFracs[delta < 0] <- (curFracs[delta < 0]) * exp(abs(delta[delta < 0]))
    simY <- rmultinom(n=1, size=curLS, prob=curFracs)
    simCounts[,ss] <- simY
  }
  
  # group 2
  for(ss in n2){
    sampleId <- sample(x=remainingSamplesB, size=1)
    remainingSamplesB <- remainingSamplesB[!remainingSamplesB %in% sampleId]
    curLS <- runif(n=1, min(colSums(counts)), max(colSums(counts)))
    curFracs <- fractions[,sampleId]
    curFracs[delta > 0] <- (curFracs[delta > 0]) * exp(abs(delta[delta > 0]))
    simY <- rmultinom(n=1, size=curLS, prob=curFracs)
    simCounts[,ss] <- simY
  }
  
  if(plot){
    plot(x=delta, y=log(rowMeans(simCounts[,n2]) / rowMeans(simCounts[,n1])),
         pch=16, cex=1/2)
    abline(0,1, col="red")
    abline(h=0, col="blue", lty=2)
  }
  
  grpSim <- factor(rep(c("A", "B"), each=nTotal/2))
  
  
  return(list(simCounts = simCounts,
              grp = grpSim,
              deId = deId))
}



evaluateSimulation <- function(simCounts, 
                               simGC, 
                               simWidth, 
                               design, 
                               deId, 
                               grSim, 
                               grpSim,
                               ruv=FALSE,
                               cqnWithLength=FALSE){
  
  grpSim <- factor(grpSim)
  design <- model.matrix(~grpSim)
  
  testEdgeR <- function(counts, design, tmm=FALSE, offset=NULL, uq=FALSE){
    d <- DGEList(counts)
    if(uq) d <- calcNormFactors(d, method = "upperquartile")
    if(tmm) d <- calcNormFactors(d)
    if(!is.null(offset)) d$offset <- offset
    d <- estimateDisp(d, design)
    fit <- glmFit(d, design)
    lrt <- glmLRT(fit, coef=2)
    lrt$table$padj <- p.adjust(lrt$table$PValue, "fdr")
    return(lrt)
  }
  
  library(DESeq2)
  testDESeq2 <- function(counts, covarDf){
    dds <- DESeqDataSetFromMatrix(counts,
                                  colData=covarDf,
                                  design=as.formula("~ group"))
    dds <- DESeq(dds)
    res <- results(dds)
    return(res)
  }
  
  ## normalize
  ## TMM
  d <- DGEList(simCounts)
  d <- calcNormFactors(d)
  cpmTMM <- edgeR::cpm(d, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)
  
  ## FQN
  fqCounts <- FQnorm(simCounts)
  cpmFQ <- edgeR::cpm(fqCounts, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
  
  ## qsmooth
  qsmoothObj <- qsmooth(object = simCounts, group_factor = grpSim)
  qsmoothCounts <- qsmoothData(qsmoothObj)
  
  ## cqn
  if(cqnWithLength){
    cqnObj <- cqn(simCounts, x=simGC, lengths=simWidth, sizeFactors=colSums(simCounts))
    cqnplot(cqnObj, xlab="GC content", lwd=2)
    cqnCounts <- 2^(cqnObj$y + cqnObj$offset)
    cpmCqn <- edgeR::cpm(cqnCounts, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
  } else {
    cqnObj <- cqn(simCounts, x=simGC, sizeFactors=colSums(simCounts),
                  lengthMethod = "fixed", lengths = 1000)
    cqnplot(cqnObj, xlab="GC content", lwd=2)
    cqnCounts <- 2^(cqnObj$y + cqnObj$offset)
    cpmCqn <- edgeR::cpm(cqnCounts, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
  }
  
  
  ## GCQN:
  gcGroups <- Hmisc::cut2(simGC, g=50)
  countsGCQN <- gcqn(simCounts, gcGroups, summary="median", round=TRUE)
  cpmGCQN <- edgeR::cpm(countsGCQN, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
  
  ## smooth GCQN
  normCountsGcqnSmooth <- gcqn_qsmooth(simCounts, gcGroups, bio=grpSim)
  
  ## EDASeq
  library(EDASeq)
  wit <- withinLaneNormalization(x=as.matrix(simCounts), y=simGC, which="full")
  countsEDASeq <- betweenLaneNormalization(wit, "full")
  cpmEDASeq <- edgeR::cpm(countsEDASeq, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
  
  ## RUVg
  if(ruv){
    library(RUVSeq)
    hk <- readRDS("~/data/genomes/hkListHumanGenomicRanges.rds")
    qh <- queryHits(findOverlaps(grSim[-deId], hk, type="within"))
    #qh <- queryHits(findOverlaps(grSim, hk, type="within"))
    negcon <- rownames(simCounts)[qh]
    ruvRes <- RUVg(as.matrix(simCounts), negcon, k=2)
    ruvgCounts <- ruvRes$normalizedCounts
    cpmRUV <- edgeR::cpm(ruvgCounts, normalized.lib.sizes=FALSE, log=TRUE, prior.count=0.25)
  }
  
  
  # none
  resNone <- testEdgeR(simCounts, design, tmm=FALSE)
  # upper-quartile
  resUQ <- testEdgeR(simCounts, design, tmm=FALSE, uq=TRUE)
  # TMM
  resTMM <- testEdgeR(simCounts, design, tmm=TRUE)
  # DESeq2 MOR
  resDESeq2 <- testDESeq2(simCounts, covarDf=data.frame(group=grpSim))
  # Cqn
  resCqn <- testEdgeR(simCounts, design, tmm=FALSE, offset=cqnObj$glm.offset)
  #FQ
  resFQ <- testEdgeR(fqCounts, design, tmm=FALSE)
  # qsmooth
  resQsmooth <- testEdgeR(qsmoothCounts, design, tmm=FALSE)
  # EDASEQ
  resEDASeq <- testEdgeR(countsEDASeq, design, tmm=FALSE)
  # GCQN
  resGCQN <- testEdgeR(countsGCQN, design, tmm=FALSE)
  # smooth GCQN
  resGCQNSmooth <- testEdgeR(normCountsGcqnSmooth, design, tmm=FALSE)
  # RUVg
  if(ruv) resRUVg <- testEdgeR(ruvgCounts, design, tmm=FALSE)
  
  truth <- rep(0, nrow(simCounts))
  truth[deId] <- 1
  truthDf <- data.frame(truth=truth)
  rownames(truthDf) <- rownames(simCounts)
  
  library(iCOBRA)
  if(ruv){
    cbd <- COBRAData(pval = data.frame(none=resNone$table$PValue,
                                       uq=resUQ$table$PValue,
                                       TMM=resTMM$table$PValue,
                                       DESeq2=resDESeq2$pvalue,
                                       qsmooth=resQsmooth$table$PValue,
                                       cqn=resCqn$table$PValue,
                                       FQ=resFQ$table$PValue,
                                       edaseq=resEDASeq$table$PValue,
                                       gcqn=resGCQN$table$PValue,
                                       gcqn_smooth=resGCQNSmooth$table$PValue,
                                       RUVg=resRUVg$table$PValue,
                                       row.names=rownames(simCounts)),
                     truth = truthDf)
  } else {
    cbd <- COBRAData(pval = data.frame(none=resNone$table$PValue,
                                       uq=resUQ$table$PValue,
                                       TMM=resTMM$table$PValue,
                                       DESeq2=resDESeq2$pvalue,
                                       qsmooth=resQsmooth$table$PValue,
                                       cqn=resCqn$table$PValue,
                                       FQ=resFQ$table$PValue,
                                       edaseq=resEDASeq$table$PValue,
                                       gcqn=resGCQN$table$PValue,
                                       gcqn_smooth=resGCQNSmooth$table$PValue,
                                       row.names=rownames(simCounts)),
                     truth = truthDf)
  }
  
  # cbd <- calculate_adjp(cbd)
  # cbd <- calculate_performance(cbd, binary_truth = "truth")
  # cbp <- prepare_data_for_plot(cbd)
  # p <- plot_fdrtprcurve(cbp, pointsize=2, yaxisrange=c(0,1),
  #                       title=paste0(nSamples," samples"))
  # print(p)
  # return(list(cbd=cbd, p=p))
  return(cbd)
}

simple_auc <- function(TPR, FPR){
  # function from https://blog.revolutionanalytics.com/2016/11/calculating-auc.html
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

calcAUC <- function(cbd){
  truth <- iCOBRA::truth(cbd)
  deId <- truth$truth
  pvalues <- iCOBRA::pval(cbd)
  perf <- apply(pvalues, 2, function(pp){
    oop <- order(pp, decreasing=FALSE)
    tpr <- cumsum(deId[oop]) / sum(deId)
    fpr <- cumsum(!deId[oop]) / sum(!deId)
    fdr <- cumsum(!deId[oop]) / 1:length(oop)
    auROC <- simple_auc(tpr, fpr)
    auFDR <- simple_auc(tpr, fdr)
    return(c(auROC = auROC,
           auFDR = auFDR))
  })
  return(list(auroc = perf["auROC",],
              aufdptpr = perf["auFDR",]))
}

# Mock evaluation
mockEvaluation <- function(simCounts,
                           grpSim,
                           simGC,
                           simWidth,
                           grSim,
                           deId = NULL,
                           ruv = FALSE,
                           plot = TRUE, 
                           ...){
  # DA on mock
  evalMock <- evaluateSimulation(simCounts = simCounts,
                                 grpSim = grpSim,
                                 deId = NULL,
                                 simGC = simGC,
                                 grSim = grSim,
                                 simWidth = width(gr),
                                 ruv=FALSE, 
                                 ...)
  # set NA to 1
  pval(evalMock)[is.na(pval(evalMock))] <- 1
  
  ## calculate three metrics: FPR, p-value uniformity, p-value uniformity variability wrt GC.
  if(plot){
    rafalib::mypar(mfrow=c(4,2))
    for(kk in 1:ncol(pval(evalMock))) hist(pval(evalMock)[,kk], 
                                           breaks=20, main=colnames(pval(evalMock))[kk])
  }
  
  # FPR: DA peaks at 5% nominal level
  fpr <- apply(pval(evalMock), 2, function(x) mean(x <= 0.05))
  fpr
  
  # p-value uniformity
  helDist <- c()
  for(kk in 1:ncol(pval(evalMock))){
    curP <- pval(evalMock)[,kk]
    x <- cbind(curP, seq(0,1,length=length(curP)))
    helDist[kk] <- statip::hellinger(x=curP, y=seq(0,1,length=length(curP)))
  }
  names(helDist) <- colnames(pval(evalMock))
  
  # p-value uniformity variability wrt GC: note confounding with mean count...
  fprGCList <- list()
  distGCList <- list()
  for(kk in 1:ncol(pval(evalMock))){
    pvals <- pval(evalMock)[,kk]
    rafalib::mypar(mfrow=c(5,4))
    fprBin <- helDistBin <- c()
    for(bb in 1:nlevels(gcGroups)){
      pvalBin <- pvals[gcGroups == levels(gcGroups)[bb]]
      hist(pvalBin, breaks=seq(0,1,length=20))
      fprBin[bb] <- mean(pvalBin <= 0.05)
      helDistBin[bb] <- statip::hellinger(x=pvalBin, y=seq(0,1,length=length(pvalBin)))
    }
    fprGCList[[kk]] <- fprBin
    distGCList[[kk]] <- helDistBin
  }
  names(fprGCList) <- names(distGCList) <- colnames(pval(evalMock))
  helDistVar <- unlist(lapply(distGCList, var))
  if(plot){
    rafalib::mypar(mfrow=c(4,2))
    lapply(fprGCList, plot, type='b')
  }
  
  return(list(fpr = fpr, # FPR in mock
       distUnifPvalMock = helDist, # distance from uniform p-value in mock
       helDistVarGC = helDistVar # variability in Hellinger distance between p-value and uniform distribution across different GC-content bins in mock.
       ))
}


signalEvaluation <- function(counts,
                             grp,
                             gc,
                             gr,
                             nTotal = 12, ...){
  simData <- simulateData(counts = counts,
                          grp = grp,
                          nTotal = 12)
  
  evalSimResults <- evaluateSimulation(simCounts = simData$simCounts,
                                       grpSim = simData$grp,
                                       deId = simData$deId,
                                       simGC = gc,
                                       simWidth = width(gr),
                                       grSim = gr, ...)
  
  ## calculate two metrics: AUC under ROC or FDR-TPR, imbalance in DA peaks wrt GC-content
  ## AUROC
  aucRes <- calcAUC(evalSimResults)
  aucRes
  
  ## imbalance in DA peaks wrt GC content
  cbd <- iCOBRA::calculate_adjp(evalSimResults)
  padj <- iCOBRA::padj(cbd)
  daPeaksGC <- apply(padj, 2, function(x){
    gcContentPeaks[which(x <= 0.05)]
  })
  rafalib::mypar(mfrow=c(3,3))
  for(kk in 1:length(daPeaksGC)){
    if(length(daPeaksGC[[kk]]) == 0) next
    hist(daPeaksGC[[kk]], breaks = 20, xlim=c(0.2, 0.8), main=names(daPeaksGC)[kk])
  }
  # measure of how different the distributions are to the truth
  trueGC <- gcContentPeaks[simData$deId]
  hist(trueGC, breaks=20) # true distribution
  decdf <- function(x, true, sample)  ecdf(true)(x) - ecdf(sample)(x)
  grid <- seq(min(trueGC), max(trueGC), length=100)
  ecdfDiff <- unlist(lapply(daPeaksGC, function(x){
    if(length(x) == 0) return(NA)
    sum(abs(decdf(x = grid,
                  true = trueGC,
                  sample = x)))
  }))
  ecdfDiff
  
  return(list(DAGCDistDiff = ecdfDiff, # how different is derived DA GC distribution from true distribution
              auroc = aucRes$auroc,
              aufdptpr = aucRes$aufdptpr))
}
