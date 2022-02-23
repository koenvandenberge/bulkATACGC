.libPaths("/accounts/campus/koenvdberge/rpkgs")
library(gcapc, lib.loc = "/accounts/campus/koenvdberge/rpkgs")
library(rtracklayer, lib.loc = "/accounts/campus/koenvdberge/rpkgs")
library(BSgenome, lib.loc = "/accounts/campus/koenvdberge/rpkgs")
library(BSgenome.Hsapiens.UCSC.hg19, lib.loc = "/accounts/campus/koenvdberge/rpkgs")
bamFiles <- list.files("/accounts/campus/koenvdberge/gcapc/data",
                       pattern = "sorted.filtered.rmDup.bam",
                       full.names = TRUE)

peakList <- list()
for(ff in 2:length(bamFiles)){
  print(paste0("sample",ff))
  cov <- read5endCoverage(bamFiles[ff])
  cov
  selChrom <- paste0("chr",c(1:22,"X","Y"))
  cov$fwd <- cov$fwd[selChrom]
  cov$rev <- cov$rev[selChrom]
  bdw <- try(bindWidth(cov, range=c(50, 500), step=50))
  if(!is(bdw, "type-error")){
    gcb <- try(gcEffects(cov, bdw, sampling=c(0.05,1), 
                         emtrace=TRUE, plot=TRUE, model='poisson'))
    if(!is(gcb, "type-error")){
      peaks <- try(gcapcPeaks(cov, gcb, bdw, plot=TRUE, permute=5L))
      if(!is(peaks, "type-error")){
        peakList[[ff]] <- peaks
        saveRDS(peaks, file=paste0("/accounts/campus/koenvdberge/gcapc/gcapcPeakList",ff,".rds"))
      } else next
    } else next
  }
  else next
}
saveRDS(peakList, file="/accounts/campus/koenvdberge/gcapc/gcapcPeakListAll.rds")