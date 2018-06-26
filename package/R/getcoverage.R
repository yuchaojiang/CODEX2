getcoverage=function (bambedObj, mapqthres) 
{
  ref <- bambedObj$ref
  bamdir <- bambedObj$bamdir
  sampname <- bambedObj$sampname
  Y <- matrix(NA, nrow = length(ref), ncol = length(sampname))
  for (i in 1:length(sampname)) {
    bamurl <- bamdir[i]
    what <- c("rname", "pos", "mapq", "qwidth")
    flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, 
                        isNotPassingQualityControls = FALSE, isFirstMateRead = TRUE)
    param <- ScanBamParam(what = what, flag = flag)
    bam <- scanBam(bamurl, param = param)[[1]]
    readlength.i=round(mean(bam$qwidth))
    if (is.nan(readlength.i)) {
      flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, 
                          isNotPassingQualityControls = FALSE)
      param <- ScanBamParam(what = what, flag = flag)
      bam <- scanBam(bamurl, param = param)[[1]]
      readlength.i=round(mean(bam$qwidth))
    }
    message("Getting coverage for sample ", sampname[i], 
            ": ", "read length ", readlength.i, ".", sep = "")
    bam.ref=GRanges(seqnames=bam$rname, ranges=IRanges(start=bam$pos, width=bam$qwidth))
    # remove reads with low mapping quality
    bam.ref=bam.ref[bam$mapq>=mapqthres]
    Y[, i] <- countOverlaps(ref, bam.ref)
  }
  rownames(Y)=paste(seqnames(ref),':',start(ref),'-',end(ref),sep='')
  colnames(Y)=sampname
  list(Y = Y)
}
