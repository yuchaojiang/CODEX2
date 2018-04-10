getcoverage=function (bambedObj, mapqthres, chr) 
{
  ref <- bambedObj$ref
  bamdir <- bambedObj$bamdir
  sampname <- bambedObj$sampname
  message("Getting coverage for chr ", chr, sep = "")
  chr.index=which(as.matrix(seqnames(ref))==chr)
  ref.chr=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index])
  Y <- matrix(NA, nrow = length(ref.chr), ncol = length(sampname))
  for (i in 1:length(sampname)) {
    bamurl <- bamdir[i]
    st <- start(ref.chr)[1]
    ed <- end(ref.chr)[length(ref.chr)]
    which <- RangesList(quack = IRanges(st - 10000, ed + 10000))
    names(which) <- as.character(chr)
    what <- c("pos", "mapq", "qwidth")
    flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, 
                        isNotPassingQualityControls = FALSE, isFirstMateRead = TRUE)
    param <- ScanBamParam(which = which, what = what, flag = flag)
    bam <- scanBam(bamurl, param = param)[[1]]
    mapqfilter <- (bam[["mapq"]] >= mapqthres)
    readlength.i=round(mean(bam[["qwidth"]]))
    if (is.nan(readlength.i)) {
      flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, 
                          isNotPassingQualityControls = FALSE)
      param <- ScanBamParam(which = which, what = what, 
                            flag = flag)
      bam <- scanBam(bamurl, param = param)[[1]]
      mapqfilter <- (bam[["mapq"]] >= mapqthres)
    }
    message("\t...sample ", sampname[i], 
            ": ", "read length ", readlength.i, ".", sep = "")
    irang <- IRanges(bam[["pos"]][mapqfilter], width = bam[["qwidth"]][mapqfilter])
    Y[, i] <- countOverlaps(ref.chr, irang)
  }
  rownames(Y)=paste(chr,':',start(ref.chr),'-',end(ref.chr),sep='')
  colnames(Y)=sampname
  list(Y = Y)
}