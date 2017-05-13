getcoverage <- function(bambedObj, mapqthres) {
    ref <- bambedObj$ref
    chr <- bambedObj$chr
    bamdir <- bambedObj$bamdir
    sampname <- bambedObj$sampname
    st <- start(ref)[1]
    ed <- end(ref)[length(ref)]
    Y <- matrix(NA, nrow = length(ref), ncol = nrow(sampname))
    readlength <- rep(NA, nrow(sampname))
    for (i in 1:nrow(sampname)) {
        bamurl <- bamdir[i]
        which <- RangesList(quack = IRanges(st - 10000, ed + 10000))
        names(which) <- as.character(chr)
        what <- c("pos", "mapq", "qwidth")
        flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE,
                            isNotPassingQualityControls = FALSE, 
                            isFirstMateRead = TRUE)
        param <- ScanBamParam(which = which, what = what, flag = flag)
        bam <- scanBam(bamurl, param = param)[[1]]
        mapqfilter <- (bam[["mapq"]] >= mapqthres)
        readlength[i] <- round(mean(bam[["qwidth"]]))
        if(is.nan(readlength[i])){
          flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE,
                              isNotPassingQualityControls = FALSE)
          param <- ScanBamParam(which = which, what = what, flag = flag)
          bam <- scanBam(bamurl, param = param)[[1]]
          mapqfilter <- (bam[["mapq"]] >= mapqthres)
          readlength[i] <- round(mean(bam[["qwidth"]]))
        }
        message("Getting coverage for sample ", sampname[i, 1], ": ", 
                "read length ", readlength[i], ".", sep = "")
        irang <- IRanges(bam[["pos"]][mapqfilter], width = 
                        bam[["qwidth"]][mapqfilter])
        Y[, i] <- countOverlaps(ref, irang)
    }
    list(Y = Y, readlength = readlength)
}