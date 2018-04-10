getmapp=function (ref) {
  mapp=rep(1, length(ref))
  for(chr in unique(seqnames(ref))){
    message("Getting mappability for chr ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index])
    if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
      chrtemp <- 23
    } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == 
             "chry") {
      chrtemp <- 24
    } else {
      chrtemp <- as.numeric(mapSeqlevels(as.character(chr), "NCBI")[1])
    }
    if (length(chrtemp) == 0) message("Chromosome cannot be found in NCBI Homo sapiens database!")
    mappability.chr <- mappability[[chrtemp]]
    mapp_ref.chr <- mapp_ref[[chrtemp]]
    for (k in 1:length(ref.chr)) {
      index <- countOverlaps(mapp_ref.chr, ref.chr[k])
      if (sum(index) > 0) {
        mapp[chr.index][k] <- mean(mappability.chr[as.logical(index)])
      }
    }
  }
  mapp=round(mapp,3)
  mapp
}