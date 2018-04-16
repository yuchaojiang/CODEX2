getgc.mm10 =function (ref) {
  seqlevelsStyle(ref)='UCSC'
  gc=rep(NA,length(ref))
  for(chr in unique(seqnames(ref))){
    message("Getting GC content for ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index])
    if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
      chrtemp <- 20
    } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == "chry") {
      chrtemp <- 21
    } else if (chr == "M" | chr == "m" | chr == "chrM" | chr == {"chrm"}){
      chrtemp <- 22
    } else {
      chrtemp <- as.numeric(mapSeqlevels(as.character(chr), "NCBI")[1])
    }
    if (length(chrtemp) == 0) message("Chromosome cannot be found in NCBI database!")
    chrm <- unmasked(Mmusculus[[chrtemp]])
    seqs <- Views(chrm, ref.chr)
    af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
    gc[chr.index] <- round((af[, "G"] + af[, "C"]) * 100, 2)
  }
  gc
}
