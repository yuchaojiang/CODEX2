getgc =function (ref) {
  gc=rep(NA,length(ref))
  for(chr in unique(seqnames(ref))){
    message("Getting GC content for ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index])
    if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
      chrtemp <- 23
    } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == 
               "chry") {
      chrtemp <- 24
    } else {
      chrtemp <- as.numeric(mapSeqlevels(as.character(chr), 
                                         "NCBI")[1])
    }
    if (length(chrtemp) == 0) message("Chromosome cannot be found in NCBI Homo sapiens database!")
    chrm <- unmasked(Hsapiens[[chrtemp]])
    seqs <- Views(chrm, ref.chr)
    af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
    gc[chr.index] <- round((af[, "G"] + af[, "C"]) * 100, 2)
  }
  gc
}
