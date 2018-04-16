getmapp = function (ref) {
  mapp <- rep(1, length(ref))
  seqlevelsStyle(ref)='UCSC'
  for(chr in unique(seqnames(ref))){
    message("Getting mappability for ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=GRanges(seqnames=chr, ranges=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index]))
    mapp.chr=rep(1, length(ref.chr))
    overlap=as.matrix(findOverlaps(ref.chr, mapp_hg19))
    for(i in unique(overlap[,1])){
      mapp.chr[i]=mean(mapp_hg19$mappability[overlap[which(overlap[,1]==i),2]])
    }
  mapp[chr.index]=mapp.chr
  }
  mapp
}