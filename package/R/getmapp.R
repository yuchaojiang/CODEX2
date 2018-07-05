getmapp = function (ref, genome = NULL) {
  if(is.null(genome)){genome = BSgenome.Hsapiens.UCSC.hg19}
  if(genome@provider_version == 'hg19'){mapp_gref = mapp_hg19}
  if(genome@provider_version == 'hg38'){mapp_gref = mapp_hg38}
  mapp <- rep(1, length(ref))
  seqlevelsStyle(ref)='UCSC'
  for(chr in unique(seqnames(ref))){
    message("Getting mappability for ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=GRanges(seqnames=chr, ranges=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index]))
    mapp.chr=rep(1, length(ref.chr))
    overlap=as.matrix(findOverlaps(ref.chr, mapp_gref))
    for(i in unique(overlap[,1])){
      index.temp=overlap[which(overlap[,1]==i),2]
      mapp.chr[i]=sum((mapp_gref$score[index.temp])*(width(mapp_gref)[index.temp]))/
        sum(width(mapp_gref)[index.temp])
      #mapp.chr[i]=mean(mapp_gref$score[index.temp])
    }
  mapp[chr.index]=mapp.chr
  }
  mapp
}