library(CODEX2)
library(BSgenome.Mmusculus.UCSC.mm10)

bed=read.table('mm10-Chr.bed',head=F)
bed=bed[,1:3] # three columns: chr, start position, end position
ref=GRanges(seqnames=bed[,1], ranges=IRanges(start=bed[,2], end=bed[,3]))
seqnames(ref)
seqlevelsStyle(ref)='UCSC'
unique(seqnames(ref))

L=100 # read length
genome=Mmusculus
chr='chr18'
ref.chr <- ref[which(seqnames(ref)==chr)]

computemapp<-function(ref.chr, L, genome, chr){
  genome.chr<-genome[[chr]]
  mapp.chr=rep(1,length(ref.chr))
  for(mappi in 1:length(mapp.chr)){
    cat(mappi,'\t')
    if(width(ref.chr)[mappi]>=L){
      dict = Views(genome.chr, start=seq((start(ref.chr))[mappi],(end(ref.chr))[mappi]-L+1,1), width=L)
    } else{
      dict = Views(genome.chr, start=start(ref.chr)[mappi]+round(width(ref.chr)[mappi]/2)-round(L/2), width=L)
    }
    if(sum(alphabetFrequency(dict, baseOnly=T)[,'other'])>0){
      mapp.chr[mappi]=0
    } else{
      pd = PDict(dict)
      ci=rep(0,length(pd))
      for(t in 1:length(genome)){
        res=matchPDict(pd,genome[[t]]); ci=ci+elementNROWS(res)
      }
      mapp.chr[mappi]=1/mean(ci)
    }
  }
  return(mapp.chr)
}

computemapp(ref.chr[1:3], L, genome, chr)
mapp.chr <- computemapp(ref.chr, L, genome, chr)
values(ref.chr) <- cbind(values(ref.chr), DataFrame(mapp.chr))  

save(ref.chr,file=paste("ref.",chr,".RData",sep=''))
