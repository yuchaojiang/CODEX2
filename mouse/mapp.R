library(CODEX)
library(WES.1KG.WUGSC) # Load Toy data from the 1000 Genomes Project.
dirPath <- system.file("extdata", package = "WES.1KG.WUGSC")
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
sampname <- as.matrix(read.table(file.path(dirPath, "sampname")))
bedFile <- file.path(dirPath, "chr22_400_to_500.bed")
chr <- 22
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "CODEX_demo", chr)
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname; chr <- bambedObj$chr

coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y; readlength <- coverageObj$readlength

mapp_codex <- getmapp(chr, ref) # This is the mappability by CODEX by default
                                # (pre-computed and stored as part of the package
                                # to save computation time)
mapp=rep(1,length(ref))
L=median(readlength) # compute the mappability in 100 base windows, the read length is 100 for this WES toy dataset

library("BSgenome.Hsapiens.UCSC.hg19")
hum <- Hsapiens$chr22

for (i in 1:length(mapp)){
  if(width(ref)[i]>=100){
    dict = Views(hum, start=seq((start(ref))[i],(end(ref))[i]-L+1,1), width=L)
  } else{
    dict = Views(hum, start=start(ref)[i]+round(width(ref)[i]/2)-round(L/2), width=L)
  }
  #pd = PDict(dict, tb.start=floor(L/2)-15, tb.end=floor(L/2)+15)
  pd = PDict(dict)
  ci=rep(0,length(pd)); res=matchPDict(pd,Hsapiens$chr1)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr2)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr3)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr4)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr5)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr6)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr7)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr8)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr9)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr10)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr11)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr12)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr13)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr14)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr15)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr16)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr17)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr18)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr19)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr20)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr21)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chr22)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chrX)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chrY)
  ci=ci+elementNROWS(res); res=matchPDict(pd,Hsapiens$chrM)
  ci=ci+elementNROWS(res)
  mapp[i]=1/mean(ci)
  if (i%%10==0){cat(i,"\n")}
}
mapp=GRanges(seqnames=chr, ranges=ref, mapp=mapp)
mapp$mapp
mapp@ranges
save(mapp,file=paste("mapp_",chr,".RData",sep=''))
