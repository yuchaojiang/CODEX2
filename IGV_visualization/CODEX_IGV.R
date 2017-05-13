projectname='nbl'
chr='1'
finalcall.chr=read.table(paste(projectname,'_',chr,'_CODEX_frac2.txt',sep=''),head=T,sep='\t')
finalcall=finalcall.chr

for(chr in c(2:22,'X')){
  cat('chr =',chr,'\n')
  finalcall.chr=read.table(paste(projectname,'_',chr,'_CODEX_frac2.txt',sep=''),head=T,sep='\t')
  finalcall=rbind(finalcall,finalcall.chr)
}

finalcall=finalcall[(finalcall$copy_no<=1.5)|(finalcall$copy_no>=2.5),] # for 'fractional' mode filtering
#finalcall=finalcall[finalcall$length_kb>=100,] # filtering out short CNVs, not recommended

tumorsampname=as.matrix(read.table('tumorname_lookup.txt'))

for(i in 1:nrow(tumorsampname)){
  cat(i,'\t')
  finalcall.temp=finalcall[which(finalcall$sample_name==tumorsampname[i,1]),]
  sampname.temp=paste(rep(tumorsampname[i,3],nrow(finalcall.temp)),'.codex',sep='')
  output=cbind(sampname.temp,finalcall.temp$chr,
               finalcall.temp$st_bp,finalcall.temp$ed_bp,
               finalcall.temp$ed_exon-finalcall.temp$st_exon+1,
               signif(pmax(log(finalcall.temp$copy_no/2,2),-4),4))
  colnames(output)=c('Sample', 'Chromosome','Start','End','Num_Probes',  
                     'Segment_Mean')
  write.table(output,file=paste('codex_segments/',tumorsampname[i,3],'.codex.seg.txt',sep=''),
              sep='\t',quote=F,col.names=T,row.names=F)
}

