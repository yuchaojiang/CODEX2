######################################################
######################################################
####                                              ####
####          get coverage, gc, mapp              ####
####                                              ####
######################################################
######################################################

library(CODEX)
bamdir=as.matrix(read.table('bamlist'))
sampname=as.matrix(read.table('sampname.txt'))
projectname='melanoma_exon'
bedFile='melanoma_exon.bed'

chr=1
# get bam directories, read in bed file, get sample names
bambedObj=getbambed(bamdir=bamdir,
                    bedFile=bedFile,
                    sampname=sampname,
                    projectname=projectname,chr)
bamdir=bambedObj$bamdir; sampname=bambedObj$sampname; ref=bambedObj$ref; projectname=bambedObj$projectname;chr=bambedObj$chr
# get raw depth of coverage
coverageObj=getcoverage(bambedObj,mapqthres=20)
Y=coverageObj$Y; readlength=coverageObj$readlength
# get gc content
gc=getgc(chr,ref)
# get mappability
mapp=getmapp(chr,ref)

ref.all=bambedObj$ref
Y.all=coverageObj$Y
gc.all=gc
mapp.all=mapp
chr.all=rep(chr,length=length(mapp))

targ.chr <- unique(as.matrix(read.table(bedFile, sep = "\t")[,1]))

for(chr in 2:23){
  if(chr==23){chr='X'}
  if(!is.element(chr,targ.chr)) next
  # get bam directories, read in bed file, get sample names
  bambedObj=getbambed(bamdir=bamdir,
                      bedFile=bedFile,
                      sampname=sampname,
                      projectname=projectname,chr)
  bamdir=bambedObj$bamdir; sampname=bambedObj$sampname; ref=bambedObj$ref; projectname=bambedObj$projectname;chr=bambedObj$chr
  # get raw depth of coverage
  coverageObj=getcoverage(bambedObj,mapqthres=20)
  Y=coverageObj$Y; readlength=coverageObj$readlength
  # get gc content
  gc=getgc(chr,ref)
  # get mappability
  mapp=getmapp(chr,ref)
  
  ref.all=c(ref.all,bambedObj$ref)
  Y.all=rbind(Y.all,coverageObj$Y)
  gc.all=c(gc.all,gc)
  mapp.all=c(mapp.all,mapp)
  chr.all=c(chr.all,rep(chr,length=length(mapp)))
}

save.image(file=paste(projectname,'_','coverage','.rda',sep=''))

######################################################
######################################################
####                                              ####
####                   normalize                  ####
####                                              ####
######################################################
######################################################

library(CODEX)
projectname='melanoma_exon'
load(paste(projectname,'_coverage.rda',sep=''))

ref.all
length(chr.all)
dim(Y.all)
length(gc.all)
length(mapp.all)

gene.all=as.matrix(read.table('melanoma_exon.bed',head=F,sep='\t')[,4])
length(gene.all)

Y=Y.all
ref=ref.all
gc=gc.all
mapp=mapp.all
gene=gene.all

qcObj=qc(Y,sampname,chr,ref,mapp,gc,cov_thresh=c(60,4000),length_thresh=c(20,2000),mapp_thresh=0.9,gc_thresh=c(20,80))
Y_qc=qcObj$Y_qc; sampname_qc=qcObj$sampname_qc; gc_qc=qcObj$gc_qc; mapp_qc=qcObj$mapp_qc; ref_qc=qcObj$ref_qc ; qcmat=qcObj$qcmat
dim(Y_qc)
length(gc_qc)
length(mapp_qc)
gene_qc=gene[which(as.logical(qcmat[,4])==TRUE)]
chr_qc=chr.all[which(as.logical(qcmat[,4])==TRUE)]
length(gene_qc)
length(sampname_qc)

sampfilter=apply(Y_qc,2,median)>=50  # need to exclude 3 samples (exon capture failure)
sampname_qc=sampname_qc[sampfilter]
Y_qc=Y_qc[,sampfilter]
rm(qcObj)

# The code commented out below has been reported by other users to lead to error.
# This code does NOT affect the latter part to run but basically outputs QC results.
# qcmat=cbind(qcmat[,1:3],gene,qcmat[,4:12])
# colnames(qcmat)[4]='gene'
# qcmat[,1]=chr.all
# write.table(qcmat,file=paste(projectname,'_qcmat','.txt',sep=''),sep='\t',quote=F,row.names=F)

normObj=normalize2(Y_qc,gc_qc,K=1:15,normal_index=29:44)   # normal index needs to be changed.
Yhat=normObj$Yhat; AIC=normObj$AIC; BIC=normObj$BIC; RSS=normObj$RSS; K=normObj$K

save.image(file=paste(projectname,'_','normalize','.rda',sep=''))

######################################################
######################################################
####                                              ####
####                    segment                   ####
####                                              ####
######################################################
######################################################

library(CODEX)
projectname='melanoma_exon'
load(paste(projectname,'_normalize.rda',sep=''))

choiceofK(AIC,BIC,RSS,K,filename=paste(projectname,'_choiceofK','.pdf',sep=''))

dim(Y_qc)
length(ref_qc)
length(sampname_qc)
dim(Yhat[[2]])
length(chr_qc)

#plot(K, RSS, type = "b", xlab = "Number of latent variables")
#plot(K, AIC, type = "b", xlab = "Number of latent variables")
#plot(K, BIC, type = "b", xlab = "Number of latent variables")

optK=which.max(BIC)
source('segment_targeted.R')
finalcall=matrix(ncol=14)
colnames(finalcall)=c('sample_name','chr','gene','cnv',
                      'st_bp','ed_bp','length_kb',
                      'st_exon','ed_exon','raw_cov',
                      'norm_cov','copy_no','lratio',
                      'mBIC')
for(genei in unique(gene_qc)){
  cat('Segmenting gene',genei,'\n')
  geneindex=which(gene_qc==genei)
  yi=Y_qc[geneindex,]
  yhati=Yhat[[optK]][geneindex,]
  refi=ref_qc[geneindex]
  chri=chr_qc[geneindex][1]
  finalcalli=segment_targeted(yi, yhati, sampname_qc, refi, genei, chri, lmax=length(geneindex), mode='fraction') 
  finalcall=rbind(finalcall,finalcalli)
}
finalcall=finalcall[-1,]
cn=(as.numeric(as.matrix(finalcall[,'copy_no'])))
cn.filter=(cn<=1.7)|(cn>=2.3)
finalcall=finalcall[cn.filter,]

length_exon=as.numeric(finalcall[,'ed_exon'])-as.numeric(finalcall[,'st_exon'])+1
finalcall=cbind(finalcall[,1:7],length_exon,finalcall[,10:14])

write.table(finalcall, file = paste( projectname,'_', optK, '_CODEX_frac.txt',
                                     sep=''), sep='\t', quote=F, row.names=F)
save.image(file=paste(projectname,'_','finalcall','.rda',sep=''))
