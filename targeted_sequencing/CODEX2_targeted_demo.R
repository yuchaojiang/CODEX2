library(CODEX2)
##########################################################
# Initialization
##########################################################
dirPath=getwd()
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
sampname <- substr(bamFile,1,9)
bedFile <- file.path(dirPath, "scWGA500kb.bed")
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "CODEX2_demo")
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname

##########################################################
# Getting GC content and mappability
##########################################################
gc <- getgc(ref)
mapp <- getmapp(ref)

##########################################################
# Getting gene names, needed for targeted sequencing, here generating gene names in silico
##########################################################
gene=rep(NA,length(ref))
for(chr in as.matrix(unique(seqnames(ref)))){
  chr.index=which(seqnames(ref)==chr)
  gene[chr.index]=paste(chr,'_gene_',ceiling(chr.index/30),sep='')
}
values(ref) <- cbind(values(ref), DataFrame(gc, mapp, gene))  

##########################################################
# Getting depth of coverage
##########################################################
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
write.csv(Y, file = paste(projectname, '_coverage.csv', sep=''), quote = FALSE)
head(Y[,1:5])

##########################################################
# Quality control
##########################################################
qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, Inf),
            length_thresh = c(20, Inf), mapp_thresh = 0.9,
            gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc
write.table(qcmat, file = paste(projectname, '_qcmat', '.txt', sep=''),
            sep = '\t', quote = FALSE, row.names = FALSE)

##########################################################
# Estimating library size factor for each sample
##########################################################
Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){prod(x)^(1/length(x))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)
plot(N, apply(Y,2,sum), xlab='Estimated library size factor', ylab='Total sum of reads')

##########################################################
# Genome-wide normalization using normalize_null
##########################################################
# If there are negative control samples, use normalize_codex2_ns()
# If there are negative control regions, use normalize_codex2_nr()
normObj.null <- normalize_null(Y_qc = Y_qc,
                               gc_qc = gc_qc,
                               K = 1:4, N = N)
Yhat <- normObj.null$Yhat
AIC <- normObj.null$AIC; BIC <- normObj.null$BIC
RSS <- normObj.null$RSS

##########################################################
# Number of latent factors
##########################################################
choiceofK(AIC, BIC, RSS, K = 1:4 , filename = "codex2_null_choiceofK.pdf")
par(mfrow = c(1, 3))
plot(1:4, RSS, type = "b", xlab = "Number of latent variables", pch=20)
plot(1:4, AIC, type = "b", xlab = "Number of latent variables", pch=20)
plot(1:4, BIC, type = "b", xlab = "Number of latent variables", pch=20)
par(mfrow = c(1,1))

##########################################################
# CBS segmentation per chromosome: optimal for WGS and WES
##########################################################
chr='chr1'
chr.index=which(seqnames(ref_qc)==chr)
finalcall.CBS <- segmentCBS(Y_qc[chr.index,],
                            Yhat, optK = which.max(BIC),
                            K = 1:5,
                            sampname_qc = sampname_qc,
                            ref_qc = ranges(ref_qc)[chr.index],
                            chr = chr, lmax = 400, mode = "integer")

##########################################################
# CBS segmentation per gene: optinmal for targeted seq
##########################################################
source('segment_targeted.R')
# Available at: https://github.com/yuchaojiang/CODEX2/blob/master/targeted_sequencing/segment_targeted.R
optK=which.max(BIC)
finalcall=matrix(ncol=14,nrow=0)
colnames(finalcall)=c('sample_name','chr','gene','cnv',
                      'st_bp','ed_bp','length_kb',
                      'st_exon','ed_exon','raw_cov',
                      'norm_cov','copy_no','lratio',
                      'mBIC')
for(genei in unique(ref_qc$gene)){
  cat('Segmenting gene',genei,'\n')
  geneindex=which(ref_qc$gene==genei)
  yi=Y_qc[geneindex,, drop=FALSE]
  yhati=Yhat[[optK]][geneindex,, drop=FALSE]
  refi=ref_qc[geneindex]
  finalcalli=segment_targeted(yi, yhati, sampname_qc, refi, genei, lmax=length(geneindex), mode='fraction') 
  finalcall=rbind(finalcall,finalcalli)
}
cn=(as.numeric(as.matrix(finalcall[,'copy_no'])))
cn.filter=(cn<=1.7)|(cn>=2.3) # removing calls with fractional copy numbers close to 2 (for heterogeneous cancer samples)
finalcall=finalcall[cn.filter,]
length_exon=as.numeric(finalcall[,'ed_exon'])-as.numeric(finalcall[,'st_exon'])+1
finalcall=cbind(finalcall[,1:7],length_exon,finalcall[,10:14])

write.table(finalcall, file = paste( projectname,'_', optK, '_CODEX_frac.txt',
                                     sep=''), sep='\t', quote=F, row.names=F)
