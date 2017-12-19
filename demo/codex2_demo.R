
# We in silico spiked in CNVs spanning exon 1580 - 1620 with a population frequency 40%.
# There are altogether 90 samples, 36 of which have the heterozygous deletion.

###################################
# Negative control samples ...
###################################

library(CODEX2)

# Y_qc and gc_qc can be obtained from the sequencing bam files using CODEX.
# See https://github.com/yuchaojiang/CODEX for details/demo codes/vignettes.
Y_qc = Y_qc_codex2
gc_qc = gc_qc_codex2
# For the case-control scenario, the normal sample index is known (samples without spike-in signals).
norm_index = norm_index_codex2

normObj=normalize_codex2_ns(Y_qc = Y_qc, gc_qc = gc_qc, 
                            K = 1:3, norm_index = norm_index)
Yhat.ns=normObj$Yhat; fGC.hat.ns=normObj$fGC.hat;
beta.hat.ns=normObj$beta.hat; g.hat.ns=normObj$g.hat; h.hat.ns=normObj$h.hat
AIC.ns=normObj$AIC; BIC.ns=normObj$BIC; RSS.ns=normObj$RSS

finalcall.codex2.ns <- segment(Y_qc, Yhat.ns, optK = which.max(BIC.ns),
                               K = 1:3, sampname_qc = paste('sample',1:ncol(Y_qc),sep=''),
                               ref_qc = IRanges(start=1:nrow(Y_qc)*100, end=1:nrow(Y_qc)*100+50),
                               chr = 18, lmax = 200, mode = "integer")
nrow(finalcall.codex2.ns[finalcall.codex2.ns[,'st_exon']=='1580',])
head(finalcall.codex2.ns[finalcall.codex2.ns[,'st_exon']=='1580',])


###################################
# Negative control regions ...
###################################

library(CODEX2)

Y_qc = Y_qc_codex2
gc_qc = gc_qc_codex2

# We can empirically identify common CNV regions by a first-pass CODEX run
# For exons residing in common CNV regions, the s.d. of normalized z-scores across all samples
# will be large.
normObj=normalize_null(Y_qc = Y_qc, gc_qc = gc_qc, K = 1:3)
z.codex = log(Y_qc/normObj$Yhat[[which.max(normObj$BIC)]])
plot(1:nrow(z.codex),apply(z.codex,1,sd))
which(apply(z.codex,1,sd)>=0.25) 
cnv_index1=which(apply(z.codex,1,sd)>=0.25) 
# This can also be provided by the user as known,
# e.g., from existing database (DGV or dbVar) or knowledge (tumor supressors or oncogenes).
cnv_index2=1580:1620 

normObj=normalize_codex2_nr(Y_qc = Y_qc, gc_qc = gc_qc, 
                            K = 1:3, cnv_index = cnv_index1)
Yhat.nr=normObj$Yhat; fGC.hat.nr=normObj$fGC.hat;
beta.hat.nr=normObj$beta.hat; g.hat.nr=normObj$g.hat; h.hat.nr=normObj$h.hat
AIC.nr=normObj$AIC; BIC.nr=normObj$BIC; RSS.nr=normObj$RSS

finalcall.codex2.nr <- segment(Y_qc, Yhat.nr, optK = which.max(BIC.nr),
                               K = 1:3, sampname_qc = paste('sample',1:ncol(Y_qc),sep=''),
                               ref_qc = IRanges(start=1:nrow(Y_qc)*100, end=1:nrow(Y_qc)*100+50),
                               chr = 18, lmax = 200, mode = "integer")
nrow(finalcall.codex2.nr[finalcall.codex2.nr[,'st_exon']=='1580',])
head(finalcall.codex2.nr[finalcall.codex2.nr[,'st_exon']=='1580',])
