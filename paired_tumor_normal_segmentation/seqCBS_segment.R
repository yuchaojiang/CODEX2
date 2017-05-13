install.packages('clue')
# install seqCBS from the tarball (don't install directly from CRAN, not updated)
install.packages('seqCBS_1.3.1.tar.gz',repos=NULL,type='source')

library(seqCBS)

optK=which.max(BIC)
dim(Y_qc)
dim(Yhat[[optK]])
length(ref_qc)
length(sampname_qc)

# If you samples come in an order of normal1, tumor1, normal2, tumor2, 
# normal3, tumor3, …, then you need to use the command below. Otherwise,
# make modifications to the tumor and normal sample index.

N0 = dim(Y_qc)[1]
nsamp = length(sampname_qc)/2
relCN = matrix(0,N0,nsamp)
tumorid = 2*(1:nsamp)
normalid = tumorid -1

for (j in 1:nsamp){
  cat('Processing matched pair',j,'\n')
  N = N0
  id1 = normalid[j]
  id2 = tumorid[j]
  normal = Y_qc[,id1]/Yhat[[optK]][,id1]*median(Yhat[[optK]][,id1])
  tumor = Y_qc[,id2]/Yhat[[optK]][,id2]*median(Yhat[[optK]][,id2])
  low = which(Yhat[[optK]][,id1]<10 | Yhat[[optK]][,id2]<10)
  if (length(low)>0){
    normal = normal[-low]
    tumor = tumor[-low]
    N = N0-length(low)
  }
  result = ScanCBS(rep(1:N,round(tumor)), rep(1:N,round(normal)))
  tau = result$tauHat
  m = length(tau)
  cn0 = result$relCN
  if (length(low)>0){
    index = (1:N0)[-low]
    tau = index[tau]
    if (tau[1]>1){
      tau[1] = 1
    }
    if (tau[m]<N0){
      tau[m] = N0
    }
  }
  dif = diff(tau)
  removeid = which(dif==0)
  if (length(removeid)>0){
    tau = tau[-removeid]
    cn0 = cn0[-removeid]
    dif = dif[-removeid]
  }
  dif[1] = dif[1]+1
  relCN[,j] = rep(cn0,dif)
}
save(relCN, file=paste("relCN_", chr, ".rda", sep=""))

library(fields)
png(paste("relCN_", chr, ".png", sep=""), width=1600, height=700, pointsize=25)
par(mar=c(4,4,1,4))
image.plot(relCN, zlim=c(0,2))
dev.off()


