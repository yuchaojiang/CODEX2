normalize_codex2_nr = function (Y_qc, gc_qc, K, cnv_index) 
{
  if (max(K) > ncol(Y_qc)) 
    stop("Number of latent Poisson factors K cannot exceed the number of \n         samples!")
  N <- colSums(Y_qc)
  Nmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = N, 
                 byrow = TRUE)
  Yhat <- list(length = length(K))
  fGC.hat <- list(length = length(K))
  beta.hat <- list(length = length(K))
  g.hat <- list(length = length(K))
  h.hat <- list(length = length(K))
  AIC <- rep(NA, length = length(K))
  BIC <- rep(NA, length = length(K))
  RSS <- rep(NA, length = length(K))
  for (ki in 1:length(K)) {
    k <- K[ki]
    message("Computing normalization with k = ", k, " latent factors ...",sep="")
    message("\tRunning CODEX in regions without common CNVs ...")
    
    
    normObj=normalize_null(Y_qc[-cnv_index,],gc_qc[-cnv_index],K=k)
    
    Yhat.codex2.null=normObj$Yhat[[1]]
    beta.hat.codex2.null=normObj$beta.hat[[1]]
    g.hat.codex2.null=normObj$g.hat[[1]]
    h.hat.codex2.null=normObj$h.hat[[1]]
    fGC.hat.codex2.null=normObj$fGC.hat[[1]]
    
    message("\tRunning CODEX2 in regions with common CNVs ...")
    
    h.hat.codex2.nr=h.hat.codex2.null
    fGC.hat.codex2.nr=matrix(ncol=ncol(Y_qc),nrow=nrow(Y_qc))
    for(j in 1:ncol(fGC.hat.codex2.nr)){
      spl <- smooth.spline(gc_qc[-cnv_index], fGC.hat.codex2.null[,j])
      fGC.hat.codex2.nr[,j] <- predict(spl, gc_qc)$y
    }
    
    g.hat.codex2.nr=matrix(nrow=nrow(Y_qc),ncol=k)
    g.hat.codex2.nr[setdiff(1:nrow(Y_qc),cnv_index),]=g.hat.codex2.null
    beta.hat.codex2.nr=rep(NA,nrow(Y_qc))
    beta.hat.codex2.nr[setdiff(1:nrow(Y_qc),cnv_index)]=beta.hat.codex2.null
    
    for(cnvi in cnv_index){
      #cat('Processing common cnv index',cnvi,'...\n')
      fGC=fGC.hat.codex2.nr[cnvi,]
      Ytemp=Y_qc[cnvi,]
      h=h.hat.codex2.nr
      
      offset.temp=log(N)+log(fGC)
      # length(Ytemp); dim(h); length(offset.temp);length(offset.temp);length(N)
      
      # need to remove zeros (homozygous deletions) which are def not null
      h=h[Ytemp>5,,drop=F]
      fGC=fGC[Ytemp>5]
      offset.temp=offset.temp[Ytemp>5]
      Ntemp=N[Ytemp>5]
      Ytemp=Ytemp[Ytemp>5]
      
      EM.seed=matrix(ncol=3,nrow=12)
      colnames(EM.seed)=c('mu','cnvfreq','log-lik')
      EM.seed[,1]=c(rep(log(1/2),4),rep(log(3/2),4),rep(log(4/2),4))
      EM.seed[,2]=rep(1:4/10,3)
      
      for(emi in 1:nrow(EM.seed)){
        mu=EM.seed[emi,1]
        cnvfreq=EM.seed[emi,2]
        
        diff=1
        numiters=1
        
        beta0=0
        g=rep(0,k)
        
        diff.final=mu.final=pi.final=beta0.final=rep(NA,500)
        g.final=matrix(ncol=k, nrow=500)
        
        while(diff>0.0001 || numiters <= 30){
          Zhat=Estep(Ytemp,h,beta0,g,mu,pi,N=Ntemp,fGC)
          curM=Mstep(Ytemp,h,Zhat,N=Ntemp,fGC,mu)
          curbeta0=curM$beta0
          curu=curM$g
          curpi=curM$pi
          curmu=curM$mu
          
          diff=max(max(abs(beta0-curbeta0)),
                   max(abs(g-curu)),
                   max(abs(pi-curpi)),
                   max(abs(mu-curmu)))
          beta0=curbeta0
          g=curu
          pi=curpi
          mu=curmu
          #cat('Iteration:',numiters-1,'\t','diff =',diff,'\n')
          
          beta0.final[numiters]=beta0
          g.final[numiters,]=g
          pi.final[numiters]=pi
          mu.final[numiters]=mu
          diff.final[numiters]=diff
          
          if(diff==0) break
          if(numiters>=500) break
          if(numiters>=3 & diff>1) break
          numiters=numiters+1
        }
        beta0=beta0.final[max(which.min(diff.final))]
        g=g.final[max(which.min(diff.final)),]
        pi=pi.final[max(which.min(diff.final))]
        mu=mu.final[max(which.min(diff.final))]
        # sum(round(Zhat))/length(Zhat);pi
        
        Yfit=round(Ntemp*fGC*exp(beta0)* exp(g %*% t(h))*(exp(mu))^(round(Zhat)))
        EM.seed[emi,3]=sum(log(dpois(Ytemp,Yfit))) # log-likelihood of the fitted model.
      }
    best.seed=which.max(EM.seed[,3])
    cnvfreq=EM.seed[best.seed,2]
    mu=EM.seed[best.seed,1]
    
    
    diff=1
    numiters=1
    
    beta0=0
    g=rep(0,k)
    
    diff.final=mu.final=pi.final=beta0.final=rep(NA,500)
    g.final=matrix(ncol=k,nrow=500)
    
    while(diff>0.0001 || numiters <= 30){
      Zhat=Estep(Ytemp,h,beta0,g,mu,pi,N=Ntemp,fGC)
      curM=Mstep(Ytemp,h,Zhat,N=Ntemp,fGC,mu)
      curbeta0=curM$beta0
      curu=curM$g
      curpi=curM$pi
      curmu=curM$mu
      
      diff=max(max(abs(beta0-curbeta0)),
               max(abs(g-curu)),
               max(abs(pi-curpi)),
               max(abs(mu-curmu)))
      beta0=curbeta0
      g=curu
      pi=curpi
      mu=curmu
      #cat('Iteration:',numiters-1,'\t','diff =',diff,'\n')
      
      beta0.final[numiters]=beta0
      g.final[numiters,]=g
      pi.final[numiters]=pi
      mu.final[numiters]=mu
      diff.final[numiters]=diff
      
      if(diff==0) break
      if(numiters>=500) break
      if(numiters>=3 & diff>1) break
      numiters=numiters+1
    }
    beta0=beta0.final[max(which.min(diff.final))]
    g=g.final[max(which.min(diff.final)),]
    pi=pi.final[max(which.min(diff.final))]
    mu=mu.final[max(which.min(diff.final))]
    
    beta.hat.codex2.nr[cnvi]=exp(beta0)
    g.hat.codex2.nr[cnvi,]=g
    }
    betahatmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = beta.hat.codex2.nr, byrow = FALSE)
    Yhat.codex2.nr <- round(fGC.hat.codex2.nr * Nmat * betahatmat * exp(g.hat.codex2.nr %*% t(h.hat.codex2.nr)), 0)
    
    Yhat[[ki]] <- Yhat.codex2.nr
    fGC.hat[[ki]] <- signif(fGC.hat.codex2.nr,3)
    beta.hat[[ki]] <- signif(beta.hat.codex2.nr,3)
    h.hat[[ki]] <- signif(h.hat.codex2.nr,3)
    g.hat[[ki]] <- signif(g.hat.codex2.nr,3)
    AIC[ki] <- 2 * sum(Y_qc * log(pmax(Yhat[[ki]], 1)) - 
                         Yhat[[ki]]) - 2 * (length(g.hat[[ki]]) + length(h.hat[[ki]]))
    BIC[ki] <- 2 * sum(Y_qc * log(pmax(Yhat[[ki]], 1)) - 
                         Yhat[[ki]]) - (length(g.hat[[ki]]) + length(h.hat[[ki]])) * log(length(Y_qc))
    RSS[ki] <- sum((Y_qc - Yhat[[ki]])^2/length(Y_qc))
    message("\tAIC", k, " = ", round(AIC[ki], 3))
    message("\tBIC", k, " = ", round(BIC[ki], 3))
    message("\tRSS", k, " = ", round(RSS[ki], 3), "\n")
  }
  list(Yhat = Yhat, fGC.hat = fGC.hat, beta.hat = beta.hat, g.hat = g.hat,
       h.hat = h.hat, AIC = AIC, BIC = BIC, RSS = RSS, K = K)
}
