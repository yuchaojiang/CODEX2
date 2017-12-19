normalize_codex2_ns = function (Y_qc, gc_qc, K, norm_index) 
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
    message("\tRunning CODEX on normal samples ...")
    normObj=normalize_null(Y_qc[,norm_index],gc_qc,K=k)
    
    Yhat.codex2.null=normObj$Yhat[[1]]
    beta.hat.codex2.null=normObj$beta.hat[[1]]
    g.hat.codex2.null=normObj$g.hat[[1]]
    h.hat.codex2.null=normObj$h.hat[[1]]
    fGC.hat.codex2.null=normObj$fGC.hat[[1]]
    
    beta.hat.codex2.ns=beta.hat.codex2.null
    g.hat.codex2.ns=g.hat.codex2.null
    h.hat.codex2.ns=matrix(nrow=ncol(Y_qc),ncol=k)
    h.hat.codex2.ns[norm_index,] =  h.hat.codex2.null
    
    fGC.hat.codex2.ns=matrix(ncol=ncol(Y_qc),nrow=nrow(Y_qc),data=NA)
    fGC.hat.codex2.ns[,norm_index]=fGC.hat.codex2.null
    
    message("\tRunning CODEX2 in case-control setting ...")
    
    for(j in which(is.na(h.hat.codex2.ns[,1]))){
      #cat(j,'\n')
      hj.temp=rep(0,k)
      diff=1
      while(diff > 0.0001){
        z=Y_qc[,j]/Nmat[,j]/beta.hat.codex2.ns/exp(g.hat.codex2.ns%*%as.matrix(hj.temp))
        spl <- smooth.spline(gc_qc, z)
        fGC.temp <- predict(spl, gc_qc)$y
        
        L <- log(Nmat[,j] * fGC.temp * beta.hat.codex2.ns)
        
        hj.ns <- glm(formula = Y_qc[, j] ~ g.hat.codex2.ns - 1, offset = L, family = poisson)$coefficients
        diff=max(abs(hj.ns-hj.temp))
        hj.temp=hj.ns
      }
      fGC.hat.codex2.ns[,j]=fGC.temp
      h.hat.codex2.ns[j,]=hj.temp
    }
    
    betahatmat.new <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = beta.hat.codex2.ns, byrow = FALSE)
    Yhat.codex2.ns <- round(fGC.hat.codex2.ns * Nmat * betahatmat.new * exp(g.hat.codex2.ns %*% t(h.hat.codex2.ns)), 0)
    
    Yhat[[ki]] <- Yhat.codex2.ns
    fGC.hat[[ki]] <- signif(fGC.hat.codex2.ns,3)
    beta.hat[[ki]] <- signif(beta.hat.codex2.ns,3)
    h.hat[[ki]] <- signif(h.hat.codex2.ns,3)
    g.hat[[ki]] <- signif(g.hat.codex2.ns,3)
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
