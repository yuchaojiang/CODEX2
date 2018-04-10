normalize_null=function (Y_qc, gc_qc, K, N) 
{
  if (max(K) > ncol(Y_qc)) 
    stop("Number of latent Poisson factors K cannot exceed the number of \n         samples!")
  Ntotal <- N
  N <- round(N/median(N)*median(colSums(Y_qc)))
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
    k = K[ki]
    normObj.codex=normalize(Y_qc, gc_qc, K=k, N=Ntotal, message = FALSE)
    Yhat.codex=normObj.codex$Yhat[[1]]
    z.codex=log(pmax(Y_qc,1)/pmax(Yhat.codex,1))
    nonCNV.index=apply(z.codex,1,sd)<0.15
    
    normObj.codex.nonCNV=normalize(Y_qc[nonCNV.index,], gc_qc[nonCNV.index],
                                   K = k, N = Ntotal, message = TRUE)
    h.hat.codex2.nr=normObj.codex.nonCNV$h.hat[[1]]
    
    
    fGC.hat.codex2.nr=matrix(ncol=ncol(Y_qc),nrow=nrow(Y_qc),data=NA)
    fGC.hat.codex2.nr[nonCNV.index,]=normObj.codex.nonCNV$fGC.hat[[1]]
    
    for(t in 1:ncol(Y_qc)){
      spl <- smooth.spline(gc_qc[nonCNV.index],fGC.hat.codex2.nr[nonCNV.index,t])
      fGC.pred = predict(spl, gc_qc)$y
      fGC.pred[fGC.pred < 0] = min(fGC.pred[fGC.pred > 0])
      fGC.hat.codex2.nr[, t] <- fGC.pred
    }
    
    beta.hat.codex2.nr=rep(NA,nrow(Y_qc))
    g.hat.codex2.nr=matrix(nrow=nrow(Y_qc),ncol=k)
    L = log(Nmat * fGC.hat.codex2.nr)
    for (s in 1:nrow(Y_qc)) {
      temp = try(glm(formula = Y_qc[s, ] ~ 
                       h.hat.codex2.nr, offset = L[s, ],family = poisson)$coefficients, silent = TRUE)
      if (is.character(temp)) {
        temp = lm(log(pmax(Y_qc[s, ], 1)) ~ h.hat.codex2.nr,offset = L[s, ])$coefficients
      }
      g.hat.codex2.nr[s,] = temp[2:length(temp)]
      beta.hat.codex2.nr[s]=exp(temp[1])
    }
    betahatmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), 
                         data = beta.hat.codex2.nr, byrow = FALSE)
    Yhat.codex2.nr <- round(fGC.hat.codex2.nr * Nmat * betahatmat * 
                              exp(g.hat.codex2.nr %*% t(h.hat.codex2.nr)), 0)
    Yhat[[ki]] <- Yhat.codex2.nr
    fGC.hat[[ki]] <- signif(fGC.hat.codex2.nr, 3)
    beta.hat[[ki]] <- signif(beta.hat.codex2.nr, 3)
    h.hat[[ki]] <- signif(h.hat.codex2.nr, 3)
    g.hat[[ki]] <- signif(g.hat.codex2.nr, 3)
    AIC[ki] <- 2 * sum(Y_qc * log(pmax(Yhat[[ki]], 1)) - 
                         Yhat[[ki]]) - 2 * (length(g.hat[[ki]]) + length(h.hat[[ki]]))
    BIC[ki] <- 2 * sum(Y_qc * log(pmax(Yhat[[ki]], 1)) - 
                         Yhat[[ki]]) - (length(g.hat[[ki]]) + length(h.hat[[ki]])) * 
      log(length(Y_qc))
    RSS[ki] <- sum((Y_qc - Yhat[[ki]])^2/length(Y_qc))
    message("\tAIC", k, " = ", round(AIC[ki], 3))
    message("\tBIC", k, " = ", round(BIC[ki], 3))
    message("\tRSS", k, " = ", round(RSS[ki], 3), "\n")
  }
  list(Yhat = Yhat, fGC.hat = fGC.hat, beta.hat = beta.hat, 
       g.hat = g.hat, h.hat = h.hat, AIC = AIC, BIC = BIC, RSS = RSS, 
       K = K)
}