normalize = function (Y_qc, gc_qc, K, N, message = NULL) 
{
  if (max(K) > ncol(Y_qc)) 
    stop("Number of latent Poisson factors K cannot exceed the number of \n         samples!")
  if (is.null(message)){message = TRUE}
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
    k <- K[ki]
    if(message){
      message("Computing normalization with k = ", k, " latent factors ...",sep="") 
    }
    maxiter <- 10
    maxhiter <- 50
    BHTHRESH <- 1e-04
    HHTHRESH <- 1e-05
    iter <- 1
    fhat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), 
                   data = 0)
    fhatnew <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc))
    betahat <- rep(1, nrow(Y_qc))
    betahatmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), 
                         data = betahat, byrow = FALSE)
    ghat <- matrix(0, nrow = nrow(Y_qc), ncol = k)
    hhat <- matrix(0, nrow = ncol(Y_qc), ncol = k)
    bhdiff <- rep(Inf, maxiter)
    fhdiff <- rep(Inf, maxiter)
    betahatlist <- list(length = maxiter)
    fhatlist <- list(length = maxiter)
    ghatlist <- list(length = maxiter)
    hhatlist <- list(length = maxiter)
    while (iter <= maxiter) {
      gcfit <- Y_qc/Nmat/betahatmat/exp(ghat %*% t(hhat))
      fhatnew <- apply(gcfit, 2, function(z) {
        spl <- smooth.spline(gc_qc, z)
        temp <- predict(spl, gc_qc)$y
        temp[temp <= 0] <- min(temp[temp > 0])
        temp
      })
      betahatnew <- apply(Y_qc/(fhatnew * Nmat * exp(ghat %*% 
                                                       t(hhat))), 1, median)
      bhdiff[iter] <- sum((betahatnew - betahat)^2)/length(betahat)
      fhdiff[iter] <- sum((fhatnew - fhat)^2)/length(fhat)
      if (fhdiff[iter] > min(fhdiff)) 
        break
      if(message){
        message("Iteration ", iter, "\t", "beta diff =", 
                signif(bhdiff[iter], 3), "\t", "f(GC) diff =", 
                signif(fhdiff[iter], 3)) 
      }
      fhat <- fhatnew
      betahat <- betahatnew
      betahatmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), 
                           data = betahat, byrow = FALSE)
      L <- log(Nmat * fhat * betahatmat)
      logmat <- log(pmax(Y_qc, 1)) - L
      logmat <- logmat - matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), 
                                data = apply(logmat, 1, mean), byrow = FALSE)
      #### Threshold logmat so we are not thrown off by the large values.
      # logmat=threshold(logmat,-3*sd(logmat),3*sd(logmat))
      svdres=svd(logmat, nu = k, nv = k)
      hhat <- svdres$v
      
      hhatnew <- hhat
      hiter <- 1
      hhdiff <- rep(Inf, maxhiter)
      while (hiter <= maxhiter) {
        for (s in 1:nrow(Y_qc)) {
          ghat[s, ] <- glm(formula = Y_qc[s, ] ~ hhat - 
                             1, offset = L[s, ], family = poisson)$coefficients
        }
        for (t in 1:ncol(Y_qc)) {
          hhatnew[t, ] <- glm(formula = Y_qc[, t] ~ ghat - 
                                1, offset = L[, t], family = poisson)$coefficients
        }
        gh <- ghat %*% t(hhatnew)
        gh <- gh - matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), 
                          data = apply(gh, 1, mean), byrow = FALSE)
        hhatnew <- svd(gh, nu = k, nv = k)$v
        hhdiff[hiter] <- sum((hhatnew - hhat)^2)/length(hhat)
        hhat <- hhatnew
        if (hhdiff[hiter] < HHTHRESH) 
          break
        if (hiter > 10 & (rank(hhdiff))[hiter] <= 3) 
          break
        hiter <- hiter + 1
      }
      fhatlist[[iter]] <- fhat
      betahatlist[[iter]] <- betahat
      ghatlist[[iter]] <- ghat
      hhatlist[[iter]] <- hhat
      if (bhdiff[iter] < BHTHRESH) 
        break
      if (iter > 5 & bhdiff[iter] > 1) 
        break
      iter <- iter + 1
    }
    optIter <- which.min(fhdiff)
    if(message){
      message(paste("Stop at Iteration ", optIter, ".", sep = ""))
    }
    fhat <- fhatlist[[optIter]]
    betahat <- betahatlist[[optIter]]
    ghat <- ghatlist[[optIter]]
    hhat <- hhatlist[[optIter]]
    betahatmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), 
                         data = betahat, byrow = FALSE)
    Yhat[[ki]] <- round(fhat * Nmat * betahatmat * exp(ghat %*% 
                                                         t(hhat)), 0)
    fGC.hat[[ki]] <- signif(fhat,3)
    beta.hat[[ki]] <- signif(betahat,3)
    h.hat[[ki]] <- signif(hhat,3)
    g.hat[[ki]] <- signif(ghat,3)
    AIC[ki] <- 2 * sum(Y_qc * log(pmax(Yhat[[ki]], 1)) - 
                         Yhat[[ki]]) - 2 * (length(ghat) + length(hhat))
    BIC[ki] <- 2 * sum(Y_qc * log(pmax(Yhat[[ki]], 1)) - 
                         Yhat[[ki]]) - (length(ghat) + length(hhat)) * log(length(Y_qc))
    RSS[ki] <- sum((Y_qc - Yhat[[ki]])^2/length(Y_qc))
    if(message){
      message("AIC", k, " = ", round(AIC[ki], 3))
      message("BIC", k, " = ", round(BIC[ki], 3))
      message("RSS", k, " = ", round(RSS[ki], 3), "\n")
    }
  }
  list(Yhat = Yhat, fGC.hat = fGC.hat, beta.hat = beta.hat, g.hat = g.hat,
       h.hat = h.hat, AIC = AIC, BIC = BIC, RSS = RSS, K = K)
}