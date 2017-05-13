segment_targeted <- function(yi, yhati, sampname_qc, refi, genei, chri, lmax, mode) {
  finalcall <- matrix(ncol = 10)
  lmax <- lmax - 1
  if(is.vector(yi)){
    yi=t(yi)
    yhati=t(yhati)
  }
  for (sampno in 1:ncol(yi)) {
    #message("Segmenting sample ", sampno, ": ", sampname_qc[sampno], ".")
    y <- yi[, sampno]
    yhat <- yhati[, sampno]
    num <- length(y)
    y <- c(y, rep(0, lmax))
    yhat <- c(yhat, rep(0, lmax))
    i <- rep(1:num, rep(lmax, num))
    j <- rep(1:lmax, num) + i
    yact <- rep(0, length(i))
    lambda <- rep(0, length(i))
    for (k in 1:num) {
      yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + 
                                                                lmax)])[-1]
      lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + 
                                                                     lmax)])[-1]
    }
    i <- i[j <= num]
    yact <- yact[j <= num]
    lambda <- lambda[j <= num]
    j <- j[j <= num]
    yact[lambda<20]=20
    lambda[lambda<20]=20
    if (mode == "integer") {
      chat <- round(2 * (yact/lambda))
    } else if (mode == "fraction") {
      chat <- 2 * (yact/lambda)
    }
    lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
    # chat[chat > 5] <- 5
    if (sum(lratio > 0) > 0) {
      if (sum(lratio > 0) >= 2) {
        finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio > 
                                                                0, ]
        finalmat <- finalmat[order(-finalmat[, 6]), ]
        s <- 1
        while (s <= (nrow(finalmat))) {
          rowstart <- finalmat[s, 1]
          rowend <- finalmat[s, 2]
          rowsel <- (finalmat[, 1] <= rowend & finalmat[, 2] >= 
                       rowstart)
          rowsel[s] <- FALSE
          finalmat <- finalmat[!rowsel, ]
          if (is.vector(finalmat)) {
            finalmat <- t(as.matrix(finalmat))
          }
          s <- s + 1
        }
      }
      if (sum(lratio > 0) == 1) {
        finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio > 
                                                                0, ]
        finalmat <- t(as.matrix(finalmat))
      }
      finalmat <- round(finalmat, digits = 3)
      loglikeij <- cumsum(finalmat[, 6])
      mBIC <- rep(NA, length(loglikeij))
      for (s in 1:nrow(finalmat)) {
        tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), 1, num)))
        P <- length(tau) - 2
        mbic <- loglikeij[s]
        mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] - 
                                       tau[1:(length(tau) - 1)]))
        mbic <- mbic + (0.5 - P) * log(num)
        mBIC[s] <- mbic
      }
      mBIC <- round(mBIC, digits = 3)
      if (mBIC[1] > 0) {
        finalmat <- cbind(rep(sampname_qc[sampno], nrow(finalmat)), 
                          rep(chri, nrow(finalmat)),rep(genei, nrow(finalmat)), finalmat)
        finalmat <- (cbind(finalmat, mBIC)[1:which.max(mBIC), ])
        finalcall <- rbind(finalcall, finalmat)
      }
    }
  }
  finalcall <- finalcall[-1, ]
  if (is.vector(finalcall)) {
    finalcall <- t(as.matrix(finalcall))
  }
  st <- start(refi)[as.numeric(finalcall[, 4])]
  ed <- end(refi)[as.numeric(finalcall[, 5])]
  cnvtype <- rep(NA, length(st))
  cnvtype[as.numeric(finalcall[, 8]) < 2] <- "del"
  cnvtype[as.numeric(finalcall[, 8]) > 2] <- "dup"
  if (nrow(finalcall) == 1) {
    finalcall <- t(as.matrix(c(finalcall[, 1:3], cnvtype, st, ed, (ed - 
                                                                     st + 1)/1000, finalcall[, 4:10])))
  } else {
    finalcall <- cbind(finalcall[, 1:3], cnvtype, st, ed, (ed - st + 
                                                             1)/1000, finalcall[, 4:10])
  }
  colnames(finalcall) <- c("sample_name", "chr", "gene" ,"cnv", "st_bp", "ed_bp", 
                           "length_kb", "st_exon", "ed_exon", "raw_cov", 
                           "norm_cov", "copy_no", "lratio", "mBIC")
  rownames(finalcall) <- rep("", nrow(finalcall))
  finalcall
}