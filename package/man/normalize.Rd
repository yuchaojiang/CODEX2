\name{normalize}
\alias{normalize}
\title{
Normalization of read depth
}
\description{
Fits a Poisson log-linear model that normalizes the read depth data from whole 
exome sequencing. Includes terms that specifically remove biases due to GC 
content, exon capture and amplification efficiency, and latent systemic 
artifacts.
}
\usage{
normalize(Y_qc, gc_qc, K, N, message)
}
\arguments{
    \item{Y_qc}{
    Read depth matrix after quality control procedure returned from 
    \code{\link{qc}}}
    \item{gc_qc}{
    Vector of GC content for each exon after quality control procedure returned
    from \code{\link{qc}}}
    \item{K}{
    Number of latent Poisson factors. Can be an integer if optimal solution has
    been chosen or a vector of integers so that AIC, BIC, and RSS are computed
    for choice of optimal k.}
    \item{N}{
    Library size factor, which is computed from the genome-wide read depth data.}
    \item{message}{
    Whether display message from normalization. Default to be TRUE.}
}
\value{
    \item{Yhat }{Normalized read depth matrix}
    \item{fGC.hat }{Estimated GC content bias matrix}
    \item{beta.hat }{Estimated exon-specific bias as a vector}
    \item{g.hat }{Estimated Poisson latent factor}
    \item{h.hat }{Estimated Poisson latent factor}
    \item{AIC }{AIC for model selection}
    \item{BIC }{BIC for model selection}
    \item{RSS }{RSS for model selection}
    \item{K }{Number of latent Poisson factors}
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\seealso{
    \code{\link{qc}},
    \code{\link{choiceofK}}
}
\examples{
Y_qc <- Y_qc_demo
gc_qc <- ref_qc_demo$gc
Y.nonzero=Y_qc[apply(Y_qc,1,function(x){!any(x==0)}),]
pseudo = apply(Y.nonzero,1,function(x){prod(x)^(1/length(x))})
N=apply(apply(Y.nonzero,2,function(x){x/pseudo}),2,median)
normObj <- normalize(Y_qc, gc_qc, K = 1:2, N = N)
Yhat <- normObj$Yhat
AIC <- normObj$AIC; BIC <- normObj$BIC; RSS <- normObj$RSS
K <- normObj$K
}
\keyword{ package }