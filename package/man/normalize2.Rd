\name{normalize2}
\alias{normalize2}
\title{
Normalization of read depth from whole exome sequencing under the case-control
setting
}
\description{
Fits a Poisson log-linear model that normalizes the read depth data from whole 
exome sequencing. Includes terms that specifically remove biases due to GC 
content, exon capture and amplification efficiency, and latent systemic 
artifacts. If the WES is designed under case-control setting, CODEX estimates
the exon-wise Poisson latent factor using only the read depths in the control
cohort, and then computes the sample-wise latent factor terms for the case 
samples by regression.
}
\usage{
normalize2(Y_qc, gc_qc, K, normal_index)
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
    \item{normal_index}{
    Indices of control samples.}
}
\value{
    \item{Yhat }{Normalized read depth matrix}
    \item{AIC }{AIC for model selection}
    \item{BIC }{BIC for model selection}
    \item{RSS }{RSS for model selection}
    \item{K }{Number of latent Poisson factors}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{qc}},
    \code{\link{choiceofK}}
}
\examples{
Y_qc <- qcObjDemo$Y_qc
gc_qc <- qcObjDemo$gc_qc
normObj <- normalize2(Y_qc, gc_qc, K = 1:5, normal_index = seq(1, 45, 2))
Yhat <- normObj$Yhat
AIC <- normObj$AIC
BIC <- normObj$BIC
RSS <- normObj$RSS
K <- normObj$K
}
\keyword{ package }