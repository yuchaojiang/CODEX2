\name{normalize_codex2_nr}
\alias{normalize_codex2_nr}
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
normalize_codex2_nr(Y_qc, gc_qc, K, cnv_index, N)
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
    \item{cnv_index}{
    Indices of exons/targets within common CNV regions. Can be empirically learnt
    from data}
    \item{N}{
    Library size factor, which is computed from the genome-wide read depth data.}
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
\examples{
Y_qc <- Y_qc_demo
gc_qc <- ref_qc_demo$gc
Y.nonzero=Y_qc[apply(Y_qc,1,function(x){!any(x==0)}),]
pseudo = apply(Y.nonzero,1,function(x){prod(x)^(1/length(x))})
N=apply(apply(Y.nonzero,2,function(x){x/pseudo}),2,median)
normObj <- normalize_codex2_nr(Y_qc, gc_qc, K = 2,
                               cnv_index = 1580:1620, N = N)
}
\keyword{ package }