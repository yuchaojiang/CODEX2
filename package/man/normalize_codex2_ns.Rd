\name{normalize_codex2_ns}
\alias{normalize_codex2_ns}
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
normalize_codex2_ns(Y_qc, gc_qc, K, norm_index)
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
    \item{norm_index}{
    Indices of control samples.}
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
normObj <- normalize_codex2_ns(Y_qc_codex2, gc_qc_codex2, K = 1,
                               norm_index = norm_index_codex2)
}
\keyword{ package }