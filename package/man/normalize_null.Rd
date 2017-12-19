\name{normalize_null}
\alias{normalize_null}
\title{
Normalization of read depth from whole exome sequencing
}
\description{
Fits a Poisson log-linear model that normalizes the read depth data from whole 
exome sequencing. Includes terms that specifically remove biases due to GC 
content, exon capture and amplification efficiency, and latent systemic 
artifacts.
}
\usage{
normalize_null(Y_qc, gc_qc, K)
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
normObj <- normalize_null(Y_qc_codex2, gc_qc_codex2, K = 1)
}
\keyword{ package }