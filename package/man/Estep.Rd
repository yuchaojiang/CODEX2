\name{Estep}
\alias{Estep}
\title{
Expectation step for CODEX2 normalization with negative control regions
}
\description{
Expectation step for CODEX2 normalization with negative control regions
}
\usage{
Estep(Ytemp, h, beta0, g, mu, pi, Ntemp, fGC)
}
\arguments{
    \item{Ytemp}{Read depth matrix}
    \item{h}{Sample specific latent factors}
    \item{beta0}{Exon-specific bias from the M-step}
    \item{g}{Exon-specific latent factors from the M-step}
    \item{mu}{Deviation on the log scale from the M-step}
    \item{pi}{CNV incidence rate from the M-step}
    \item{Ntemp}{Total number of reads}
    \item{fGC}{GC content bias matrix}
}
\value{
    \item{Zhat }{Matrix of optimized missing data to be fed into M-step}
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\keyword{ package }