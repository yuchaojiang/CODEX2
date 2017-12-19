\name{Mstep}
\alias{Mstep}
\title{
Maximization step for CODEX2 normalization with negative control regions
}
\description{
Maximization step for CODEX2 normalization with negative control regions
}
\usage{
Mstep(Ytemp, h, Zhat, Ntemp, fGC, mu)
}
\arguments{
    \item{Ytemp}{Read depth matrix}
    \item{h}{Sample specific latent factors}
    \item{Zhat}{Matrix of optimized missing data}
    \item{Ntemp}{Total number of reads}
    \item{fGC}{GC content bias matrix}
    \item{mu}{Deviation on the log scale}
}
\value{
    \item{beta0 }{Estimated exon-specific bias for the exon of interest}
    \item{g }{Estimated exon-specific latent factor for the exon of interest}
    \item{pi }{CNV incidence rate}
    \item{mu }{Deviation (log scale) of CNV normalized signal from null}
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\keyword{ package }