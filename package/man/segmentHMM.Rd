\name{segmentHMM}
\alias{segmentHMM}
\title{
Hidden Markov model for segmentation of CNV events and genotyping
}
\description{
Hidden Markov model for segmentation of CNV events and genotyping, using 
normalized read depth from whole exome sequencing.  This algorithm incorporates distance between exons as information. 
}
\usage{
segmentHMM(Y_qc, Yhat, optK, K, sampname_qc, ref_qc, chr, mode)
}
\arguments{
    \item{Y_qc}{
    Raw read depth matrix after quality control procedure returned from 
    \code{\link{qc}}}
    \item{Yhat}{Normalized read depth matrix returned from 
    \code{\link{normalize}}}
    \item{optK}{
    Optimal value \code{K} returned from \code{\link{choiceofK}}}
    \item{K}{
    Number of latent Poisson factors. Can be an integer if optimal solution has
    been chosen or a vector of integers so that AIC, BIC, and RSS are computed
    for choice of optimal k.}
    \item{sampname_qc}{
    Vector of sample names after quality control procedure returned from 
    \code{\link{qc}}}
    \item{ref_qc}{
    IRanges object of genomic positions of each exon after quality control 
    procedure returned from \code{\link{qc}}}
    \item{chr}{
    Chromosome number returned from \code{\link{getbambed}}}
    \item{mode}{	
    Can be either "integer" or "fraction", which respectively correspond to format of the     returned copy numbers.  This does not affect the algorithm. }
}
\value{
    Final callset of CNVs with genotyping results.
}
\author{
 Yuchao Jiang \email{yuchaoj@email.unc.edu} 
 
 Gene Urrutia \email{urutia@live.unc.edu} 
}
\examples{
Y_qc <- Y_qc_demo
Yhat <- normObj_demo$Yhat
BIC <- normObj_demo$BIC
K <- normObj_demo$K
sampname_qc <- paste('sample',1:ncol(Y_qc), sep='')
ref_qc <- ref_qc_demo
finalcall <- segmentHMM(Y_qc, Yhat, optK = K[which.max(BIC)], K = K, sampname_qc,
    ranges(ref_qc), chr = 20, mode = "integer")
}
\keyword{ package }

