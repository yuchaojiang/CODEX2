\name{segmentCBS}
\alias{segmentCBS}
\title{
Recursive segmentation algorithm for CNV detection and genotyping
}
\description{
Recursive segmentation algorithm for CNV detection and genotyping, using 
normalized read depth from whole exome sequencing. Updated from CODEX.
}
\usage{
segmentCBS(Y_qc, Yhat, optK, K, sampname_qc, ref_qc, chr, lmax, mode)
}
\arguments{
    \item{Y_qc}{
    Raw read depth matrix after quality control procedure
    }
    \item{Yhat}{Normalized read depth matrix
    }
    \item{optK}{
    Optimal value \code{K} returned }
    \item{K}{
    Number of latent Poisson factors. Can be an integer if optimal solution has
    been chosen or a vector of integers so that AIC, BIC, and RSS are computed
    for choice of optimal k.}
    \item{sampname_qc}{
    Vector of sample names after quality control procedure}
    \item{ref_qc}{
    IRanges object of genomic positions of each exon after quality control 
    procedure}
    \item{chr}{
    Chromosome number}
    \item{lmax}{
    Maximum CNV length in number of exons returned.}
    \item{mode}{
    Can be either "integer" or "fraction", which respectively correspond to 
    format of the returned copy numbers.}
}
\value{
    Final callset of CNVs with genotyping results.
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\examples{
library(CODEX2)
Y_qc <- qcObjDemo$Y_qc
Yhat <- normObjDemo$Yhat
BIC <- normObjDemo$BIC
K <- normObjDemo$K
sampname_qc <- qcObjDemo$sampname_qc
ref_qc <- qcObjDemo$ref_qc
chr <- bambedObjDemo$chr
finalcall <- segment(Y_qc, Yhat, optK = K[which.max(BIC)], K = K, sampname_qc,
    ref_qc, chr, lmax = 200, mode = "integer")
finalcall
}
\keyword{ package }