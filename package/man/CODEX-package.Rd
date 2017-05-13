\name{CODEX-package}
\alias{CODEX-package}
\alias{CODEX}
\docType{package}
\title{
A Normalization and Copy Number Variation Detection Method for Whole Exome 
Sequencing
}
\description{
CODEX is a normalization and copy number variation calling procedure for whole 
exome DNA sequencing data. CODEX relies on the availability of multiple samples 
processed using the same sequencing pipeline for normalization, and does not 
require matched controls. The normalization model in CODEX includes terms that 
specifically remove biases due to GC content, exon length and targeting and 
amplification efficiency, and latent systemic artifacts. CODEX also includes a 
Poisson likelihood-based recursive segmentation procedure that explicitly models
the count-based exome sequencing data.
}
\details{
\tabular{ll}{
Package: \tab CODEX\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2015-01-13\cr
License: \tab GPL-2\cr
}
CODEX takes as input the bam files/directories for whole exome sequencing 
datasets and bed files for exonic positions, returns raw and normalized coverage
for each exon, and calls copy number variations with genotyping results.
}
\author{
Yuchao Jiang <yuchaoj@wharton.upenn.edu>, Nancy R. Zhang
}
\keyword{ package }