\name{CODEX2-package}
\alias{CODEX2-package}
\alias{CODEX2}
\docType{package}
\title{
Full-Spectrum Copy Number Variation Detection by High-Throughput DNA Sequencing
}
\description{
High-throughput DNA sequencing enables detection of copy number variations (CNVs) on the genome-wide scale with finer resolution compared to array-based methods, but suffers from biases and artifacts that lead to false discoveries and low sensitivity. We describe CODEX2, a statistical framework for full-spectrum CNV profiling that is sensitive for variants with both common and rare population frequencies and that is applicable to study designs with and without negative control samples. We demonstrate and evaluate CODEX2 on whole-exome and targeted sequencing data, where biases are the most prominent. CODEX2 outperforms existing methods and, in particular, significantly improves sensitivity for common CNVs.
}
\details{
\tabular{ll}{
Package: \tab CODEX2\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2017-12-18\cr
License: \tab GPL-2\cr
}
CODEX2 takes as input the bam files/directories for whole exome sequencing 
datasets and bed files for exonic positions, returns raw and normalized coverage
for each exon, and calls copy number variations with genotyping results.
}
\author{
Yuchao Jiang <yuchaoj@email.unc.edu>, Nancy R. Zhang
}
\keyword{ package }