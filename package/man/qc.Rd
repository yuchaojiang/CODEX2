\name{qc}
\alias{qc}
\title{
Quality control procedure for depth of coverage
}
\description{
Applies a quality control procedure to the depth of coverage matrix both 
sample-wise and exon-wise before normalization.
}
\usage{
qc(Y, sampname, chr, ref, mapp, gc,cov_thresh,length_thresh,mapp_thresh,
  gc_thresh)
}
\arguments{
    \item{Y}{
    Original read depth matrix returned from \code{\link{getcoverage}}}
    \item{sampname}{
    Vector of sample names returned from \code{\link{getbambed}}}
    \item{chr}{
    Chromosome.}
    \item{ref}{
    IRanges object specifying exonic positions returned from 
    \code{\link{getbambed}}}
    \item{mapp}{
    Vector of mappability for each exon returned from \code{\link{getmapp}}}
    \item{gc}{
    Vector of GC content for each exon returned from \code{\link{getgc}}}
    \item{cov_thresh}{
    Vector specifying the upper and lower bound of exonic median coverage 
    threshold for QC. 20-4000 recommended.}
    \item{length_thresh}{
    Vector specifying the upper and lower bound of exonic length threshold for
    QC. 20-2000 recommended.}
    \item{mapp_thresh}{
    Scalar variable specifying exonic mappability threshold for QC. 0.9 
    recommended.}
    \item{gc_thresh}{
    Vector specifying the upper and lower bound of exonic GC content threshold 
    for QC. 20-80 recommended.}
}
\details{
It is suggested that analysis by CODEX be carried out in a batch-wise fashion if
multiple batches exist. CODEX further filters out exons that: have extremely low
coverage--median read depth across all samples less than 20 or greater than 
4000; are extremely short--less than 20 bp; are extremely hard to map--
mappability less than 0.9; have extreme GC content--less than 20 or greater than
80. The above filtering thresholds are recommended and can be user-defined to be
adapted to different sequencing protocols.}
\value{
    \item{Y_qc}{Updated \code{Y} after QC}
    \item{sampname_qc}{Updated \code{sampname} after QC}
    \item{gc_qc}{Updated \code{gc} after QC}
    \item{mapp_qc}{Updated \code{mapp} after QC}
    \item{ref_qc}{Updated \code{ref} after QC}
    \item{qcmat}{Matrix specifying results of exon-wise QC procedures}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{getbambed}},
    \code{\link{getgc}},
    \code{\link{getmapp}}
}
\examples{
Y <- coverageObjDemo$Y
sampname <- bambedObjDemo$sampname
chr <- bambedObjDemo$chr
ref <- bambedObjDemo$ref
gc <- gcDemo
mapp <- mappDemo
cov_thresh <- c(20, 4000)
length_thresh <- c(20, 2000)
mapp_thresh <- 0.9
gc_thresh <- c(20, 80)
qcObj <- qc(Y, sampname, chr, ref, mapp, gc, cov_thresh, length_thresh, 
    mapp_thresh, gc_thresh)
Y_qc <- qcObj$Y_qc
sampname_qc <- qcObj$sampname_qc
gc_qc <- qcObj$gc_qc
mapp_qc <- qcObj$mapp_qc
ref_qc <- qcObj$ref_qc
qcmat <- qcObj$qcmat
}
\keyword{ package }