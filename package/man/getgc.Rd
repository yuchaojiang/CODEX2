\name{getgc}
\alias{getgc}
\title{
Get GC content for each exonic target
}
\description{
    Computes GC content for each exon. Will be later used in QC procedure and 
    normalization.
}
\usage{
getgc(chr, ref)
}
\arguments{
    \item{chr}{
    Chromosome returned from \code{\link{getbambed}}}
    \item{ref}{
    IRanges object returned from \code{\link{getbambed}}}
}
\value{
Vector of GC content for each exon.
}
\references{
Team TBD. BSgenome.Hsapiens.UCSC.hg19: Full genome sequences for Homo sapiens 
(UCSC version hg19). R package version 1.3.99.
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{getbambed}},
    \code{\link{qc}},
    \code{\link{normalize}}
}
\examples{
ref <- IRanges(st = 51207851, end = 51207982)
gc <- getgc(chr = 22, ref)
}
\keyword{ package }