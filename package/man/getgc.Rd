\name{getgc}
\alias{getgc}
\title{
Get GC content
}
\description{
    Computes GC content for each exon/target/window. Will be later used in QC procedure and normalization.
}
\usage{
getgc(ref)
}
\arguments{
    \item{ref}{
    GRanges object returned from \code{\link{getbambed}}}
    \item{genome}{Genome assembly, with a default of BSgenome.Hsapiens.UCSC.hg19}
}
\value{
Vector of GC content for each exon.
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\seealso{
    \code{\link{getbambed}},
    \code{\link{normalize}}
}
\examples{
ref <- GRanges(seqname = 22, ranges = IRanges(st = 51207851, end = 51207982))
gc <- getgc(ref, genome = BSgenome.Hsapiens.UCSC.hg19)
}
\keyword{ package }