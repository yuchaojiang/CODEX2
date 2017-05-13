\name{getmapp}
\alias{getmapp}
\title{
Get mappability for each exonic target
}
\description{
Computes mappability for each exon. To save running time, take values from 
pre-computed results. Will be later used in QC procedure.
}
\usage{
getmapp(chr, ref)
}
\arguments{
    \item{chr}{
    Chromosome returned from \code{\link{getbambed}}}
    \item{ref}{
    IRanges object returned from \code{\link{getbambed}}}
}
\details{
To calculate the exonic mappability, we first construct consecutive reads of 
length 90 that are one base pair apart along the exon. The sequences are taken 
from the hg19 reference. We then find possible positions across the genome that
the reads can map to allowing for two mismatches. We compute the mean of the 
probabilities that the overlapped reads map to the target places where they are
generated and use this as the mappability of the exon.
}
\value{
Vector of mappability for each exon.
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
\code{\link{getbambed}},
\code{\link{qc}}
}
\examples{
ref <- IRanges(st = 51207851, end = 51207982)
mapp <- getmapp(chr = 22, ref)
}
\keyword{ package }