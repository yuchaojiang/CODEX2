\name{mapp_ref}
\alias{mapp_ref}
\docType{data}
\title{
Position reference for pre-computed mappability results.
}
\description{
List consisting of IRanges objects specifying exonic positions whose 
mappabilities are pre-computed across the genome.
}
\usage{data(mapp_ref)}
\value{
List of length 24 with genomic positions of pre-computed mappability of each
chromosome.
}
\details{
Genomic positions for pre-computed mappabilities. Method used is detailed in
\code{\link{getmapp}}.
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
\code{\link{getmapp}}
}
\examples{
# mappability exon reference of chromosome 1
mapp_ref[[1]]
}
\keyword{datasets}