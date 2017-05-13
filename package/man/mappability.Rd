\name{mappability}
\alias{mappability}
\docType{data}
\title{
Pre-computed mappabilities
}
\description{
The results of pre-computed mappabilities to save running time.
}
\usage{data(mappability)}
\value{
List of length 24 with pre-computed mappability of each chromosome.
}
\details{
Pre-computed mappabilities. Method used is detailed in \code{\link{getmapp}}.
}
\seealso{
\code{\link{getmapp}}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\examples{
# mappability of chromosome 1
head(round(mappability[[1]], 2))
}
\keyword{datasets}