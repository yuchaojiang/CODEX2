\name{normObj_demo}
\alias{normObj_demo}
\docType{data}
\title{
Normalization results for toy dataset for CODEX2.
}
\description{
Pre-stored normalization results for demonstration.
}
\usage{data(normObj_demo)}
\value{
List of normalized read counts.
}
\details{
Pre-computed from the toy dataset.
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\examples{
head(normObj_demo[[2]])  # normalized read depth with 2 latent factors
}
\keyword{datasets}