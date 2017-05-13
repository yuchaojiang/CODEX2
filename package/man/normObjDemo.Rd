\name{normObjDemo}
\alias{normObjDemo}
\docType{data}
\title{
Demo data pre-stored for normObj.
}
\description{
Pre-stored normObj data for demonstration purposes.
}
\usage{data(normObjDemo)}
\value{
normObj demo data (list) pre-computed.
}
\details{
Pre-computed using whole exome sequencing data of 46 HapMap samples.
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\examples{
Yhat <- normObjDemo$Yhat
AIC <- normObjDemo$AIC
BIC <- normObjDemo$BIC
RSS <- normObjDemo$RSS
K <- normObjDemo$K
}
\keyword{datasets}