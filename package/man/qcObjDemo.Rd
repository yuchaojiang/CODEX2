\name{qcObjDemo}
\alias{qcObjDemo}
\docType{data}
\title{
Demo data pre-stored for qcObj.
}
\description{
Pre-stored qcObj data for demonstration purposes.
}
\usage{data(qcObjDemo)}
\value{
qcObj demo data (list) pre-computed.
}
\details{
Pre-computed using whole exome sequencing data of 46 HapMap samples.
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\examples{
Y_qc <- qcObjDemo$Y_qc
sampname_qc <- qcObjDemo$sampname_qc
gc_qc <- qcObjDemo$gc_qc
mapp_qc <- qcObjDemo$mapp_qc
ref_qc <- qcObjDemo$ref_qc
}
\keyword{datasets}