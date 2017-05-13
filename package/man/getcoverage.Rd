\name{getcoverage}
\alias{getcoverage}
\title{
Get depth of coverage from whole exome sequencing
}
\description{
Gets depth of coverage for each exon across all samples from whole exome 
sequencing files.
}
\usage{
getcoverage(bambedObj, mapqthres)
}
\arguments{
    \item{bambedObj}{
    Object returned from \code{\link{getbambed}}
    }
    \item{mapqthres}{
    Mapping quality threshold hold of reads.
    }
}
\value{
    \item{Y}{Read depth matrix}
    \item{readlength}{Vector of read length for each sample}
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
    \code{\link{getbambed}}
}
\examples{
library(WES.1KG.WUGSC)
dirPath <- system.file("extdata", package = "WES.1KG.WUGSC")
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
sampnameFile <- file.path(dirPath, "sampname")
sampname <- as.matrix(read.table(sampnameFile))
chr <- 22
bambedObj <- getbambed(bamdir = bamdir, bedFile = file.path(dirPath, 
    "chr22_400_to_500.bed"), sampname = sampname,
    projectname = "CODEX_demo", chr)
bamdir <- bambedObj$bamdir
sampname <- bambedObj$sampname
ref <- bambedObj$ref
projectname <- bambedObj$projectname
chr <- bambedObj$chr
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
readlength <- coverageObj$readlength
}
\keyword{ package }