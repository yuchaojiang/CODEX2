\name{choiceofK}
\alias{choiceofK}
\title{
Determine the number of latent factors K.
}
\description{
Determines the number of latent variables K via AIC, BIC, and deviance 
reduction. A pdf file containing all three plots is generated.
}
\usage{
choiceofK(AIC, BIC, RSS, K, filename)
}
\arguments{
    \item{AIC}{ vector of AIC for each K returned from \code{\link{normalize}}}
    \item{BIC}{ vector of BIC for each K returned from \code{\link{normalize}}}
    \item{RSS}{ vector of RSS for each K returned from \code{\link{normalize}}}
    \item{K}{ vector of K returned from \code{\link{normalize}}}
    \item{filename}{ Filename of the output plot of AIC and RSS}
}
\value{
pdf file with three plots: AIC, BIC, and percentage of variance explained versus
the number of latent factors.
}
\details{
AIC: Akaike information criterion, used for model selection; BIC: Bayesian 
information criterion, used for model selection; RSS: residue sum of squares, 
used to plot scree plot and determine the 'elbow'.
}
\author{
Yuchao Jiang \email{yuchaoj@wharton.upenn.edu}
}
\seealso{
\code{\link{normalize}},
\code{\link{segment}}
}
\examples{
AIC <- normObjDemo$AIC
BIC <- normObjDemo$BIC
RSS <- normObjDemo$RSS
K <- normObjDemo$K
projectname <- bambedObjDemo$projectname
chr <- bambedObjDemo$chr
choiceofK(AIC, BIC, RSS, K, filename = paste(projectname, "_", chr, 
    "_choiceofK", ".pdf", sep = ""))}
\keyword{ package }