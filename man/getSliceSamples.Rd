% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/SStop.R
\name{getSliceSamples}
\alias{getSliceSamples}
\title{getSliceSamples}
\usage{
getSliceSamples(LL, x0, nBurnin = 10 * 1024, nSample = 10000,
  verbose = FALSE, useFactor = TRUE, useTrueFactors = FALSE)
}
\arguments{
\item{LL}{the loglikelihood function}

\item{x0}{a vector specifying the start for the chain}
}
\value{
a matrix with columns as parameters and rows as post-burnin samples
}
\description{
Function to run automated-factor slice sampling
}
\details{
todo: more about provenance and optional arguments...
}
\examples{
rosen <- function(x) -(1-x[1])^2 - 100*(x[2] - x[1]^2)^2
runss <- getSliceSamples(LL=arosen,x0=c(0,0),useFactor=TRUE)
corplot(runss)
}
\author{
Pete Dodd
}

