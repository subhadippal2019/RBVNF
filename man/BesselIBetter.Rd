% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BesselFunctions.r
\name{BesselIBetter}
\alias{BesselIBetter}
\title{computes the ratio of two BesselI frunctions.}
\usage{
BesselIBetter(x, nu, ifLog = FALSE)
}
\arguments{
\item{x:}{A positive real number or a vector of real numbers. The value at which the Ratio of Bessel to be calculated.}

\item{nu:}{A positive real number, the order of BesselI}

\item{ifLog:}{if ifLog=TRUE  is specified then log of the Bessel function will be retunred}
}
\value{
BesselI evaluated at x>0.
\deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{x!}} for \eqn{x = 0, 1, 2, \ldots}.
}
\description{
computes the ratio of two BesselI frunctions.
}
\examples{
#library(Bessel)
BesselI(x=10, nu=5,  ifLog=TRUE)
BesselI(x=10000000000, nu=7,  ifLog=TRUE)
}
