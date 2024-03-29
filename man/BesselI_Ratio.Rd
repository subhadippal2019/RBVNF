% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MPGMeanSd.R
\name{BesselI_Ratio}
\alias{BesselI_Ratio}
\title{Generates samples from Modfified Ploya Gamma (MPG) distribution.}
\usage{
BesselI_Ratio(x, nu1, nu2 = NULL)
}
\arguments{
\item{x:}{A positive real number.}

\item{nu1:}{A positive real number, order of the BesselI in the denominator.}

\item{nu2:}{A positive real number, order of the BesselI in the numerator.}
}
\value{
BesselI_Ratio(x, nu1, nu2):= BesselI(x, nu1)/BesselI(x, nu2)
}
\description{
Generates samples from Modfified Ploya Gamma (MPG) distribution.
}
\examples{
library(Bessel)
BesselI_Ratio(x=2, nu1=2, nu2=1)
}
