

#' computes the ratio of two BesselI frunctions.
#'
#' @param x:  A positive real number or a vector of real numbers. The value at which the Ratio of Bessel to be calculated.
#' @param nu: A positive real number, the order of BesselI
#' @param ifLog: if ifLog=TRUE  is specified then log of the Bessel function will be retunred
#' @return BesselI evaluated at x>0.
#' \deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{x!}} for \eqn{x = 0, 1, 2, \ldots}.
#' @examples
#' #library(Bessel)
#' BesselI(x=10, nu=5,  ifLog=TRUE)
#' BesselI(x=10000000000, nu=7,  ifLog=TRUE)
#' @export
BesselIBetter<-function(x, nu, ifLog=FALSE ){
  val= besselI(x = x, nu = nu, expon.scaled = TRUE)
  log_val=ifelse(abs(val)>0, log(val),   -.5*(log(2*pi)+log(x)  )) + x  #large sample approximation of I
  #log(besselI( x = x, nu = nu, expon.scaled = TRUE  ) )
  if(!ifLog){
    return(exp(log_val))
  }
  return(log_val)
}






#' computes the ratio of two BesselI frunctions.
#'
#' @param x:  A positive real number or a vector of real numbers. The value at which the Ratio of Bessel to be calculated.
#' @param nuNum: A positive real number, the order of BesselI in the numerator
#' @param nuDeno: A positive real number, the order of BesselI in the denominator
#' @param ifLog: if ifLog=TRUE  is specified then log of the ratio of the Bessel function will be retunred
#'  Number of Terms added in the series represantion of the density.
#' @return Ratio of Bessel function
#' \deqn{   \frac{ I_{nuNum}(x) }{ I_{nuDeno}(x) }  }
#' @examples
#' BesselIRatio(x=10, nuNum=5,nuDeno=7,  ifLog=TRUE)
#' BesselIRatio(x=10000000000, nuNum=5, nuDeno=7,  ifLog=TRUE)
#' @export
BesselIRatio<-function(x, nuNum, nuDeno, ifLog=FALSE){
  #besselI(x= x, nu = nu+1, expon.scaled = TRUE)/besselI(x= x, nu = nu, expon.scaled = TRUE)
  log_val=BesselIBetter(x = x, nu = nuNum, ifLog = TRUE) - BesselIBetter(x = x, nu = nuDeno, ifLog = TRUE)
  if(!ifLog){
    return(exp(log_val))
  }
  return(log_val)
}


#' computes the ratio of two BesselI frunctions.
#'
#' @param x:  A positive real number or a vector of real numbers. The value at which the Ratio of Bessel to be calculated.
#' @param nu: A positive real number, the order of BesselI
#' @param ifLog: if ifLog=TRUE  is specified then log of the Bessel function will be retunred
#' @return BesselI evaluated at x>0.
#' \deqn{ \frac{ I_{nu+1}(x) }{ I_{nu}(x) }  }.
#' @examples
#' BesselIR(x=10, nu=6,  ifLog=TRUE)
#' @export
BesselIR<-function(x, nu, ifLog=FALSE ){
val = BesselIRatio(x = x, nuNum = nu+1, nuDeno = nu, ifLog = ifLog )
return(val)
}


