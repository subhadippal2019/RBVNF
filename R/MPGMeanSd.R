




#' Generates samples from Modfified Ploya Gamma (MPG) distribution.
#'
#' @param x: A positive real number.
#' @param nu1: A positive real number, order of the BesselI in the denominator.
#' @param nu2: A positive real number, order of the BesselI in the numerator.
#' @return BesselI_Ratio(x, nu1, nu2):= BesselI(x, nu1)/BesselI(x, nu2)
#' @examples
#' library(Bessel)
#' BesselI_Ratio(x=2, nu1=2, nu2=1)
#' @export
BesselI_Ratio<-function(x, nu1, nu2=NULL){
        if(is.null(nu2)){ nu2=nu1+1}
        return(exp(log(BesselI(z=x, nu=nu2,expon.scaled = TRUE )) - log(BesselI(z=x, nu=nu1,expon.scaled = TRUE ))))
}






#' Generates samples from Modfified Ploya Gamma (MPG) distribution.
#'
#' @param x: A positive real number.
#' @param nu1: A positive real number, order of the BesselI in the denominator.
#' @param nu2: A positive real number, order of the BesselI in the numerator.
#' @return BesselK_Ratio(x, nu1, nu2):= BesselI(x, nu1)/BesselI(x, nu2)
#' @examples
#' library(Bessel)
#' BesselK_Ratio(x=2, nu1=2, nu2=1)
#' @export
BesselK_Ratio<-function(x, nu1, nu2=NULL){
  if(is.null(nu2)){ nu2=nu1+1}
  return(exp(log(BesselK(z=x, nu=nu2,expon.scaled = TRUE )) - log(BesselK(z=x, nu=nu1,expon.scaled = TRUE ))))
}



#ggplot()+stat_function(fun =BesselK_Ratio, args = list(nu1=nu_1, nu2=nu_2 ))+xlim(.00001, 10)





#' Generates samples from Modfified Ploya Gamma (MPG) distribution.
#'
#' @param alpha: A positive real number. It is a parameter of the distribution.
#' @param nu: A positive real number, order of the MPG distribution.
#' @return Mean (mean) and Standard Deviation (sd) of the MPG(alpha, nu) distribution.
#' @examples
#' library(Bessel)
#' MPG_mean_sd(alpha=5, nu=0)
#' data=rMPG(n=1000, alpha=5, nu=0)
#' #' data=rMPG(n=1000, alpha=50, nu=0)
#' mean(data)
#' sd(data)
#' @export
MPG_mean_sd<-function(alpha=1,nu){

              RatioI<-function(x, nu1, nu2=NULL){
                if(is.null(nu2)){ nu2=nu1+1}
                return(exp(log(BesselI(z=x, nu=nu2,expon.scaled = TRUE )) - log(BesselI(z=x, nu=nu1,expon.scaled = TRUE ))))
              }

              #nu=.5
              #alpha=.01
              mean<-RatioI(alpha, nu1 = nu)/(2*alpha)

              #sd=   sqrt( (R(alpha, nu1 = nu) *(R(alpha, nu1 = nu) -R(alpha, nu1 = nu+1) )   )) /(2*alpha)

              sd=sqrt((  (RatioI(alpha, nu1 = nu))^2 -   RatioI(alpha, nu1 = nu, nu2 = nu+2) ))/(2*alpha)

              return(list(mean=mean,sd=sd ))

}


#' computes the CDF frunction for the Modfified Ploya Gamma (MPG) distribution.
#'
#' @param x:  A positive real number or a vector of real numbers. The value at which the density to be calculated.
#' @param alpha: A positive real number.
#' @param nu: A positive real number, order of the MPG distribution.
#' @param Num_of_Terms: Integer greater than 50.
#'  Number of Terms added in the series represantion of the density.
#' @return Density of the MPG(alpha, nu) evaluated at x
#' @examples
#' library(Bessel)
#' pMPG(.05, alpha=5, nu=0)
#' pMPG(c(.01, .05,.1) , alpha=5, nu=0)
#' plot(function(x){pMPG(x, alpha=5, nu=0, iflog = FALSE)}, xlim=c(.005, .2))
#' @export
pMPG<-function(x,alpha,nu,iflog=FALSE, Num_of_Terms=200){

  K=Num_of_Terms
  x=x*(x>.0001)+.0001*(x<=.0001)
  j_nu_0=bessel_zero_Jnu(nu,1:K);
  J_nuPlus1=besselJ(j_nu_0,(nu+1));
  pMPG_single<-function(t1, a, nu, iflog=iflog){
    return(1-log_survival_nu(t=t1,a=a,nu=nu,j_nu_0=j_nu_0,J_nuPlus1=J_nuPlus1,iflog=iflog))
  }
  val=apply(matrix(x, ncol=1), MARGIN = 1, FUN = function(xvec){pMPG_single(xvec,a=alpha, nu=nu, iflog=iflog )})
  # Need to incorporate the gig approximaton
  return(val)
}











################################

GIG_mean_sd<-function(alpha=1,nu, alt=TRUE){

  RatioK<-function(x, nu1, nu2=NULL){
    if(is.null(nu2)){ nu2=nu1+1}
    return(exp(log(BesselK(z=x, nu=nu2,expon.scaled = TRUE )) - log(BesselK(z=x, nu=nu1,expon.scaled = TRUE ))))
  }

  #nu=.5
  #alpha=.01

  nu1=-nu-1; chi= .5; psi=2*alpha^2
  if( alt){nu1= -nu-1.5}
  mean<-RatioK(alpha, nu1 = nu1)/(2*alpha)

  #sd=   sqrt( (R(alpha, nu1 = nu) *(R(alpha, nu1 = nu) -R(alpha, nu1 = nu+1) )   )) /(2*alpha)

  sd=sqrt((  RatioK(alpha, nu1 = nu1, nu2 = nu1+2)-(RatioK(alpha, nu1 = nu1))^2  ))/(2*alpha)

  return(list(mean=mean,sd=sd ))

}

