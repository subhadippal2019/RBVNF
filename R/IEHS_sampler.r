

#' Generates samples from Inverse Expontially Tilted Half Cauchy (IEHC) distribution.
#'
#' @param n:  Number of samples to be generated.
#' @param a: A positive real number, IEHCshape parameter in latex.
#' @param b: A positive real number. IEHCscale parameter in latex.
#' @return Samples from the IEHC(n, a, b) distribution.
#' @examples
#' data=rIEHC(n=1000, a=1, b=.1)
#' @export
rIEHC<-function(n=1, a, b){
        val=rInvIEHC(n = n, a = a, b = b )
        sample=sqrt(b/val)
        return(sample)
}

#############################################################################
#############################################################################

#' @export
rInvIEHC<-function(n=1, a, b){
      #browser()
      m=IEHC_m_opt(a= a, b = b)
      x_sample_prop<-rgamma(n = n, shape = (a-1)/2+b/(m+b), scale = 1)
      #print(m);print(b);print((a-1)/2+b/(m+b))
      log_accept<-Inv_IEHC_density_unscaled(x =x_sample_prop, a=a, b=b, iflog = TRUE)   - Inv_IEHC_proposal(x = x_sample_prop, a=a, b=b, m = m, iflog = TRUE)
      x_sample=x_sample_prop[log(runif(n))<log_accept]
      size=length(x_sample)

              while(size<n){
                        x_sample_prop<-rgamma(n = n-size, shape = (a-1)/2+b/(m+b), scale = 1)
                        log_accept<-Inv_IEHC_density_unscaled(x =x_sample_prop, a = a, b= b, iflog = TRUE)  - Inv_IEHC_proposal(x = x_sample_prop, a = a, b = b, m = m, iflog = TRUE)
                        x_sample_new=x_sample_prop[log(runif(n-size))<log_accept]
                        x_sample=c(x_sample,x_sample_new )
                        size=length(x_sample)
              }

      return(x_sample)
}

#############################################################################
#############################################################################

Inv_IEHC_proposal<-function(x, a,b,  m=NULL, iflog=TRUE){
      if(is.null(m)){  m=IEHC_m_opt(a = a, b = b) }
        val= (m/(m+b))* log(m)+ ((a-1)/2+b/(m+b)-1)* log(x)-x-log(b+m)#val=m^(m/(m+A)) * x^((p-1)/2+A/(m+A)-1)* exp(-x)/(A+m)
      if( !iflog){  val=exp(val)  }
      return(val)
}

#############################################################################
#############################################################################
Inv_IEHC_density_unscaled<-function(x,a, b, iflog=TRUE){
        val= ((a-1)/2)* log(x)-x-log(b+x)
         #x^((p-1)/2)*exp(-x)/((A+x))
        if( !iflog){  val=exp(val)  }
        return(val)
}



#############################################################################
#############################################################################

# Need to improve this when b is small assume a>=1 for safe results
IEHC_m_opt<-function(a,b=1, lower_search_lim=.0000000000001, conv_tol_level=.00000000001){
      g<-function(x){ return(log(x)-digamma((a-1)/2+b/(b+x))) }
      m_opt=uniroot(f = g, interval = c(lower_search_lim, 2*a+2/b), tol  = conv_tol_level )
      return(m_opt$root)
}

#############################################################################
#############################################################################
