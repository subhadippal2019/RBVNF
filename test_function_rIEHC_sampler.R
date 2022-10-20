
epsilon=.001
I1=integrate(function(x){exp(-epsilon*x)/(1+x)}, lower = 0, upper = Inf)
I2=integrate(function(x){x*exp(-epsilon*x)/(1+x)}, lower = 0, upper = Inf)

I2$value/I1$value


sam<-replicate(n = 100000,expr = sample_h_u_multiple(epsilon = epsilon, a = .8, b = 2))
mean(sam)
sam1<-replicate(n = 100000, expr = sample_h_u_multiple_alt(epsilon = epsilon, del = .001))
mean(sam1)

sam2<-replicate(n = 100000, expr = (sample_h_u_multiple_alt1(n = 1, epsilon = epsilon,ifLog = FALSE )))
mean(sam2)






Inv_IEHC_density_unscaled<-function(x,A, p, iflog=TRUE){
  #if(any(x)<0){return(0)}
  val= ((p-1)/2)* log(x)-x-log(A+x)
  #x^((p-1)/2)*exp(-x)/((A+x))
  if( !iflog){  val=exp(val)  }
  return(val)
}


Inv_IEHC_density_unscaled_mean<-function(x,A, p, iflog=TRUE){
  #if(any(x)<0){return(0)}
  val= ((p-1)/2)* log(x)-x-log(A+x)+log(x)
  #x^((p-1)/2)*exp(-x)/((A+x))
  if( !iflog){  val=exp(val)  }
  return(val)
}

A=1000;p=2
I1=integrate(function(x){Inv_IEHC_density_unscaled_mean(x, A=A, p=p, iflog = FALSE)}, lower = 0, upper = Inf)
I2=integrate(function(x){Inv_IEHC_density_unscaled(x, A=A, p=p, iflog = FALSE)}, lower = 0, upper = Inf)

I1$value/I2$value




f_proposal<-function(x, A,m=NULL){

  #if(is.character(m)){m=(sqrt(gamma^2+8*(alpha-1)*beta)-gamma)/(4*beta)}
  if(is.null(m)){m=1}

  val=m^(m/(m+A)) * x^(-m/(m+A))* exp(-x)/(A+m)
  return(val)

}


f_unscaled<-function(x,A){
  if(sum(x)<0){return(0)}
  val=exp(-x)/((A+x))
  return(val)
}


eps=.00000000001

m_opt<-function(A){
  g<-function(x){

    return(log(x)-digamma(A/(A+x)))
  }
  m_opt=uniroot(f = g, interval = c(eps, 10+2/A), tol = .00000000001 )
  return(m_opt$root)
}


f<-function(m){
  Am=A*m
  frac=Am/(1+Am)
  val=frac*log(m)+lgamma(frac)-log(Am+1)
  return(val)
}

A=.01
m_op<-m_opt(A )

I1=integrate(function(x){f_unscaled(x, A)}, lower = 0, upper = Inf)
I2=integrate(function(x){f_proposal(x, A, m_op)}, lower = 0, upper = Inf)

I1$value/I2$value









f_proposal<-function(x, A,p=1, m=NULL){

  #if(is.character(m)){m=(sqrt(gamma^2+8*(alpha-1)*beta)-gamma)/(4*beta)}
  if(is.null(m)){m=1}

  val=m^(m/(m+A)) * x^((p-1)/2+A/(m+A)-1)* exp(-x)/(A+m)
  return(val)

}


f_unscaled<-function(x,A, p=1){
  if(sum(x)<0){return(0)}
  val=x^((p-1)/2)*exp(-x)/((A+x))
  return(val)
}


eps=.00001

m_opt<-function(A,p=1){
  g<-function(x){

    return(log(x)-digamma((p-1)/2+A/(A+x)))
  }
  m_opt=uniroot(f = g, interval = c(eps, 10+2/A), tol = .00000000001 )
  return(m_opt$root)
}



A=1000
p=2
m_op<-m_opt(A , p)

I1=integrate(function(x){f_unscaled(x, A, p=p)}, lower = 0, upper = Inf)
I2=integrate(function(x){f_proposal(x, A,p=p,  m_op)}, lower = 0, upper = Inf)

I1$value/I2$value













IEHC_density_unscaled<-function(x,a, b, iflog=TRUE){
  #if(any(x)<0){return(0)}
  val= -b/x^2-a*log(x)-log(1+x^2)
  #x^((p-1)/2)*exp(-x)/((A+x))
  if( !iflog){  val=exp(val)  }
  return(val)
}


IEHC_density_unscaled_mean<-function(x,a, b, iflog=TRUE){
  #if(any(x)<0){return(0)}
  val= -b/x^2-a*log(x)-log(1+x^2)+log(x)
  #x^((p-1)/2)*exp(-x)/((A+x))
  if( !iflog){  val=exp(val)  }
  return(val)
}

b=1;a=10
I1=integrate(function(x){IEHC_density_unscaled_mean(x, a=a, b=b, iflog = FALSE)}, lower = 0, upper = Inf)
I2=integrate(function(x){IEHC_density_unscaled(x, a=a, b=b, iflog = FALSE)}, lower = 0, upper = Inf)

I1$value/I2$value

x<-rIEHC(n = 100000, a=a,b=b);mean(x)

