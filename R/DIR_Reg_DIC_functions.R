


help(besselI)

log_BesselI<-function(x, nu){
  val=log(besselI(x, nu=nu, expon.scaled = TRUE)) + x
  return(val)
}



log_Likelihood_lasso<-function(Y, X, beta,lasso_lambda){
  nu= dim(Y)[2]/2-1
  #browser()
  #beta=beta_EM
    XBeta<- X%*%beta
    XBeta_Norm= apply(X = XBeta,MARGIN = 1, FUN = function(xx){return(norm(as.matrix(xx)))} )
   #TR_Yt_XBeta=  sum(diag(t(Y)%*%XBeta)) # Sum of the exponent term
  log_Likelihood= sum(diag(t(Y)%*%XBeta))+
                  sum(nu*log(XBeta_Norm)- log_BesselI(x=XBeta_Norm, nu=.5))-
                   lasso_lambda*sum(abs(beta))
  return(log_Likelihood)
  }


Calculate_DIC_i<-function(lst){
  Mc_len<- lst$RunDetails$MCSamplerSize
  val=0*(1:Mc_len)-1
  for(iter  in 1:Mc_len){
    val[iter]= -2*log_Likelihood_lasso(Y=lst$RunDetails$Data$Y,
                                 X=lst$RunDetails$Data$X,
                                 beta= lst$MC$Mc_Beta[iter, , ],
                                 lasso_lambda=lst$RunDetails$lasso_lambda)
  }
return(val)
}


#' @export
Calculate_DIC<-function(lst){
  val= Calculate_DIC_i(lst)
  beta_post_avg= apply(lst$MC$Mc_Beta, MARGIN = c(2,3), mean) # theta_bar
  D_theta_bar= -2*log_Likelihood_lasso(Y=lst$RunDetails$Data$Y,
                                       X=lst$RunDetails$Data$X,
                                       beta= beta_post_avg,
                                       lasso_lambda=lst$RunDetails$lasso_lambda)


  pd2=var(val)/2
  pd1= mean(val)- D_theta_bar

    DIC_12=mean(val)+ pd2
    DIC_11=mean(val)+ pd1


    DIC_21= D_theta_bar+ 2*pd1
  DIC_22= D_theta_bar+ 2*pd2
  return(cbind(DIC_11=DIC_11, DIC_12=DIC_12, DIC_21=DIC_21, DIC_22=DIC_22))
}
