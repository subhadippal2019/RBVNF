# Standard Optimization algorithm instead of EM


n=100  # NUmber of the samples
p=10  # NUmber of the regression covariates
d=2
data_lst = Data_generator_vnf_reg(n=n, p=p, d=d, concentration_factor = 1, beta_factor = 10)
X=data_lst$X
Y=data_lst$Y


# Example structure using hmclearn
library(hmclearn)

# 1. Define Log-Posterior
logPOSTERIOR <- function(theta) {
  return(log_Likelihood_for_HMC(Y = Y, X = X, beta_vec = theta))
}

# 2. Define Gradient of Log-Posterior
glogPOSTERIOR <- function(theta) {
  return(c(gradient_log_Likelihood_for_HMC(Y = Y, X = X, beta_vec = theta)))
}



log_Likelihood_for_HMC<-function(beta_vec, Y, X){
  nu= dim(Y)[2]/2-1
  p=dim(X)[2]
  #browser()
  #beta=beta_EM
  beta=matrix(beta_vec, nrow=p)
  XBeta<- X%*%beta
  XBeta_Norm= apply(X = XBeta,MARGIN = 1, FUN = function(xx){return(norm(as.matrix(xx)))} )
  #TR_Yt_XBeta=  sum(diag(t(Y)%*%XBeta)) # Sum of the exponent term
  log_Likelihood= sum(diag(t(Y)%*%XBeta))+
    sum(nu*log(XBeta_Norm)- log_BesselI(x=XBeta_Norm, nu=.5))
  #-  lasso_lambda*sum(abs(beta))
  return(log_Likelihood)
}

gradient_log_Likelihood_for_HMC<-function(beta_vec, Y, X){
  nu= dim(Y)[2]/2-1
  p=dim(X)[2]
  n=dim(Y)[1]
  d= dim(Y)[2]
  #browser()
  #beta=beta_EM
  beta=matrix(beta_vec, nrow=p)
  XBeta<- X%*%beta
  XBeta_Norm= apply(X = XBeta,MARGIN = 1, FUN = function(xx){return(norm(as.matrix(xx)))} )
  R_nu_x_by_x=BesselI_Ratio(x =XBeta_Norm ,nu1 = nu+1, nu2 = nu)/XBeta_Norm

  #browser()
  SUM_vec_y_i_KRON_x_i= (kronecker(X[1,],Y[1,]))
  Sum_I_kron_x_x_t_betavec=R_nu_x_by_x[1]*kronecker(diag(d), (X[1, ])%*%t(X[1, ]))%*%beta_vec
  for(i in 2:n){
    SUM_vec_y_i_KRON_x_i=SUM_vec_y_i_KRON_x_i+ (kronecker(X[i,],Y[i,]))
    Sum_I_kron_x_x_t_betavec=Sum_I_kron_x_x_t_betavec+R_nu_x_by_x[i]*kronecker(diag(d), (X[i, ])%*%t(X[i, ]))%*%beta_vec
  }
  gradient_log_Likelihood_for_HMC=c(SUM_vec_y_i_KRON_x_i)-c(Sum_I_kron_x_x_t_betavec)

  # Need to Add prior Infor
  return(gradient_log_Likelihood_for_HMC)
}



#install.packages("maxLik")
library(maxLik)

#Initial Values for the loop
# beta_init<- solve(t(X)%*%X + diag(p)/prior$beta_sigma_sq)%*%secondPartOptEq   }
beta_init_vec<- c(solve(t(X)%*%X)%*%t(X)%*%Y )

res <- maxLik(logLik = logPOSTERIOR , grad = glogPOSTERIOR,  start=c(beta_init_vec), method = "CG") # use 'wrong' start values



start_beta= secondPartOptEq<-(t(X)%*%Y+ prior$beta_mean/prior$beta_sigma_sq)





Beta_est=apply(lst_BLASSO_Beta_MCMC$MC$Mc_Beta_n_large_p, MARGIN = c(2,3), FUN = mean)















theta_init<- RBVNF::EM_Dir_regression_optimizer_V1(Y = Y, X = X)
#c(solve(t(X)%*%X)%*%t(X)%*%Y)
#c(data_lst$beta)
# 3. Run HMC
fit <- hmc(
  N = 1000,
  theta.init =c(theta_init),# .5+c(rep(0, p*d)),
  epsilon = 0.008,
  L = 10,
  logPOSTERIOR = logPOSTERIOR,
  glogPOSTERIOR = glogPOSTERIOR,
  param = list(X = X, Y=Y), parallel=FALSE, chains=1
)
diagplots(fit, burnin=1, comparison.theta=c(data_lst$beta))







