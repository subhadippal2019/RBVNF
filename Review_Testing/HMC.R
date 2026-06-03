
n=1000  # NUmber of the samples
p=2   # NUmber of the regression covariates
d=3
data_lst = Data_generator_vnf_reg(n=n, p=p, d=d, concentration_factor = 1, beta_factor = 10)
X=data_lst$X
Y=data_lst$Y


# Example structure using hmclearn
library(hmclearn)

# 1. Define Log-Posterior
logPOSTERIOR <- function(theta, X=X, Y=Y) {
  return(log_Likelihood_for_HMC(Y = Y, X = X, beta_vec = theta))
}

# 2. Define Gradient of Log-Posterior
glogPOSTERIOR <- function(theta, X=X, Y=Y) {
  return(c(gradient_log_Likelihood_for_HMC(Y = Y, X = X, beta_vec = theta)))
}



log_Likelihood_for_HMC<-function(Y, X, beta_vec){
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

gradient_log_Likelihood_for_HMC<-function(Y, X, beta_vec){
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
#
# f<- hmc(N = 1000,
#          theta.init = c(rep(0, 4), 1),
#          epsilon = 0.01,
#          L = 10,
#          logPOSTERIOR = linear_posterior,
#          glogPOSTERIOR = g_linear_posterior,
#          varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#          param=list(y=y, X=X), parallel=FALSE, chains=1)
# diagplots(fit, burnin=100, comparison.theta=c(beta_vec))
#
#


# Y=data_lst$Y
# X=data_lst$X
# beta_vec=c(data_lst$beta)
#
#
# beta_mat=matrix(beta_vec, nrow=p)
#
#
# X_beta=(X%*%beta_mat)
#
# nu=d/2-1
# log_BesselI_X_beta=apply(X_beta, MARGIN = 1, FUN = function(x){(RBVNF::log_BesselI(norm(x), nu=nu))})
# log_norm_X_beta= apply(X_beta, MARGIN = 1, FUN = function(x){(log(norm(x)))})
# #RBVNF::log_BesselI(norm(X_beta), nu=nu)
# sum_xt_beta_y=sum(diag(t(Y)%*%X%*%beta_mat))
# #sum_xt_beta_y=#. ss=0;for(i in 1:100){ss=ss+sum(X_beta[i,]*Y[i,])}
# log_Likelihood_lasso










# Linear regression example
set.seed(522)
X <- cbind(1, matrix(rnorm(300), ncol=3))
betavals <- c(0.5, -1, 2, -3)
y <- X%*%betavals + rnorm(100, sd=.2)
f <- hmc(N = 1000,
         theta.init = c(rep(0, 4), 1),
         epsilon = 0.01,
         L = 10,
         logPOSTERIOR = linear_posterior,
         glogPOSTERIOR = g_linear_posterior,
         varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
         param=list(y=y, X=X), parallel=FALSE, chains=1)
diagplots(f, burnin=300, comparison.theta=c(betavals, 2*log(.2)))

