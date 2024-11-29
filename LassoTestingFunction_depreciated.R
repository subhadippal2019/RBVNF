# Gaussian
x = matrix(rnorm(100 * 20), 100, 20)
y = rnorm(100)
fit1 = glmnet(x = x, y = y, family = "gaussian", lambda=1)
print(fit1)
coef(fit1, s = 0.01)  # extract coefficients at a single value of lambda
predict(fit1, newx = x[1:10, ], s = c(0.01, 0.005))  # make predictions


#Test Lasso
data_lst = Data_generator_vnf_reg(n=200, p=p, d=d, concentration_factor = 1, beta_factor = 5)

data_lst = Data_generator_vnf_reg_sparse(n = 200, p = 20, d = 2, NumOfNonZeroBeta = c(6, 5, 10), SetUp = 2)
Y = data_lst$Y;X=data_lst$X;

beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)


beta_EM_Lasso=EM_BLASSO_Dir_regression_optimizer_V1(Y=Y, X=X, beta_init = NULL,   EM_tolerence = .00001, lasso_lambda = .002)


beta_EM_Lasso_LAM=EM_BLASSO_Dir_regression_optimizer_V1_EMLAMBDA(Y=Y, X=X, beta_init = NULL,   EM_tolerence = .00001, lasso_lambda = .01)



lst_BLASSO_Beta_MCMC=MCMC_BLASSO_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =1000, lasso_lambda = .5, Sample_lasso_lambda = c(10, 100))
lst=lst_BLASSO_Beta_MCMC
Beta_est=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
abs(Beta_est1)-abs(Beta_est)




lst_Beta_MCMC=MCMC_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =1000)
lst1=lst_Beta_MCMC
Beta_est1=apply(lst1$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
Beta_est1-Beta_est

i=1;j= 1
Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))
Beta_est=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
Beta_sd=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = sd)
Beta_est1<- matrix(paste0(  round(c(Beta_est),2),"(", round(c(Beta_sd),2),")& "), nrow=10)
Beta_est2= cbind(paste(colnames(X),"&"), Beta_est1, paste("\\\\"))
paste0(Beta_est2, collapse="//")
cat(Beta_est2)





Y_vec=c(Y)
X_kron=kronecker(diag(d), X)

Standard_Lasso_fit_vanila = glmnet(x =X_kron , y = Y_vec, family = "gaussian",  alpha = 1)

