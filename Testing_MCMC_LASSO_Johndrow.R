
n=1000  # NUmber of the samples
p=10    # NUmber of the regression covariates
d=2  # Number of direcions in the direcional data

Num_of_nonzero_beta= round(p*.4)
Min_Non_Zero_beta= 10
Max_Non_Zero_beta=20
#### bbeta is a matrix of dimension p\times d
#bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
sigma_square=1
tau_square=1000

data_lst<-Data_generator_vnf_reg_sparse(n=n,
                                        p=p,
                                        d=d,
                                        SetUp = 3,
                                        NumOfNonZeroBeta=c(Num_of_nonzero_beta,
                                                           Min_Non_Zero_beta,
                                                           Max_Non_Zero_beta
                                        )
)

#lst_BLASSO_Beta_MCMC=MCMC_BLASSO_Dir_regression_sampler_V1(Y=data_lst$Y,
# X=data_lst$X,
#  prior=NULL,
#   beta_init = NULL,
#   MCSamplerSize =1000,
#   lasso_lambda =  0.01,
#   Sample_lasso_lambda = NULL)

lst_BLASSO_Beta_MCMC  =               MCMC_BLASSO_Dir_regression_sampler_V1( Y=data_lst$Y,
                                                                             X=data_lst$X,
                                                                             prior=NULL,
                                                                             beta_init = NULL,
                                                                             MCSamplerSize =1000,
                                                                             n_less_than_p_sampler ="JOHNDROW",
                                                                             lasso_lambda_spec = list(
                                                                               Type="SAMPLE",
                                                                               lasso_lambda=0.01,
                                                                               hyper_lambda_selector= NULL
                                                                             )
)
#lst1=lst
lst=lst_BLASSO_Beta_MCMC
Beta_est=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
data_lst$beta
cbind(c(t(Beta_est)), c(t(data_lst$beta)))
#Plot_MCMC_Diag_Triplet(lst$MC$lasso_lambda_all,y_lab_text = bquote(lambda))


#lst=MCMC_Dir_regression_sampler_V1(Y=data_lst$Y, X=data_lst$X, prior=NULL, beta_init = NULL, MCSamplerSize =5500)

#MC_SIM_Reg_Dir_Data<-list(MC_lst=lst_BLASSO_Beta_MCMC, Sim_Data_lst=data_lst)


