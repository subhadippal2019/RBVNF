











n=500  # NUmber of the samples
p=20    # NUmber of the regression covariates
d=2   # Number of direcions in the direcional data

Num_of_nonzero_beta= round(p*.10)
Min_Non_Zero_beta= 1
Max_Non_Zero_beta=10
#### bbeta is a matrix of dimension p\times d
#bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
sigma_square=1
tau_square=1000


#######################################################################################################
#######################################################################################################
#######################################################################################################

for(SimNum in 56:85){


  #data_lst = Data_generator_vnf_reg(n=n, p=p, d=d, concentration_factor = 1, beta_factor = 10)

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
                                                                                 MCSamplerSize =5500,
                                                                                 lasso_lambda_spec = list(
                                                                                                            Type="SAMPLE",
                                                                                                            lasso_lambda=0.01,
                                                                                                            hyper_lambda_selector= NULL
                                                                                                          )
                                                                                                    )
    #lst=lst_BLASSO_Beta_MCMC
  #Beta_est=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
  #Plot_MCMC_Diag_Triplet(lst$MC$lasso_lambda_all,y_lab_text = bquote(lambda))


  #lst=MCMC_Dir_regression_sampler_V1(Y=data_lst$Y, X=data_lst$X, prior=NULL, beta_init = NULL, MCSamplerSize =5500)

  MC_SIM_Reg_Dir_Data<-list(MC_lst=lst_BLASSO_Beta_MCMC, Sim_Data_lst=data_lst)


  save(MC_SIM_Reg_Dir_Data, file=paste0("/Users/subhadippal/Desktop/Lasso_Simulation_RBVNF/MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",d, "_SimNUmber_",SimNum, ".RData" ))
  #save(MC_SIM_Reg_Dir_Data, file=paste0("G:\\My Drive\\MyRPackages\\MC_SIM_Reg_Dir_Data\\MC_SIM_Reg_Dir_Data_n_750_d_eq_",d, "_SimNUmber_",SimNum, ".RData" ))

  print(paste0("sim Number", SimNum, "completed"))
  rm(list=c("MC_SIM_Reg_Dir_Data","lst_BLASSO_Beta_MCMC" , "data_lst" ))
}



