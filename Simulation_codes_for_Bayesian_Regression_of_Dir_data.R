

n=1000 # NUmber of the samples
p=10  # NUmber of the regression covariates
d=2 # Number of direcions in the direcional data


#### bbeta is a matrix of dimension p\times d
#bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
sigma_square=1
tau_square=1000


#######################################################################################################
#######################################################################################################
#######################################################################################################

for(SimNum in 1:50){


  data_lst = Data_generator_vnf_reg(n=n, p=p, d=d, concentration_factor = 1, beta_factor = 10)

  lst=MCMC_Dir_regression_sampler_V1(Y=data_lst$Y, X=data_lst$X, prior=NULL, beta_init = NULL, MCSamplerSize =5500)

  MC_SIM_Reg_Dir_Data<-list(MC_lst=lst, Sim_Data_lst=data_lst)


  save(MC_SIM_Reg_Dir_Data, file=paste0("G:\\My Drive\\MyRPackages\\MC_SIM_Reg_Dir_Data\\MC_SIM_Reg_Dir_Data_d_eq_",d, "_SimNUmber_",SimNum, ".RData" ))

  print(paste0("sim Number", SimNum, "completed"))
  rm(list=c("MC_SIM_Reg_Dir_Data","lst" , "data_lst" ))

}




n=1000 # NUmber of the samples
p=10  # NUmber of the regression covariates
d=3 # Number of direcions in the direcional data


#### bbeta is a matrix of dimension p\times d
#bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
sigma_square=1
tau_square=1000


#######################################################################################################
#######################################################################################################
#######################################################################################################

for(SimNum in 1:50){


  data_lst = Data_generator_vnf_reg(n=n, p=p, d=d, concentration_factor = 1, beta_factor = 10)

  lst=MCMC_Dir_regression_sampler_V1(Y=data_lst$Y, X=data_lst$X, prior=NULL, beta_init = NULL, MCSamplerSize =5500)

  MC_SIM_Reg_Dir_Data<-list(MC_lst=lst, Sim_Data_lst=data_lst)


  save(MC_SIM_Reg_Dir_Data, file=paste0("G:\\My Drive\\MyRPackages\\MC_SIM_Reg_Dir_Data\\MC_SIM_Reg_Dir_Data_d_eq_",d, "_SimNUmber_",SimNum, ".RData" ))

  print(paste0("sim Number", SimNum, "completed"))
  rm(list=c("MC_SIM_Reg_Dir_Data","lst" , "data_lst" ))

}
