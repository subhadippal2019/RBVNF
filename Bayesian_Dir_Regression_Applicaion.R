#####################################################################################################
#####################################################################################################

n=1000 # NUmber of the samples
p=4  # NUmber of the regression covariates
d=5 # Number of direcions in the direcional data


#### bbeta is a matrix of dimension p\times d
#bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
sigma_square=1
tau_square=1000

#######################################################################################################
#######################################################################################################
#######################################################################################################

data_lst = Data_generator_vnf_reg(n=n, p=p, d=d, concentration_factor = 1, beta_factor = 5)

Y = data_lst$Y;X=data_lst$X;
#beta=data_lst$beta;
#Theta=data_lst$Theta

#######################################################################################################
#######################################################################################################
#######################################################################################################
########################### The MCMC Sampler ##########################################################
#######################################################################################################
tau_square=10000
prior=list()
prior$beta_0=matrix( 0, nrow=p, ncol=d)
prior$Sigma_0=tau_square*diag(p)        ############
prior$sigma_square_a=1
prior$sigma_square_b=1




# i=1
# sum_x_i_theta_i_est=X[i,]%*%t(Y[i, ])*0
# sum_var_mat= X[1,]%*%t(X[1,])*0
# for(i in 1:dim(Y)[1]){
#   sum_var_mat=sum_var_mat+ X[i,]%*%t(X[i,])
#   sum_x_i_theta_i_est=sum_x_i_theta_i_est+X[i,]%*%t(Y[i, ])
# }
#
# solve(sum_var_mat)%*%sum_x_i_theta_i_est
#
#
#
#
# sum_var_mat= X[1,]%*%t(X[1,])*0
# sum_Xi_theta_i=(as.matrix(X[i,]))%*%Theta[i,]*0
# for(i in 1:n){
#   sum_var_mat=sum_var_mat+ X[i,]%*%t(X[i,])
#   sum_Xi_theta_i=sum_Xi_theta_i+(as.matrix(X[i,]))%*%Theta[i,]
# }
# # beta_var= solve(sum_var_mat/sigma_square + solve(prior$Sigma_0))
# # beta_mean= beta_var %*% (     (sum_Xi_theta_i)/sigma_square + solve(prior$Sigma_0)%*%prior$beta_0        )
# beta_var= solve(sum_var_mat/sigma_square + solve(prior$Sigma_0))
# beta_mean= beta_var %*% (     (sum_Xi_theta_i)/sigma_square + solve(prior$Sigma_0)%*%prior$beta_0        )

#Theta_init=Theta;




#beta_init =data_lst$beta #solve(t(X)%*%X)%*%t(X)%*%Y







#

######################################################################################################
###################################### Data Generator ################################################
######################################################################################################
#source("C:\\Users\\subha\\Dropbox\\projects\\Regression of Directional data\\RCode\\functions.r")

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

data_lst = Data_generator_vnf_reg(n=n, p=p, d=d, concentration_factor = 1, beta_factor = 5)

Y = data_lst$Y;X=data_lst$X;

beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)
beta_init=beta_EM
lst=MCMC_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = beta_init,  MCSamplerSize =100)


i=2;j= 1
Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))


(data_lst$beta[i,j])

#(data_lst$beta)/matrix(lst$MC$Mc_Beta[1000, ], nrow=p)
#matrix(apply(lst$beta_all, MARGIN = c(2,3), FUN = mean), nrow=p)
apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
(data_lst$beta)
