

####################################################################################################################
####################################################################################################################

#ComputeD_ET<-function(betaMat_old, X, nu){
# intput_arg  <-  apply(X%*%betaMat_old, MARGIN = 1, FUN = norm)
#log_ET      <-  BesselIR(x=intput_arg, nu = nu, ifLog = TRUE)- log(2)- log(intput_arg)
#return(exp(log_ET))
#}

####################################################################################################################
####################################################################################################################

EM_BLASSO_Dir_regression_optimizer_V1_EMLAMBDA<-function(Y,
                                                X,
                                                beta_init=NULL,
                                                EM_tolerence=.0001,
                                                lasso_lambda=.1, # Lasso Penalty Parameter
                                                Max_EM_iter=10000,
                                                Convergence_loss=FALSE,
                                                scale_factor=2){
  n=dim(Y)[1]; p=dim(X)[2]; d=dim(Y)[2]; nu=d/2-1
  prior=NULL
  if(is.null(prior)){
    # tau_square=10000
    prior=list()
    prior$beta_mean=matrix( 0, nrow=p, ncol=d)
    #prior$beta_Sigma_0=tau_square*diag(p)        ############
    prior$beta_sigma_sq= 10000
    prior$sigma_square_shape=1
    prior$sigma_square_scale=1
    #prior$Sigma_Psi= matrix( 0, nrow=p, ncol=p)  #
    #prior$Sigma_nu=4  #
  }

  #nu= dim(Y)[2]/2-1


  secondPartOptEq<-(t(X)%*%Y+ prior$beta_mean/prior$beta_sigma_sq)

  #Initial Values for the loop
  if(is.null(beta_init)){    beta_init<- solve(t(X)%*%X + diag(p)/prior$beta_sigma_sq)%*%secondPartOptEq }
  betaMat_old= beta_init;EMiter=1;diff=1000
  diff_vec=c()
  #browser()
  ###################################################################################################
  while(   ((EMiter<Max_EM_iter)&&(diff>EM_tolerence))  ){
    ET<-ComputeD_ET(betaMat_old = betaMat_old, X = X, nu = nu)
    SQRT_D_ET_X= sqrt(ET)*X #sweep(X, MARGIN=1,sqrt(ET), `*`) #diag(sqrt(ET))%*% X #=
    #SQRT_D_ET_X= diag(sqrt(ET))%*% X


    I_d= diag(d)
    I_KRON_SQRT_D_ET_X= kronecker(I_d, SQRT_D_ET_X, FUN = "*")



    #NEG_SQRT_D_ET= diag(1/sqrt(ET))
    #I_KRON_NEG_SQRT_D_ET_Y= kronecker(I_d, NEG_SQRT_D_ET, FUN = "*")%*% as.vector(Y)/2
    I_KRON_NEG_SQRT_D_ET_Y= as.vector((1/sqrt(ET))*Y/2) # = NEG_SQRT_D_ET%*%Y/2
    I_KRON_NEG_SQRT_D_ET_Y_1= (1/sqrt(ET))*Y/2

    #mStandard_Lasso_fit = glmnet(x =SQRT_D_ET_X , y = I_KRON_NEG_SQRT_D_ET_Y_1, family = "mgaussian", lambda=lasso_lambda, alpha =1 )
    #mbetaMat_opt<- as.matrix(do.call(cbind, mStandard_Lasso_fit$beta))


    #Standard_Lasso_fit_hm=lasso_coord_desc(X=I_KRON_SQRT_D_ET_X,y=I_KRON_NEG_SQRT_D_ET_Y,beta=c(beta_init),lambda=lasso_lambda,tol=1e-6,maxiter=1000 )
    #betaMat_opt_sf_hm

    Standard_Lasso_fit = glmnet(x =I_KRON_SQRT_D_ET_X , y = I_KRON_NEG_SQRT_D_ET_Y, family = "gaussian", lambda=lasso_lambda, alpha = 1)
    #Standard_Lasso_fit1 = biglasso(X = as.big.matrix(I_KRON_SQRT_D_ET_X) ,y =  I_KRON_NEG_SQRT_D_ET_Y, family = "gaussian", lambda=lasso_lambda, alpha = 1)
    betaMat_opt<- matrix(Standard_Lasso_fit$beta, ncol=d)
    #solve(2*t(X)%*%diag(ET)%*%X + diag(p)/prior$beta_sigma_sq)%*%secondPartOptEq


    Sum_E_Tau_ij_sq = sum(abs(betaMat_opt)/lasso_lambda)+ p*d/lasso_lambda^2
          #Lasso_beta_EM=MCMC_BLASSO_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = betaMat_opt,  MCSamplerSize =200, lasso_lambda =lasso_lambda)
          #Sum_E_Tau_ij_sq=sum(apply(Lasso_beta_EM$MC$Tau_ij_sq_all, MARGIN = 2, mean))
    lasso_lambda= sqrt(scale_factor*p*d/(Sum_E_Tau_ij_sq)) #### sqrt(2*p*d/(Sum_E_Tau_ij_sq))


    print(lasso_lambda)
    #Calculating the difference
    diff<-max(abs(betaMat_opt-betaMat_old))#  norm(betaMat_opt-betaMat_old)
    diff_vec=c(diff_vec, diff)
    betaMat_old<-betaMat_opt
    ###########################
    EMiter=EMiter+1
    print(EMiter)
  }
  #end While ########################################################################################

  if(Convergence_loss){ return_lst<-list(betaMat_opt=betaMat_opt,conv=diff_vec) }
  if(!Convergence_loss){ return_lst=betaMat_opt }

  return( return_lst)
}



####################################################################################################################
####################################################################################################################

