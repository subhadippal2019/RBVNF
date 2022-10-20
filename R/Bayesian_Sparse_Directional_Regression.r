


#' @export
Mat_inv<-function(M){
  return(chol2inv(chol(M)))
  #return(solve(M))
  #return(pd.solve(M))
}

#Data_sim<-GererateData(n=100, p=500, SigmaTrue = .1, NumOfNonZeroBeta = 10)
#
# GererateData<-function(n=70,p=100,SigmaTrue=2, NumOfNonZeroBeta=5 ){
#   #Data Generation #n=70#p = 100  #%###### I have added
#   SigmaTrue = 2;
#   BetaTrue = rep(0, p)
#   #BetaTrue(1:23) = 2.^(-(-2:.25:3.5))
#   #BetaTrue[1:NumOfNonZeroBeta] = 2^(-seq(-2,3.5, by = .25))
#   BetaTrue[1:NumOfNonZeroBeta] = 2^(-seq(-3,3.5,length.out =NumOfNonZeroBeta ))
#   TauTrue = 1;
#
#
#   X = matrix(rnorm(n*p, 0,1), ncol=p, nrow=n);
#   y = X%*%BetaTrue+SigmaTrue*rnorm(n,0,1);
#   SimData=list(X=X, y=y, BetaTrue=BetaTrue, TauTrue=TauTrue, SigmaTrue=SigmaTrue)
#   return(SimData)
# }




#' @export
Data_generator_vnf_reg_sparse<-function(n=100,p=30, d=3,NumOfNonZeroBeta=c(5, 4, 6),  beta=NULL, X=NULL){
  #### beta is a  bbeta is a matrix of dimension p\times d matrix.
  if(!is.null(beta)){ p=dim(beta)[1];   d=dim(beta)[2];   }
  if(!is.null( X )){  n=dim(X)[1];   }

  if( is.null(beta)){
    beta =matrix(0 ,nrow=p, ncol=d )
            for(jj in 1:d){
                beta[1:NumOfNonZeroBeta[jj], jj] = (2+.5/jj)^(-seq(-3,3.5,length.out =NumOfNonZeroBeta[jj] ))
            }
    }
  if( is.null(X)){ X = array(  rnorm(n*p, 0, 1), dim=c(n,p)    ) }
  #browser()rm(list=ls())

  ##################  browser()


  Y= matrix(replicate(n*d,  NA), nrow=n)
  for(i in 1:n){
    x_i=X[i,] #### p regression variables for the i^{th} observarons
    theta_i = (     t(x_i) %*% beta     )
    #Theta[i, ]=theta_i
    y_i =  rvmf_alt(n=1, theta=as.vector(theta_i ))
    Y[i, ]=y_i
  }
  data_lst=list(Y=Y, X=X, beta=beta)
  return(data_lst)
}






#' @export
MCMC_Dir_regression_sampler_sparse_V1<-function(Y, X, prior, beta_init_vec, Sigma_init=NULL, MCSamplerSize=50, K=100,eps_accuracy=.00000001){

  n=dim(Y)[1]; p=dim(X)[2]; d=dim(Y)[2]; nu=d/2-1
  #######################################################################################################
  #######################################################################################################
  ####################################Initial Value######################################################

  beta=beta_init_vec
  if(is.null(Sigma_init)){
    tau_square=10000;
    Sigma=tau_square*diag(p)
  }
  #Theta=Theta_init

  if(is.null(prior)){
    # tau_square=10000
    prior=list()
    prior$beta_mean=matrix( 0, nrow=p, ncol=d)
    prior$beta_Sigma_0=tau_square*diag(p)        ############
    prior$sigma_square_shape=1
    prior$sigma_square_scale=1
    #prior$Sigma_Psi= matrix( 0, nrow=p, ncol=p)  #
    #prior$Sigma_nu=4  #
  }





  #K = 100

  j_nu_0=gsl::bessel_zero_Jnu(nu,1:K);
  J_nuPlus1=BesselJ(j_nu_0,(nu+1));


  ########################### Storing all Values ########################################################
  #######################################################################################################
  beta_all=array(data = NA, dim=c(MCSamplerSize, p, d) )
  lambda_all=array(data = NA, dim=c(MCSamplerSize, p, d) )
  tau_all=NULL
  #Sigma_all<-array(data = NA, dim=c(MCSamplerSize, p, p) )
  #sigma_square_all=NULL
  T_aug_all=NULL
  #######################################################################################################
  #######################################################################################################

  sigma_square=1000
  tau_vec=replicate(n = d, 1);  lambda_mat=matrix( 1, nrow = p, ncol = d)



 # browser()
  for(iter in 1:MCSamplerSize){

    ############################### Sampling of shrinkage parameters #################################
    ##################################################################################################

    for(tau_id in 1:d){
      b_par<-.5*sum((beta[,tau_id ])^2/(lambda_mat[,tau_id])^2)
      tau_vec[tau_id]<-rIEHC(n = 1, a =p , b =b_par )
    }

    for(lamda_row_id in 1:p){ #j
      for(lamda_col_id in 1:d){#k
        b_par<-.5*( beta[lamda_row_id,lamda_col_id ])^2 / tau_vec[lamda_col_id]
        if(b_par<.0000000001){b_par=.0000000001}
        lambda_mat[lamda_row_id, lamda_col_id]<-rIEHC(n = 1, a =1 , b =b_par   )
      }
    }



    #######################################################################################################
    ######### Sample of the Augmented Variable ############################################################

    T_aug<-Sample_all_T_Aug(n=n, nu=nu, X=X, beta=beta,K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1 , eps_accuracy = eps_accuracy  )


    ######################################################################################################
    ####################Sample of beta  ###########################################################
    ######################################################################################################

    Diag_T<- diag(T_aug)

    # need tau_vec and lambda_mat
    D_lambda_tau=(diag(tau_vec)%x% diag(p))%*% diag(c(lambda_mat))
    #beta_post_Sigma<-solve(2*t(X)%*% Diag_T %*% X +solve(prior$beta_Sigma_0)  )
    #beta_post_mean<- beta_post_Sigma %*% (t(X)%*%Y+ solve(prior$beta_Sigma_0)%*%prior$beta_mean)
    #Sigma=sigma_square*diag(p)
   # beta_post_Sigma <-  solve( 2*t(X)%*% Diag_T %*% X )
    iMat_d<-diag(rep(1,d)); y_vec=c(Y)
    beta_vec_Sigma<-   solve( iMat_d%x%( 2*t(X)%*% Diag_T %*% X) + solve(D_lambda_tau))
    beta_post_mean<- beta_vec_Sigma %*% (iMat_d%x%t(X))%*%y_vec


    #beta_post_mean<- beta_post_Sigma %*% (t(X)%*%Y+ solve(Sigma)%*%prior$beta_mean)
    beta_vec<-mvrnorm(n=1,beta_post_mean,beta_vec_Sigma)
    beta<-matrix(beta_vec, nrow = p, ncol = d)
    #beta=rMatrixNormal(M=beta_post_mean, V1 =beta_post_Sigma , V2 =iMat_d )






    #################################################################################
    #################################################################################

    ######################################################################################################
    ######################################################################################################
    ########################## Storing All Variables #####################################################
    ######################################################################################################
    beta_all[iter, ,   ]=beta
    lambda_all[iter, ,   ]=lambda_mat
    #Theta_all[iter, , ]= Theta
    #sigma_square_all=c(sigma_square_all, sigma_square)
    #sigma_square_all[iter]=sigma_square
    T_aug_all=rbind(T_aug_all, T_aug)
    tau_all=rbind(tau_all, tau_vec)


    ######################################################################################################
    ######################################################################################################
    ########################## Print Statement ###########################################################
    ######################################################################################################
    if(iter%%100==0){
      print(paste0("MC_Iter=", iter,"completed"))
    }

  } #### The End of the MCMC

  lst=list();
  lst$data=list(X=X, Y=Y)
  lst$beta_all=beta_all
  #lst$Theta_all=Theta_all
  #lst$sigma_square_all=sigma_square_all
  lst$T_aug_all=T_aug_all
  lst$tau_all=tau_all
  lst$lambda_all=lambda_all
  return(lst)
}









#
#
# Del_f<-function(y,theta, nu){
#   (1/norm(theta))* besselI(norm(theta), nu=nu+1)/ besselI(norm(theta), nu=nu)*theta -y
# }
#
#
#
# Theta_init<-function(Y){
#
#
#   d=dim(Y)[2]
#   norm_theta=norm(Y[i, ])
#   besselI(x = norm_theta, nu = d/2-1)
#   #for solving f=0
#   for( NR_iter in 1:100){
#     theta_new<-theta- theta_current-solve(Del_f)%*% Del_f
#     Theta_store=rbind(Theta_store, theta_new )
#   }
#
#
# }
#
#



#' @export
MCMC_Dir_regression_HS_n_less_p<-function(Y, X, prior, beta_init_vec, Sigma_init=NULL, MCSamplerSize=50, K=100,eps_accuracy=.00000001){

  n=dim(Y)[1]; p=dim(X)[2]; d=dim(Y)[2]; nu=d/2-1
  #######################################################################################################
  #######################################################################################################
  ####################################Initial Value######################################################

  beta=beta_init_vec
  if(is.null(Sigma_init)){
    tau_square=10000;
    Sigma=tau_square*diag(p)
  }
  #Theta=Theta_init

  if(is.null(prior)){
    # tau_square=10000
    prior=list()
    prior$beta_mean=matrix( 0, nrow=p, ncol=d)
    prior$beta_Sigma_0=tau_square*diag(p)        ############
    prior$sigma_square_shape=1
    prior$sigma_square_scale=1
    #prior$Sigma_Psi= matrix( 0, nrow=p, ncol=p)  #
    #prior$Sigma_nu=4  #
  }





  #K = 100

  j_nu_0=gsl::bessel_zero_Jnu(nu,1:K);
  J_nuPlus1=BesselJ(j_nu_0,(nu+1));


  ########################### Storing all Values ########################################################
  #######################################################################################################
  beta_all=array(data = NA, dim=c(MCSamplerSize, p, d) )
  lambda_all=array(data = NA, dim=c(MCSamplerSize, p, d) )
  tau_all=NULL
  #Sigma_all<-array(data = NA, dim=c(MCSamplerSize, p, p) )
  #sigma_square_all=NULL
  T_aug_all=NULL
  #######################################################################################################
  #######################################################################################################

  sigma_square=1000
  tau_vec=replicate(n = d, 1);  lambda_mat=matrix( 1, nrow = p, ncol = d)



  #browser()
  for(iter in 1:MCSamplerSize){

    ############################### Sampling of shrinkage parameters #################################
    ##################################################################################################

    for(tau_id in 1:d){
      b_par<-.5*sum((beta[,tau_id ])^2/(lambda_mat[,tau_id])^2)
      tau_vec[tau_id]<-rIEHC(n = 1, a =p , b =b_par )
    }

    for(lamda_row_id in 1:p){ #j
      for(lamda_col_id in 1:d){#k
        b_par<-.5*( beta[lamda_row_id,lamda_col_id ])^2 / tau_vec[lamda_col_id]
        if(b_par<.0000000001){b_par=.0000000001}
        lambda_mat[lamda_row_id, lamda_col_id]<-rIEHC(n = 1, a =1 , b =b_par   )
      }
    }



    #######################################################################################################
    ######### Sample of the Augmented Variable ############################################################

    T_aug<-Sample_all_T_Aug(n=n, nu=nu, X=X, beta=beta,K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1 , eps_accuracy = eps_accuracy  )


    ######################################################################################################
    ####################Sample of beta  ###########################################################
    ######################################################################################################

    Diag_T<- diag(T_aug); X_star=sqrt(2)*diag(sqrt(T_aug))%*% X

    # need tau_vec and lambda_mat
    D_lambda_tau=(diag(tau_vec)%x% diag(p))%*% diag(c(lambda_mat))
    #beta_post_Sigma<-solve(2*t(X)%*% Diag_T %*% X +solve(prior$beta_Sigma_0)  )
    #beta_post_mean<- beta_post_Sigma %*% (t(X)%*%Y+ solve(prior$beta_Sigma_0)%*%prior$beta_mean)
    #Sigma=sigma_square*diag(p)
    # beta_post_Sigma <-  solve( 2*t(X)%*% Diag_T %*% X )
    #y_vec_star= ( diag(d) %x% diag(2*sqrt(T_aug) ) )%*%c(Y)
    #Phi= diag(d) %x% X_star;
    #u=rnorm(p*d)*sqrt(diag(D_lambda_tau));
    #v=Phi%*%u+rnorm(n*d);
    #v_star=(Mat_inv(Phi%*%D_lambda_tau%*%t(Phi) + diag(n*d)  ))  %*%((y_vec_star-v)); # alpha=y/sqrt(sigma_sq)
    #beta_vec=(u+ D_lambda_tau %*%(diag(d) %x%t(X))%*%v_star)  # v_star=w



    #################



    iMat_d<-diag(rep(1,d)); y_vec=c(Y)
    Phi=  (iMat_d%x% X_star)
    y_vec_star= ( diag(d) %x% diag(.5/sqrt(T_aug) ) )%*%c(Y)

    U_1=rnorm(p*d)*sqrt(diag(D_lambda_tau));
    V_1=Phi%*%U_1+rnorm(d*n);
    V_2=(Mat_inv( Phi%*% D_lambda_tau %*%t(Phi) +diag(d*n) ))  %*%( y_vec_star -V_1); # alpha=y/sqrt(sigma_sq)
    beta_vec=(U_1+D_lambda_tau%*%t(Phi)%*%V_2)  # v_star=w



    #beta_vec_Sigma<-   solve( iMat_d%x%( 2*t(X)%*% Diag_T %*% X) + solve(D_lambda_tau))
    #beta_post_mean<- beta_vec_Sigma %*% (iMat_d%x%t(X))%*%y_vec


    #beta_post_mean<- beta_post_Sigma %*% (t(X)%*%Y+ solve(Sigma)%*%prior$beta_mean)
    #beta_vec<-mvrnorm(n=1,beta_post_mean,beta_vec_Sigma)
    beta<-matrix(beta_vec, nrow = p, ncol = d)
    #beta=rMatrixNormal(M=beta_post_mean, V1 =beta_post_Sigma , V2 =iMat_d )






    #################################################################################
    #################################################################################

    ######################################################################################################
    ######################################################################################################
    ########################## Storing All Variables #####################################################
    ######################################################################################################
    beta_all[iter, ,   ]=beta
    lambda_all[iter, ,   ]=lambda_mat
    #Theta_all[iter, , ]= Theta
    #sigma_square_all=c(sigma_square_all, sigma_square)
    #sigma_square_all[iter]=sigma_square
    T_aug_all=rbind(T_aug_all, T_aug)
    tau_all=rbind(tau_all, tau_vec)


    ######################################################################################################
    ######################################################################################################
    ########################## Print Statement ###########################################################
    ######################################################################################################
    if(iter%%100==0){
      print(paste0("MC_Iter=", iter,"completed"))
    }

  } #### The End of the MCMC

  lst=list();
  lst$data=list(X=X, Y=Y)
  lst$beta_all=beta_all
  #lst$Theta_all=Theta_all
  #lst$sigma_square_all=sigma_square_all
  lst$T_aug_all=T_aug_all
  lst$tau_all=tau_all
  lst$lambda_all=lambda_all
  return(lst)
}









