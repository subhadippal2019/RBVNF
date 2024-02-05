




#####################################################################################################
########## Generic functions ########################################################################
#####################################################################################################
norm<-function(x){ return(sqrt( sum( x^2   ) )) }

r_Inv_Gamma<-function(n=1, shape, rate){
  return(1/rgamma(n=n, shape=shape, rate=rate))
}

Trace<-function(A){
  return(sum(diag(A)))
}


#####################################################################################################
############################## rvmf in the alternative parameterization #############################
#####################################################################################################


#' @export
rvmf_alt<-function(n=1, theta){
  #theta=as.vector(as.vector)
  #if(length(dim(theta)))
  kappa=sqrt(sum(theta^2))
 return( Rfast::rvmf(n=1, mu=theta/kappa, k=kappa))
}



#' @export
Data_generator_vnf_reg<-function(n=100,p=NULL, d=NULL,  beta=NULL, X=NULL, concentration_factor=1, beta_factor=1){
    #### beta is a  bbeta is a matrix of dimension p\times d matrix.
      if(!is.null(beta)){ p=dim(beta)[1];   d=dim(beta)[2];   }
      if(!is.null( X )){  n=dim(X)[1];   }

      #if( is.null(beta)){ beta =beta_factor*matrix(rnorm(p*d) ,nrow=p, ncol=d ) } # normally distributed X
      if( is.null(beta)){ beta =beta_factor*matrix( runif(p*d, -1, 1) ,nrow=p, ncol=d ) }  # Uniform(-10,10)
      if( is.null(X)){ X = array(  rnorm(n*p, 0, 1), dim=c(n,p)    ) }


    ##################  browser()


    Y= concentration_factor*matrix(replicate(n*d,  NA), nrow=n)
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




########################################################################################################
########################################################################################################
################################## The MCMC Engine #####################################################
########################################################################################################
########################################################################################################
#' @export
Sample_all_T_Aug<- function(n, nu, X, beta,K , j_nu_0 , J_nuPlus1  , eps_accuracy  ){
  T_aug=replicate( n,0)
  for(i in 1:n   ){
    x_i=X[i,] #### p regression variables for the i^{th} observarons
    theta_i = (     t(x_i) %*% beta     )
    kappa_i=norm_vec(theta_i)#sqrt( sum( theta_i^2   ) )
    #T_aug[i]=PG_randomSample(a = kappa_i,nu=nu, K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1 , eps_accuracy = eps_accuracy  )
   #PG_randomSample(a = kappa_i,nu=nu, K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1 , eps_accuracy = eps_accuracy  )
   T_aug[i]= rMPG(n=1, alpha=kappa_i, nu=nu)
    }
 # browser()
  return(T_aug)
}



#' Performs a Bayesian regression for the directional responses.
#'
#' @param Y:  n times d  matrix containing n directional response of dimension parameter d.
#' norm of each row of the matrix Y is 1.
#' @param X: A n*p design matrix containing the data of the p variables.
#' If it is required, the intercept should be included in the first column.
#' @param prior: Prior specification for the posterior. Default (prior=NULL) uses the very flat prior.
#' @param beta_init: Initial value of rregression coefficient. The default (beta_innit=NULL) utilizes the EM algorithm to get a initial value
#' @param MCSamplerSize: Number of MCMC samples required to be generated.
#' @return MCMC samples from for  Posterior of mode and concentration parameter of Vonmises distribution.
#' @examples
#' library(Rfast)
#' #library(benchmarkme)
#' library(Bessel)
#' library(gsl)
#' library(ghyp)
#' n=1000 # NUmber of the samples
#' p=4  # NUmber of the regression covariates
#' d=5 # Number of direcions in the direcional data
#' data_lst = Data_generator_vnf_reg(n=n, p=p, d=d, concentration_factor = 1, beta_factor = 5)
#' Y = data_lst$Y;X=data_lst$X;
#' lst=MCMC_Dir_regression_sampler_V1(Y=Y, X=X,  MCSamplerSize =100)
#' i=1;j= 1
#' library(ggplot2)
#' library(cowplot) #
#'Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))
#' @export
MCMC_Dir_regression_sampler_V1<-function(Y, X, prior=NULL, beta_init=NULL, Sigma_init=NULL, MCSamplerSize=50, K=100,eps_accuracy=.00000001){

  n=dim(Y)[1]; p=dim(X)[2]; d=dim(Y)[2]; nu=d/2-1
  #######################################################################################################
  #######################################################################################################
  ####################################Initial Value######################################################


  if(is.null(beta_init)){
    print("Default Procedure using EM is being used to obtain initial value of the regression coefficients that will be  used to start the MCMC Data Augmentation Algorithm. Iteration number of EM algorithm is being printed untill convergence." )
    beta_init=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,EM_tolerence = .00001)
  }
  else {

  if(!is.matrix(beta_init)){
    beta_init=NULL
    warning("beta_init: The initial value of beta is not matrix. The specification will be ignored and the function will initialize with the Default procedure.  ")
  }

  if(is.matrix(beta_init)){
    if(any(dim(beta_init) != c(p,d))){
      warning("beta_init: The initial value of beta is a matrix but not of a appropriate dimension. The specification will be ignored and the function will initialize with the Default procedure.  ")
      beta_init=NULL
    }
    if(any(!is.numeric(beta_init))){
      warning("beta_init: The initial value of beta is a matrix but not numeric. The specification will be ignored and the function will initialize with the Default procedure.  ")
      beta_init=NULL
    }
  }
}



  if(is.null(Sigma_init)){
        tau_square=10000
        Sigma=tau_square*diag(p)
  }
  ####################################Initial Assignment Step############################################
  beta=beta_init
  #######################################################################################################
  #######################################################################################################
  #################################### Prior Specification ##############################################



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



  Start_Time=Sys.time()
  #K = 100
  j_nu_0=gsl::bessel_zero_Jnu(nu,1:K);
  J_nuPlus1=BesselJ(j_nu_0,(nu+1));


  ########################### Storing all Values ########################################################
  #######################################################################################################
  beta_all=array(data = NA, dim=c(MCSamplerSize, p, d) )
  #Sigma_all<-array(data = NA, dim=c(MCSamplerSize, p, p) )
  #sigma_square_all=NULL
  T_aug_all=NULL
  #######################################################################################################
  #######################################################################################################

  sigma_square=1000
  print(" Initial value and prior information obtained successfully.  The MCMC samples are being generated. This step may take significnt amount of time depending on the MCMC sample size to be Generated.   " )
  #browser()
  for(iter in 1:MCSamplerSize){
    #######################################################################################################
    ######### Sample of the Augmented Variable ############################################################

    T_aug<-Sample_all_T_Aug(n=n, nu=nu, X=X, beta=beta,K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1 , eps_accuracy = eps_accuracy  )


    ######################################################################################################
    ####################Sample of beta  ###########################################################
    ######################################################################################################

              Diag_T<- diag(T_aug)
              #beta_post_Sigma<-solve(2*t(X)%*% Diag_T %*% X +solve(prior$beta_Sigma_0)  )
              #beta_post_mean<- beta_post_Sigma %*% (t(X)%*%Y+ solve(prior$beta_Sigma_0)%*%prior$beta_mean)
              Sigma=sigma_square*diag(p)
              beta_post_Sigma<-solve(2*t(X)%*% Diag_T %*% X +solve(Sigma)  )
              beta_post_mean<- beta_post_Sigma %*% (t(X)%*%Y+ solve(Sigma)%*%prior$beta_mean)
              iMat_d<-diag(rep(1,d))
              beta=rMatrixNormal(M=beta_post_mean, V1 =beta_post_Sigma , V2 =iMat_d )

   ######################################################################################################
   ####################Sample of Sigma  ###########################################################
   ######################################################################################################
              #post_sigma_square_shape<-prior$sigma_square_shape+p*d/2
              #post_Sigma_Psi<-prior$Sigma_Psi+ (beta-prior$beta_mean)%*%t(beta-prior$beta_mean)
              #post_Sigma_nu<-prior$Sigma_nu+d
              #post_sigma_square_scale<-prior$sigma_square_scale+.5*Trace(t(beta-prior$beta_mean)%*%(beta-prior$beta_mean))
              #post_Sigma_nu<-prior$Sigma_nu+d
              #sigma_square<- r_Inv_Gamma(n = 1, shape =post_sigma_square_shape, rate =  post_sigma_square_scale)
                 #print(sigma_square)
                 # Sigma<-MCMCpack::riwish(v=post_Sigma_nu, S=post_Sigma_Psi)

    ######################################################################################################
    ######################################################################################################
    ########################## Storing All Variables #####################################################
    ######################################################################################################
    beta_all[iter, ,   ]=beta
    #Theta_all[iter, , ]= Theta
    #sigma_square_all=c(sigma_square_all, sigma_square)
    #sigma_square_all[iter]=sigma_square
    T_aug_all=rbind(T_aug_all, T_aug)


    ######################################################################################################
    ######################################################################################################
    ########################## Print Statement ###########################################################
    ######################################################################################################
    if(iter%%100==0){
           print(paste0("MC_Iter=", iter,"completed"))
    }

  } #### The End of the MCMC





  ## Codes for storing the pivotal quantities

  Run_Time=Sys.time()-Start_Time
  #Sys_info=get_sys_details()


  ###################### put this chunk inside a function
  function_name=as.character(match.call()[1])
  function_def0<-  capture.output(print(get( function_name)));
  function_def0[1]=paste0( function_name,"<-",function_def0[1])
  function_def=paste( function_def0 , collapse = "\n")
  #function_def1=dput(get(function_name))
  #########

  Sys_info=Sys.info()

  # MC<-list(Mc_Beta=Store_Beta,
  # Mc_sigma_sq=Store_sigma_sq,
  # Mc_xi=Store_xi,
  # Mc_eta=Store_eta,
  # Mc_Nu=Store_Nu)


  MC<-list(         Mc_Beta=  beta_all,
                    T_aux_var=T_aug_all)


  RunDetails<- list(Data=list(Y=Y, X=X),
                 prior=prior,
                 beta_init=beta_init,
                 Methodtype=paste0("Data Augmentation Gibbs for regression of directional vectors in d= ", d, "dimensional Euclidean Space, i.e. Sphere d-1"),
                 Run_Time=Run_Time,
                 System_info=Sys_info,
                 call=match.call(),
                 function_def=function_def,
                 additionalInfo="Utilize #cat(lst$RunDetails$function_def) to see the function definition properly.\\ MC contains the main Markov chain on beta and augmented variables.\\ K is a technical parameter.It is the Truncation term of a MPG density series. ",
                 MCSamplerSize=MCSamplerSize,
                 K=K,
                 eps_accuracy=eps_accuracy)



  lst<-list(MC=MC,
            RunDetails=RunDetails
            )


  #Y, X, prior, beta_init_vec, Sigma_init=NULL, MCSamplerSize=50, K=100,eps_accuracy=.00000001


#cat(lst$function_def)


  #lst=list();
  #lst$beta_all=beta_all
  #lst$Theta_all=Theta_all
  #lst$sigma_square_all=sigma_square_all
  #lst$T_aug_all=T_aug_all
  return(lst)
}












### The old Formulation
#' @export
MCMC_Dir_regression_sampler<-function(Y, X, prior, Theta_init, beta_init_vec, mc_length=50){

  n=dim(Y)[1]; p=dim(X)[2]; d=dim(Y)[2]
  #######################################################################################################
  #######################################################################################################
  ####################################Initial Value######################################################

  beta=beta_init_vec
  Theta=Theta_init

  ########################### Storing all Values ########################################################
  #######################################################################################################
  beta_all=NULL
  sigma_square_all=NULL
  Theta_all = array(NA, dim=c(mc_length, n, d))
  Trt_all=NULL
  #######################################################################################################
  #######################################################################################################

  sigma_square=1
  for(iter in 1:mc_length){
                          #######################################################################################################
                          ######### Sample of the Augmented Variable ############################################################
                          Trt=replicate( n,0)
                          for(i in 1:dim(Theta)[1]   ){
                              kappa_post=norm(Theta[i,])#sqrt( sum( theta_i^2   ) )
                              Trt[i]=PG_randomSample(a =kappa_post  ,nu=d/2-1  );
                          }
                          #######################################################################################################
                          ######################################### sample of theta_i    ########################################
                          #######################################################################################################

                          for(i in 1:dim(Theta)[1]){
                              theta_i_hat<-t((sigma_square*Y[i,]+t(X[i,])%*%beta)/(1+2*Trt[i]*sigma_square))
                              theta_i_var= diag(d)*(sigma_square/(1+2*Trt[i]*sigma_square)) ### Need to change this line once the general variance components compes
                              theta_i= rmvnorm(n = 1, mean = theta_i_hat, sigma =theta_i_var )#rnorm( 1, theta_i_hat,  theta_i_var)
                              Theta[i, ]= theta_i
                          }

                          ######################################################################################################
                          ####################Sample of sigma Square ###########################################################
                          ######################################################################################################
                          #sum_thet_i_minus_x_i_beta=0
                          #for( i in 1:n){
                          #  sum_thet_i_minus_x_i_beta=sum_thet_i_minus_x_i_beta+(Theta[i, ]-X[i,]%*%beta)%*%t(Theta[i, ]-X[i,]%*%beta)
                          #}

                          #post_sigma_square_a=n*d+prior$sigma_square_a
                          #post_sigma_square_b=sum_thet_i_minus_x_i_beta/2+prior$sigma_square_b

                          #sigma_square=r_Inv_Gamma(n=1, shape=post_sigma_square_a, rate=post_sigma_square_b)
                          ######################################################################################################
                          ####################Sample of Beta ###################################################################
                          ######################################################################################################

                          #sigma_square=1; tau_square=1
                          ## Sample of
                          # sum_var_mat =
                            #n=dim(X)[1]
                            #p=dim(X)[2]
                            sum_var_mat= X[1,]%*%t(X[1,])*0
                            sum_Xi_theta_i=(as.matrix(X[i,]))%*%Theta[i,]*0
                            for(i in 1:n){
                              sum_var_mat=sum_var_mat+ X[i,]%*%t(X[i,])
                              sum_Xi_theta_i=sum_Xi_theta_i+(as.matrix(X[i,]))%*%Theta[i,]
                            }
                           # beta_var= solve(sum_var_mat/sigma_square + solve(prior$Sigma_0))
                           # beta_mean= beta_var %*% (     (sum_Xi_theta_i)/sigma_square + solve(prior$Sigma_0)%*%prior$beta_0        )
                            beta_var= solve(sum_var_mat/sigma_square + solve(prior$Sigma_0))
                            beta_mean= beta_var %*% (     (sum_Xi_theta_i)/sigma_square + solve(prior$Sigma_0)%*%prior$beta_0        )


                            # sum_var_mat
                            #vec_beta_var= solve(kronecker( diag(d),  solve(prior$Sigma_0) ) + kronecker(sum_var_mat, diag(d) )  )
                            vec_beta_var= solve(kronecker(solve(prior$Sigma_0), diag(d) ) + kronecker(diag(d), sum_var_mat )  )
                            #browser()
                            vec_beta=mvtnorm::rmvnorm(n=1, mean = as.vector(beta_mean), sigma =  vec_beta_var)
                            sampled_beta=matrix(vec_beta, nrow=p)




                            ######################################################################################################
                            ######################################################################################################
                            ########################## Storing All Variables #####################################################
                            ######################################################################################################
                            beta_all=rbind(beta_all,vec_beta)
                            #Theta_all[iter, , ]= Theta
                            sigma_square_all[iter]=sigma_square
                            Trt_all=rbind(Trt_all, Trt)


                            ######################################################################################################
                            ######################################################################################################
                            ########################## Print Statement ###########################################################
                            ######################################################################################################
                            print(paste0("MC_Iter=", iter,"completed"))

                      } #### The End of the MCMC

  lst=list();
  lst$beta_all=beta_all
  #lst$Theta_all=Theta_all
  lst$sigma_square_all=sigma_square_all
  lst$Trt_all=Trt_all
  return(lst)
}






Plot_Estimated_output_vs_Observed<-function(lst, data_lst){

  Y_est=X%*%beta_est
  Y_est_norm=t(apply(Y_est, 1, function(x){ x/norm(x) }))

  par(mfrow=c(d, 1))
  for( i  in 1:d){
  plt=plot( data_lst$Y[,i ], Y_est_norm[,i])
  }

}


#library(gsl)

Del_f<-function(y,theta, nu){
 (1/norm(theta))* besselI(norm(theta), nu=nu+1)/ besselI(norm(theta), nu=nu)*theta -y
}



Theta_init<-function(Y){


  d=dim(Y)[2]
  norm_theta=norm(Y[i, ])
  besselI(x = norm_theta, nu = d/2-1)
  #for solving f=0
  for( NR_iter in 1:100){
    theta_new<-theta- theta_current-solve(Del_f)%*% Del_f
    Theta_store=rbind(Theta_store, theta_new )
  }


}











