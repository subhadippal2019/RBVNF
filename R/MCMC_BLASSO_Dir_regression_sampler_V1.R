

sample_Tau_ij<-function(gig_chi, gig_psi){
  gig_chi= as.matrix(gig_chi)
  Tau_ij_all<- apply(X = gig_chi,
                      MARGIN = 1,
                      FUN = function(chi_x){
                                    ghyp::rgig(n = 1, lambda = 1/2, chi = chi_x, psi = gig_psi)
                                    }
                      )
return(Tau_ij_all)
}




Generate_Prior_Lambda_HyperGamma<-function(Lambda_lower, Lambda_upper, info_weight=1 ){
  if(info_weight<.5){info_weight=.5}
  approx_mean        <- ((Lambda_lower^2+Lambda_upper^2)/2)
  approx_sd          <- ((Lambda_upper^2-Lambda_lower^2)/(4*info_weight))
  Prior_Lambda_Hyper <- list(shape=approx_mean^2/approx_sd, rate=approx_mean/approx_sd)
  return(Prior_Lambda_Hyper)
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
#' lst=MCMC_BLASSO_Dir_regression_sampler_V1(Y=Y, X=X,  MCSamplerSize =100)
#' i=1;j= 1
#' library(ggplot2)
#' library(cowplot) #
#'Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))
#' @export
MCMC_BLASSO_Dir_regression_sampler_V1<-function(Y, X, prior=NULL,
                                                      beta_init=NULL,
                                                      Sigma_init=NULL,
                                                      MCSamplerSize=50,
                                                      K=100,
                                                      eps_accuracy=.00000001,
                                                      lasso_lambda=.01,
                                                      Sample_lasso_lambda=NULL # c(lambda_min, lambda_max, info_weight=1)
                                                      ){

  n=dim(Y)[1]; p=dim(X)[2]; d=dim(Y)[2]; nu=d/2-1
  names_X=colnames(X)
  #######################################################################################################
  #######################################################################################################
  ####################################Initial Value######################################################

  if(!is.null(Sample_lasso_lambda)){
          Prior_Lambda_Hyper<-Generate_Prior_Lambda_HyperGamma(Lambda_lower  = Sample_lasso_lambda[1],
                                                               Lambda_upper =Sample_lasso_lambda[2],
                                                               info_weight =Sample_lasso_lambda[3]  )
          lambda_prior_rate= Prior_Lambda_Hyper$rate
          lambda_prior_shape=Prior_Lambda_Hyper$shape
          Sample_lasso_lambda=TRUE
  }
  if(is.null(Sample_lasso_lambda)){
    Sample_lasso_lambda=FALSE}


  if(is.null(beta_init)){
    print("Default Procedure using EM is being used to obtain initial value of the regression coefficients that will be  used to start the MCMC Data Augmentation Algorithm. Iteration number of EM algorithm is being printed untill convergence." )
    #beta_init=EM_BLASSO_Dir_regression_optimizer_V1(Y=Y, X=X, beta_init = NULL,EM_tolerence = .001, Max_EM_iter = 500)
    beta_init=EM_BLASSO_Dir_regression_optimizer_V1(Y=Y, X=X, beta_init = NULL,EM_tolerence = .001, Max_EM_iter = 500, lasso_lambda = lasso_lambda)
    #EM_BLASSO_Dir_regression_optimizer_V1(Y=Y, X=X, beta_init = NULL, lasso_lambda = max(xx$lambda.1se),   EM_tolerence = .00001)
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
  Tau_ij_sq_all_store=NULL
  lasso_lambda_all=lasso_lambda
  #######################################################################################################
  #######################################################################################################

  sigma_square=1000
  print(" Initial value and prior information obtained successfully.  The MCMC samples are being generated. This step may take significnt amount of time depending on the MCMC sample size to be Generated.   " )
  #browser()
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


    #gig_chi=beta_jk^2; gig_psi= lambda^2

    #Tau_all=

    Tau_ij_sq_all= sample_Tau_ij(gig_chi=c(beta^2), gig_psi= lasso_lambda^2) #ghyp::rgig(n = 1, lambda = 1/2, chi = 1, psi = gig_psi)


    #GIGrgv::rgig(n = 10, lambda = 1, chi = 1, psi = 1)
    iMat_d<-diag(rep(1,d))
    D_tau_inv=diag(c(1/Tau_ij_sq_all))
    beta_post_Sigma<-solve( kronecker(iMat_d, 2*t(X)%*% Diag_T %*% X) +  D_tau_inv )
    beta_post_mean<- beta_post_Sigma %*% ( c(t(X)%*%Y ))

    beta_vec<-MASS::mvrnorm(n=1, mu =beta_post_mean , Sigma = beta_post_Sigma)

    beta=matrix(beta_vec, ncol=d)
   # beta=rMatrixNormal(M=beta_post_mean, V1 =beta_post_Sigma , V2 =iMat_d )

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

    if(Sample_lasso_lambda){
      lasso_lambda_post_rate= sum(Tau_ij_sq_all)/2+lambda_prior_rate
      lasso_lambda_post_shape= p*d+lambda_prior_shape
      lasso_lambda_sq=rgamma(n = 1, shape =lasso_lambda_post_shape , rate = lasso_lambda_post_rate )
      lasso_lambda= sqrt(lasso_lambda_sq)
      lasso_lambda_all=c(lasso_lambda_all,lasso_lambda )
    }


    ######################################################################################################
    ######################################################################################################
    ########################## Storing All Variables #####################################################
    ######################################################################################################
    beta_all[iter, ,   ]=beta
    #Theta_all[iter, , ]= Theta
    #sigma_square_all=c(sigma_square_all, sigma_square)
    #sigma_square_all[iter]=sigma_square
    T_aug_all=rbind(T_aug_all, T_aug)
    Tau_ij_sq_all_store= rbind(Tau_ij_sq_all_store, c(Tau_ij_sq_all))


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
                    T_aux_var=T_aug_all,
                    Tau_ij_sq_all= Tau_ij_sq_all_store,
                    lasso_lambda_all=lasso_lambda_all
                    )


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
                    eps_accuracy=eps_accuracy,
                    lasso_lambda=lasso_lambda_all)



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

