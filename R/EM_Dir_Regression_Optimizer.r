

####################################################################################################################
####################################################################################################################

      ComputeD_ET<-function(betaMat_old, X, nu){
        intput_arg  <-  apply(X%*%betaMat_old, MARGIN = 1, FUN = norm)
        log_ET      <-  BesselIR(x=intput_arg, nu = nu, ifLog = TRUE)- log(2)- log(intput_arg)

        return(exp(log_ET))
      }

####################################################################################################################
####################################################################################################################

      EM_Dir_regression_optimizer_V1<-function(Y,
                                               X,
                                               prior=NULL,
                                               beta_init=NULL,
                                               EM_tolerence=.0001,
                                               Max_EM_iter=10000){

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

              nu= dim(Y)[2]/2-1
              secondPartOptEq<-(t(X)%*%Y+ prior$beta_mean/prior$beta_sigma_sq)

              #Initial Values for the loop
              if(is.null(beta_init)){    beta_init<- solve(t(X)%*%X + diag(p)/prior$beta_sigma_sq)%*%secondPartOptEq   }
              betaMat_old= beta_init;EMiter=1;diff=1000

            ###################################################################################################
            while(   ((EMiter<Max_EM_iter)&&(diff>EM_tolerence))  ){

                            ET<-ComputeD_ET(betaMat_old = betaMat_old, X = X, nu = nu)
                            betaMat_opt<- solve(2*t(X)%*%diag(ET)%*%X + diag(p)/prior$beta_sigma_sq)%*%secondPartOptEq

                      #Calculating the difference
                      diff<-max(abs(betaMat_opt-betaMat_old))#  norm(betaMat_opt-betaMat_old)
                      betaMat_old<-betaMat_opt
                      ###########################
                      EMiter=EMiter+1
                      print(EMiter)
            }
            #end While ########################################################################################

      return(betaMat_opt)
      }



####################################################################################################################
####################################################################################################################

