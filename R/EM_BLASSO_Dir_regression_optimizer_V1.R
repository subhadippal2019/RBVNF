

####################################################################################################################
####################################################################################################################

#ComputeD_ET<-function(betaMat_old, X, nu){
 # intput_arg  <-  apply(X%*%betaMat_old, MARGIN = 1, FUN = norm)
  #log_ET      <-  BesselIR(x=intput_arg, nu = nu, ifLog = TRUE)- log(2)- log(intput_arg)
  #return(exp(log_ET))
#}

####################################################################################################################
####################################################################################################################

#' @export
EM_BLASSO_Dir_regression_optimizer_V1<-function(Y,
                                         X,
                                         beta_init=NULL,
                                         EM_tolerence=.0001,
                                         lasso_lambda=.1, # Lasso Penalty Parameter
                                         Max_EM_iter=10000,
                                         Convergence_loss=FALSE,
                                         if_Print=TRUE){
  prior=NULL
  n=dim(Y)[1]; p=dim(X)[2]; d=dim(Y)[2]; nu=d/2-1
  #browser()
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


    #Calculating the difference
    diff<-max(abs(betaMat_opt-betaMat_old))#  norm(betaMat_opt-betaMat_old)
    diff_vec=c(diff_vec, diff)
    betaMat_old<-betaMat_opt
    ###########################
    if(sum(betaMat_opt^2)<0.000000001){ EMiter=Max_EM_iter}
    EMiter=EMiter+1
    if(if_Print){print(EMiter)}
  }
  #end While ########################################################################################

  if(Convergence_loss){ return_lst<-list(betaMat_opt=betaMat_opt,conv=diff_vec) }
  if(!Convergence_loss){ return_lst=betaMat_opt }

  return( return_lst)
}



####################################################################################################################
####################################################################################################################
##### Codes for the Cross Validation ######


Predict_Y <- function(Beta_hat, X_new){
  Y_unnorm  <- X_new%*% (Beta_hat)
  Y_hat     <- t(apply(Y_unnorm, MARGIN =1,  function(y){return(y/sqrt(sum(y^2)))} ))
  return(Y_hat)
}

Compute_Error<-function(Y_true, Y_hat){
  n= dim(Y_true)[1]
  #browser()
  #Y_hat= Predict_Y(Beta_hat, X_new=X_true)
  RMSE= sqrt(sum((Y_true - Y_hat)^2)/n)
  return(RMSE)
}



#########################
#########################


#### New more efficient function
Cross_validation_error_all_lambda<-function(Y=Y, X=X, cv_lasso_lambda=NULL, k_fold=10,  if_Print=FALSE, Max_EM_iter=10000, strategy_1=TRUE){
  n= dim(Y)[1]
  col_Num=dim(Y)[2]
  X_col_Num= dim(X)[2]
  k=k_fold
  N_lambda=length(cv_lasso_lambda)
  Err=matrix(data=NA,nrow=N_lambda, ncol=k  )
  Y_input=Y; X_input=X



  for(i in 1:k ){
    #sel_rows<-(c(1:n))[-i]
    sample_per_set=as.integer(n/k)
    start=sample_per_set*(i-1)+1
    end=sample_per_set*(i)
    sel_rows<-start:end
    Y_train= Y[sel_rows, ]; X_train= X[sel_rows, ]
    Y_test= matrix(Y[-sel_rows, ],ncol = col_Num); X_test= matrix(X[-sel_rows, ], ncol=X_col_Num)
    beta_init=NULL
    #browser()
    for(lambda_id in 1:N_lambda){
      if(lambda_id>1){
        if(strategy_1){
          beta_init=Beta_est
        }
      }
      Beta_est<-EM_BLASSO_Dir_regression_optimizer_V1(Y=Y_train,
                                                      X=X_train,
                                                      beta_init=beta_init,
                                                      EM_tolerence=.0001,
                                                      lasso_lambda=cv_lasso_lambda[lambda_id], # Lasso Penalty Parameter
                                                      Max_EM_iter=Max_EM_iter,
                                                      Convergence_loss=FALSE,
                                                      if_Print=if_Print
      )

      Y_hat= Predict_Y(Beta_hat =Beta_est, X_new = X_test )
      Err[lambda_id, i]= Compute_Error(Y_true =Y_test,Y_hat=Y_hat  )
    }
    print(paste0("Cross Validated Error computation is completed for all lambda and for the Testing Fold Number= ",i ))

  }
  rownames(Err)=paste0("Lambda_", 1:length(cv_lasso_lambda))
  colnames(Err)=paste0("Fold_", 1:k)
  #browser()
  return((Err))
}



#system.time( Cross_validation_error_all_lambda(Y=Y, X=X,cv_lasso_lambda=c(0.01, 0.0100) ,k_fold = 2, if_Print=FALSE) )







#########################
#########################







#' @export
Cross_validation_error<-function(Y=Y, X=X, lasso_lambda=.001, k_fold=4,replication=1,  if_Print=FALSE, Max_EM_iter=10000){
  n= dim(Y)[1]
  col_Num=dim(Y)[2]
  X_col_Num= dim(X)[2]
  k=k_fold
  Err=matrix(data=NA,nrow=replication, ncol=k  )
  Y_input=Y; X_input=X

  for(repIter in 1:replication){
    if(repIter>1){
         RandPerm= sample(1:n, replace=FALSE)
         Y=Y_input[RandPerm, ]; X= X_input[RandPerm, ]
    }

              for(i in 1:k ){
                #sel_rows<-(c(1:n))[-i]
                sample_per_set=as.integer(n/k)
                start=sample_per_set*(i-1)+1
                end=sample_per_set*(i)
                sel_rows<-start:end
                Y_train= Y[sel_rows, ]; X_train= X[sel_rows, ]
                Y_test= matrix(Y[-sel_rows, ],ncol = col_Num); X_test= matrix(X[-sel_rows, ], ncol=X_col_Num)

                Beta_est<-EM_BLASSO_Dir_regression_optimizer_V1(Y=Y_train,
                           X=X_train,
                           beta_init=NULL,
                           EM_tolerence=.0001,
                           lasso_lambda=lasso_lambda, # Lasso Penalty Parameter
                           Max_EM_iter=Max_EM_iter,
                           Convergence_loss=FALSE,
                           if_Print=if_Print)

               Y_hat= Predict_Y(Beta_hat =Beta_est, X_new = X_test )
               Err[repIter, i]= Compute_Error(Y_true =Y_test,Y_hat=Y_hat  )
              }
  }
  return((Err))
}











#X_new= (data_lst$X[1:5, ]);Y_true= (data_lst$Y[1:5, ])


########### Function for Cross Validation #####################

#' @export
EM_BLASSO_Dir_regression_optimizer_V0.cv<-function(Y,
                                                   X,
                                                   beta_init=NULL,
                                                   EM_tolerence=.0001,
                                                   cv_lasso_lambda=NULL, # Lasso Penalty Parameter for cross validation
                                                   Max_EM_iter=10000,
                                                   Convergence_loss=FALSE,
                                                   if_Print=TRUE,
                                                   cv_k_fold=10,
                                                   cv_lambda_n=100,
                                                   epsilon_lambda_range_min=.0001){


  #browser()
  if(is.null(cv_lasso_lambda)){
    #lso<-cv.glmnet(x=kronecker(diag(d), X), y=c(Y))
    #Max_lambda<-max(lso$lambda)
    xx_kron<-kronecker(diag(d), X);y_vec=c(Y)
    #Max_lambda<-min( apply(Y, MARGIN = 2, FUN = function(yy){ max(abs(colSums(X*yy)))/length(yy)}))
    Max_lambda<- max(abs(colSums(xx_kron*y_vec)))/length(y_vec)
    #cv_lasso_lambda=seq(0, Max_lambda, length=cv_lambda_n)
    cv_lasso_lambda=round(exp(seq(log(Max_lambda*epsilon_lambda_range_min), log(Max_lambda), length.out = cv_lambda_n)), digits = 10)
  }
  #browser()
  cvm=cvsd=cvup=cvlo=0*cv_lasso_lambda
  for(ii in 1:length(cv_lasso_lambda)){
    CV_error_for_all_folds<-Cross_validation_error(Y = Y, X = X, k_fold = cv_k_fold, lasso_lambda = cv_lasso_lambda[ii],Max_EM_iter = Max_EM_iter,  replication = 1)
    cvm[ii]=mean(CV_error_for_all_folds)
    cvsd[ii]= sd(CV_error_for_all_folds)
    cvup[ii]=cvm[ii]+ 2*cvsd[ii]
    cvlo[ii]= cvm[ii]- 2*cvsd[ii]
    if(if_Print){ print(ii)}
  }

  which_min_id<-which.min(cvm)
  lambda.min=cv_lasso_lambda[which_min_id]

  cvm_min_plus_1se=  cvm[which_min_id] +  cvsd[which_min_id]
  #lambda.min.1se=max(cv_lasso_lambda[val<val_min_plus_1se])
  lambda.min.1se<-(cv_lasso_lambda[which((cvm<cvm_min_plus_1se)*(cvm>=cvm[which_min_id])==1)])
  #browser()
  #plot(cv_lasso_lambda, cvm)

  cv_lst=list(
    lambda=cv_lasso_lambda,
    cvm=cvm,
    cvsd=cvsd,
    cvup=cvup,
    cvlo=cvlo,
    lambda.min=lambda.min,
    lambda.1se=lambda.min.1se
  )
  return(cv_lst)
}

####################################
####################################


#' @export
EM_BLASSO_Dir_regression_optimizer_V1.cv<-function(Y,
                                                   X,
                                                   beta_init=NULL,
                                                   EM_tolerence=.0001,
                                                   cv_lasso_lambda=NULL, # Lasso Penalty Parameter for cross validation
                                                   Max_EM_iter=10000,
                                                   Convergence_loss=FALSE,
                                                   if_Print=TRUE,
                                                   cv_k_fold=10,
                                                   cv_lambda_n=100,
                                                   epsilon_lambda_range_min=.0001){


  #browser()
  if(is.null(cv_lasso_lambda)){
    #lso<-cv.glmnet(x=kronecker(diag(d), X), y=c(Y))
    #Max_lambda<-max(lso$lambda)
    xx_kron<-kronecker(diag(d), X);y_vec=c(Y)
    #Max_lambda<-min( apply(Y, MARGIN = 2, FUN = function(yy){ max(abs(colSums(X*yy)))/length(yy)}))
    Max_lambda<- max(abs(colSums(xx_kron*y_vec)))/length(y_vec)
    #cv_lasso_lambda=seq(0, Max_lambda, length=cv_lambda_n)
    cv_lasso_lambda=round(exp(seq(log(Max_lambda*epsilon_lambda_range_min), log(Max_lambda), length.out = cv_lambda_n)), digits = 10)
  }
  #browser()
  cvm=cvsd=cvup=cvlo=0*cv_lasso_lambda


cv_Err_Table<-Cross_validation_error_all_lambda(Y=Y, X=X,cv_lasso_lambda=cv_lasso_lambda ,k_fold = cv_k_fold,Max_EM_iter = Max_EM_iter,  if_Print=FALSE)

cvm<-apply(cv_Err_Table, MARGIN = 1, FUN = mean)
cvsd<-apply(cv_Err_Table, MARGIN = 1, FUN = sd)
cvup=cvm+ 2*cvsd
cvlo= cvm- 2*cvsd

  which_min_id<-which.min(cvm)
  lambda.min=cv_lasso_lambda[which_min_id]

  cvm_min_plus_1se=  cvm[which_min_id] +  cvsd[which_min_id]
  #lambda.min.1se=max(cv_lasso_lambda[val<val_min_plus_1se])
  lambda.min.1se<-(cv_lasso_lambda[which((cvm<cvm_min_plus_1se)*(cvm>=cvm[which_min_id])==1)])
  #browser()
  #plot(cv_lasso_lambda, cvm)

   cv_lst=list(
    lambda=cv_lasso_lambda,
    cvm=cvm,
    cvsd=cvsd,
    cvup=cvup,
    cvlo=cvlo,
    lambda.min=lambda.min,
    lambda.1se=lambda.min.1se,
    cv_Error_Table=cv_Err_Table
  )
  return(cv_lst)
}






###################################
###################################

#' @export
plot.cv.DirReg<-function(cvobj,sign.lambda=1,...){
  #cvobj=x
  xlab = expression(Log(lambda))
  #  xlab="log(Lambda)"
  if(sign.lambda<0)xlab=paste("-",xlab,sep="")
  plot.args=list(x=sign.lambda*log(cvobj$lambda),y=cvobj$cvm,ylim=range(cvobj$cvup,cvobj$cvlo),xlab=xlab,ylab="Cross Validated Error",type="n")
  new.args=list(...)
  if(length(new.args))plot.args[names(new.args)]=new.args
  do.call("plot",plot.args)
  error.bars(sign.lambda*log(cvobj$lambda),cvobj$cvup,cvobj$cvlo,width=0.01,
             col="cornflowerblue")
  #   col="antiquewhite2")

  points(sign.lambda*log(cvobj$lambda),cvobj$cvm,pch=20,
         col="purple")
  #axis(side=3,at=sign.lambda*log(cvobj$lambda),tick=FALSE,line=0)
  abline(v=sign.lambda*log(cvobj$lambda.min),lty=3,col="darkgreen")
  abline(v=sign.lambda*log(max(cvobj$lambda.1se)),lty=3, col="darkgreen")
  invisible()
}



error.bars <-
  function(x, upper, lower, width = 0.02, ...)
  {
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
  }




