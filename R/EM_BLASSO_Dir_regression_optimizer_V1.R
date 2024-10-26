

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
  names_X=colnames(X)
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
    #I_KRON_NEG_SQRT_D_ET_Y_1= (1/sqrt(ET))*Y/2

    #mStandard_Lasso_fit = glmnet(x =SQRT_D_ET_X , y = I_KRON_NEG_SQRT_D_ET_Y_1, family = "mgaussian", lambda=lasso_lambda, alpha =1 )
    #mbetaMat_opt<- as.matrix(do.call(cbind, mStandard_Lasso_fit$beta))

  #if(any(is.nan(I_KRON_SQRT_D_ET_X))){browser()}
  # if(any(is.nan(I_KRON_NEG_SQRT_D_ET_Y))){browser()}
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
    if(sum(betaMat_opt^2)<0.0000001){ EMiter=Max_EM_iter}
    EMiter=EMiter+1
    if(if_Print){print(EMiter)}
  }
  #end While ########################################################################################
  rownames(betaMat_opt)=names_X; colnames(betaMat_opt)= paste0("Y_", 1:dim(Y)[2])
  if(Convergence_loss){ return_lst<-list(betaMat_opt=betaMat_opt,conv=diff_vec) }
  if(!Convergence_loss){ return_lst=betaMat_opt }

  return( return_lst)
}



####################################################################################################################
####################################################################################################################
##### Codes for the Cross Validation ######


Predict_Y <- function(Beta_hat, X_new){
  Y_unnorm  <- X_new%*% (Beta_hat)
  Y_hat     <- t(apply(Y_unnorm, MARGIN =1,  function(y){return(y/(sqrt(sum(y^2)))+1e-20)} ))
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
    Y_train= Y[-sel_rows, ]; X_train= X[-sel_rows, ]
    Y_test= matrix(Y[sel_rows, ],ncol = col_Num); X_test= matrix(X[sel_rows, ], ncol=X_col_Num)
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

      Y_hat= Predict_Y(Beta_hat = Beta_est, X_new = X_test )
      Err[lambda_id, i]= Compute_Error(Y_true = Y_test, Y_hat = Y_hat  )
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
                Y_train= Y[-sel_rows, ]; X_train= X[-sel_rows, ]
                Y_test= matrix(Y[sel_rows, ],ncol = col_Num); X_test= matrix(X[sel_rows, ], ncol=X_col_Num)
                print(paste0("dimension of X in Training Set=", dim(X_train)))
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
                                                   epsilon_lambda_range_min=.0001,
                                                   lambda_Range_Type=2){


  #browser()

  if(is.null(cv_lasso_lambda)){
        if(lambda_Range_Type!=1){
              lso<-cv.glmnet(x=kronecker(diag(d), X), y=c(Y))
              #cv_lasso_lambda<-(lso$lambda)
              Max_lambda=max(lso$lambda)
        }
         if(lambda_Range_Type==1){
                    xx_kron<-kronecker(diag(d), X);y_vec=c(Y)
                    #Max_lambda<-mean( apply(Y, MARGIN = 2, FUN = function(yy){ max(abs(colSums(X*yy)))/length(yy)}))
                   Max_lambda<- 2*max(abs(colSums(xx_kron*y_vec)))/length(y_vec)
                   #cv_lasso_lambda=seq(0, Max_lambda, length=cv_lambda_n)
         }
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
              lambda.1se=max(lambda.min.1se),
              cv_Error_Table=cv_Err_Table
              )

  return(cv_lst)
}






###################################
###################################

#' @export
plot.cv.Dir_Lasso_Reg<-function(cvobj,sign.lambda=1,...){
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




########## Advanced veriosn of the plot need ggplot

#' @export
plot.cv.Dir_Lasso_Reg_gg<-function(lst,
                                   bar_col="darksalmon",
                                   alpha_bar=.8,
                                   marked_region_col="darkseagreen1",
                                   alpha_marked_region=.8,
                                   Arrow_Mark_Lambda_Min=TRUE,
                                   color_theme=NULL){

  # Need library(ggplot)
  # Need library(ggfx)
  if(color_theme==1){
    bar_col="darksalmon";alpha_bar=.9;  marked_region_col="darkseagreen1";  alpha_marked_region=.7;Arrow_Mark_Lambda_Min=TRUE
  }
  if(color_theme==2){
    bar_col="darkorchid";alpha_bar=.9;  marked_region_col="darkolivegreen2";  alpha_marked_region=.7;Arrow_Mark_Lambda_Min=TRUE
  }
  if(color_theme==1){
    bar_col="darksalmon";alpha_bar=.9;  marked_region_col="darkseagreen1";  alpha_marked_region=.7;Arrow_Mark_Lambda_Min=TRUE
  }
  if(color_theme==3){
    bar_col="lightskyblue";alpha_bar=.9;  marked_region_col="lightseagreen";  alpha_marked_region=.7;Arrow_Mark_Lambda_Min=TRUE
  }
  xx=lst

  #browser()

  log_lambda=log(xx$lambda)
  diff<-(log_lambda[2]-log_lambda[1])/4

  x_min_all=log_lambda-diff
  x_max_all=log_lambda+diff

  y_min_all= xx$cvlo
  y_max_all=xx$cvup
  min_y_min= min(y_min_all, na.rm = TRUE)*.9
  max_y_max=max(y_max_all, na.rm = TRUE)*1.2
  x_marked_min<-log(xx$lambda.min)-2*diff
  x_marked_max<-log(max(xx$lambda.1se))+2*diff

  ########################## Arrow Specification #####  ##############################################
              arrow1_x_start<-log(xx$lambda.min);    arrow1_x_end<-log(xx$lambda.min)
              arrow1_y_start<-max_y_max*.98;    arrow1_y_end<-min(y_min_all)*.82

              arrow2_x_start<-log(max(xx$lambda.1se));    arrow2_x_end<-log(max(xx$lambda.1se))
              arrow2_y_start<-max_y_max*.98;    arrow2_y_end<-min(y_min_all)*.82

  ########################################   ####################   #################################
  #browser()
  df<-data.frame(x=log_lambda, y= xx$cvm)
  p <- ggplot(df, aes(x = x, y = y))
  p=p+with_outer_glow(annotate("rect", xmin =x_marked_min , xmax =x_marked_max , ymin = min_y_min, ymax = max_y_max,
                               alpha = alpha_marked_region,fill =marked_region_col, col="white",linewidth=.7 ),sigma = 2, x_offset = 0, y_offset = 0, colour = "black")
  #with_outer_glow
  for(i in 1:length(y_min_all)){
           p<-p+with_outer_glow(annotate("rect", xmin = x_min_all[i], xmax = x_max_all[i], ymin = y_min_all[i], ymax = y_max_all[i],
                                  alpha = alpha_bar,fill = bar_col, col=bar_col, linewidth=.2),sigma = 2, x_offset = 0, y_offset = 0, colour = "black")
  }

  p<-p+ with_shadow(geom_line(linewidth=.2, col="white"),sigma = 2, x_offset = 0, y_offset = 0, colour = "black" )
  p<-p+with_shadow(geom_point(col="black", size=1.3), sigma = 2, x_offset = 0, y_offset = 0, colour = "white")
  #p<-p+geom_point(col="black", size=.5)
  p<-p+labs(y = "Cross-Validation Error", x = expression(Log(lambda)))

      #browser()
      if(Arrow_Mark_Lambda_Min){
                    p<-p+with_shadow(annotate("segment", x =arrow1_x_start , y = arrow1_y_start*.9, xend = arrow1_x_end, yend = arrow1_y_end,
                                  size = .3, linejoin = "mitre",arrow = arrow(type = "closed", length = unit(0.01, "npc")), col="white")
                                  ,sigma = 2, x_offset = 0, y_offset = 0, colour = "black")
                    p<-p+with_shadow(annotate("segment", x =arrow2_x_start , y = arrow2_y_start*.9, xend = arrow2_x_start, yend = arrow2_y_end,
                                  size = .3, linejoin = "mitre",arrow = arrow(type = "closed", length = unit(0.01, "npc")), col="white")
                                  ,sigma = 2, x_offset = 0, y_offset = 0, colour = "black")

                    p<-p+with_shadow(annotate("segment", x =arrow1_x_start , y = arrow1_y_start, xend = arrow1_x_end, yend = arrow1_y_start*.9,
                                  size = 4.2, linejoin = "mitre", arrow = arrow(type = "closed", length = unit(0.01, "npc")), col="darkolivegreen2")
                                  ,sigma = 2, x_offset = 0, y_offset = 0, colour = "black")
                    p<-p+with_shadow(annotate("segment", x =arrow2_x_start , y = arrow2_y_start, xend = arrow2_x_end, yend = arrow2_y_start*.9,
                                  size = 4.2, linejoin = "mitre",arrow = arrow(type = "closed", length = unit(0.01, "npc")), col="darkolivegreen2")
                                  ,sigma = 2, x_offset = 0, y_offset = 0, colour = "black")


                    p<-p+annotate("text", x=arrow1_x_start,y=arrow1_y_start, label =  expression(Log(lambda["min"])), color = "black", angle = 90, hjust = 1.5, size = 2.75, fontface = "bold")
                    p<-p+annotate("text", x=arrow2_x_start,y=arrow1_y_start, label =  expression(Log(lambda["1.SE"])), color = "black", angle = 90, hjust = 1.5, size = 2.75, fontface = "bold")

                    #p<-p+annotate("text", x=arrow1_x_start,y=arrow1_y_start, label =  expression(Log(lambda["min"])), color = "black", angle = 0, hjust = .1, size = 2.5, fontface = "bold")
                    #p<-p+annotate("text", x=arrow2_x_start,y=arrow1_y_start, label =  expression(Log(lambda["1SE"])), color = "black", angle = 0, hjust = .1, size = 2.5, fontface = "bold")

      }
  return(p)
}


