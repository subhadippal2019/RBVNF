
# n=750   ; p=30       ; d=3
# n=750   ; p=30       ; d=2
# n=750   ; p=30       ; d=5
# n=200   ; p=30       ; d=3


n=1000    # Number of the samples
p=30      # Number of the regression covariates
d=5      # Number of directions in the directional data

#Test Lasso
data_lst  =   Data_generator_vnf_reg(n=2000, p=p, d=d, concentration_factor = 1, beta_factor = 5)
data_lst  <-  Data_generator_vnf_reg_sparse(n=200, p=p, d=d,SetUp = 2, NumOfNonZeroBeta=c(4, 1, 10))



Y = data_lst$Y;X=data_lst$X;
#cv.glmnet(x = X, y = Y[,1])

beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)


beta_EM_Lasso=EM_BLASSO_Dir_regression_optimizer_V1(Y=Y, X=X, beta_init = NULL,   EM_tolerence = .00001, lasso_lambda = .002)

beta_EM_Lasso_LAM=EM_BLASSO_Dir_regression_optimizer_V1_EMLAMBDA(Y=Y, X=X, beta_init = NULL,   EM_tolerence = .00001, lasso_lambda = 0.025, scale_factor = 1) # EM recomendation is scale_factor=2


###################
#Cross_validation_error(Y = Y, X = X, k_fold = 5, lasso_lambda = .045, Max_EM_iter = 2000)

#lso<-cv.glmnet(x=kronecker(diag(d), X), y=c(Y))
#lso$lambda.min

#cv_set=seq(0.001, 0.38, length=100)
#val=0*cv_set
#for(ii in 1:length(cv_set)){
#val[ii]=mean(Cross_validation_error(Y = Y, X = X, k_fold = 10, lasso_lambda = cv_set[ii],Max_EM_iter = 2000,  replication = 1))*1000
#print(ii)
#}



#plot(cv_set, val, type="l")
#cv_set[which.min(val)]


system.time(xx<-EM_BLASSO_Dir_regression_optimizer_V1.cv(Y=data_lst$Y,
                                                         X=data_lst$X,
                                                         beta_init = NULL,
                                                         Max_EM_iter=1000,
                                                         cv_k_fold = 10,
                                                         cv_lambda_n = 10,
                                                         epsilon_lambda_range_min = .0001,
                                                         lambda_Range_Type = 2
                                                         )
                                                      )

plot.cv.Dir_Lasso_Reg(xx)
plt_cv_lambda<-plot.cv.Dir_Lasso_Reg_gg(xx, color_theme = 2)

lst_em_lasso_object<-list(crossValidation_xx=xx,datadata=data_lst,plt_cv_lambda=plt_cv_lambda, n=200, p=30, Num_nonZero=4, d=5 )
save(lst_em_lasso_object, file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/workspaces/lst_em_lasso_object_n=200_p_30_d_5_nonzero_4.RData")


ggsave(lst_em_lasso_object$plt_cv_lambda, width = 15, height = 7,
       file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/DirReg_WriteUpShared/fig/lst_em_lasso_object_n=200_p_30_d_5_nonzero_4.pdf")



load("/Users/subhadippal/Dropbox/projects/Regression of Directional data/workspaces/lst_em_lasso_object_n=200_p_30_nonzero_4.RData")
beta_EM_Lasso=EM_BLASSO_Dir_regression_optimizer_V1(Y=data_lst$Y, X=data_lst$X, beta_init = NULL, lasso_lambda = max(xx$lambda.1se),   EM_tolerence = .00001)

beta_EM_Lasso1=EM_BLASSO_Dir_regression_optimizer_V1(Y=data_lst$Y, X=data_lst$X, beta_init = NULL, lasso_lambda = max(xx$lambda.min),   EM_tolerence = .00001)






lst_BLASSO_Beta_MCMC=MCMC_BLASSO_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =1000,
                                                           lasso_lambda_spec = list(
                                                             Type="SAMPLE",
                                                             lasso_lambda=max(xx$lambda.1se),
                                                             hyper_lambda_selector= list(Lambda_lower=xx$lambda.min, Lambda_upper=xx$lambda.1se,info_weight=10 ) ))
lst=lst_BLASSO_Beta_MCMC







lst_BLASSO_Beta_MCMC_sample_Lambda=MCMC_BLASSO_Dir_regression_sampler_V1(Y=Y,
                                                           X=X,
                                                           prior=NULL,
                                                           beta_init = NULL,
                                                           MCSamplerSize =1000,
                                                           lasso_lambda_spec = list(
                                                                        Type="SAMPLE",
                                                                         lasso_lambda=max(xx$lambda.1se),
                                                                         hyper_lambda_selector= list(Lambda_lower=xx$lambda.min, Lambda_upper=xx$lambda.1se,info_weight=10 ) ))


lst=lst_BLASSO_Beta_MCMC_sample_Lambda
Beta_est=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
Plot_MCMC_Diag_Triplet(lst$MC$lasso_lambda_all,y_lab_text = bquote(lambda))
abs(Beta_est1)-abs(Beta_est)

i=10;j= 2
Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))
Beta_est=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
Beta_sd=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = sd)
Beta_est1<- matrix(paste0(  round(c(Beta_est),2),"(", round(c(Beta_sd),2),")& "), nrow=10)
Beta_est2= cbind(paste(colnames(X),"&"), Beta_est1, paste("\\\\"))
paste0(Beta_est2, collapse="//")
cat(Beta_est2)









#### New more efficient function
Cross_validation_error_all_lambda_BCUp<-function(Y=Y, X=X, cv_lasso_lambda=c(0.01, 0.011, 0.012), k_fold=10,  if_Print=FALSE, Max_EM_iter=10000, strategy_1=TRUE){
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

      Y_hat= Predict_Y(Beta_hat =Beta_est, X_new = X_test )
      Err[lambda_id, i]= Compute_Error(Y_true =Y_test,Y_hat=Y_hat  )
      }


  }
  return((Err))
}



#system.time( Cross_validation_error_all_lambda(Y=Y, X=X,cv_lasso_lambda=c(0.01, 0.0100) ,k_fold = 2, if_Print=FALSE) )








