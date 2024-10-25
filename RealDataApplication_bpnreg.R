install.packages("bpnreg")


######
Effect_of_var<-function(beta_EM, sel_var ){

  baseLine_dir=beta_EM[1,]/ norm(beta_EM[1,])
  beta_sel_dir<-beta_EM[sel_var, ]/norm(beta_EM[sel_var, ])


  val_return<-as.numeric(acos(t(baseLine_dir)%*%beta_sel_dir))
  return(val_return)
}
# for Motor Data
library(bpnreg)

fit.Motor <- bpnr(pred.I = Phaserad ~ 1 + Cond, data = Motor, its = 100, burn = 10, n.lag = 3)
coef_circ(fit.Motor)
coef_circ(fit.Motor, type = "categorical")


#Constructing Response and Design Matrix
Y= cbind(cos(Motor$Phaserad), sin(Motor$Phaserad))
X= cbind(replicate(dim(Y)[1], 1),exp=ifelse(Motor$Cond=='exp', 1, 0), semi.imp= ifelse(Motor$Cond=='semi.imp', 1, 0 ))



n=dim(Y)[1] # Number of the samples
p=dim(X)[2] # Number of the regression covariates
d=dim(Y)[2] # Number of directions in the direcional data


#### bbeta is a matrix of dimension p\times d
#bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
sigma_square=1
tau_square=1000


  beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)
  beta_init=beta_EM
  lst=MCMC_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =1000)

  i=3;j= 2
  Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))
  apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
  apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = sd)







  #################################################################################
  ################################## Second DataSet ( Maps data from bpnreg )  ###############################
  #################################################################################
  fit.Maps <- bpnme(pred.I = Error.rad ~ Maze + Trial.type + L.c + (1|Subject),
                    data = Maps, its = 100, burn = 1, n.lag = 1)
  coef_circ(fit.Maps)

  #Constructing Response and Design Matrix
  Y= cbind(cos(Maps$Error.rad), sin(Maps$Error.rad))
  X= cbind(replicate(dim(Y)[1], 1),as.integer(Maps$Maze)-1, as.integer(Maps$Trial.type)-1, Maps$L.c)



  n=dim(Y)[1] # NUmber of the samples
  p=dim(X)[2] # NUmber of the regression covariates
  d=dim(Y)[2] # Number of direcions in the direcional data
  #### bbeta is a matrix of dimension p\times d
  #bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
  sigma_square=1
  tau_square=1000

  beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)
  beta_init=beta_EM
  lst=MCMC_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =10000)

  i=1;j= 1
  Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))
  apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
  apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = sd)










  ########################################################
  ######## Example 2######################################
  #########################################################
  file="/Users/subhadippal/Dropbox/projects/Old projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg.csv"
  vcg_data = read.csv(file)


  X1 = vcg_data[,3:5]
  X2 = vcg_data[,6:8]

  X1_M = vcg_data[,12:14]
  X2_M = vcg_data[,15:17]

  ##### Entire data
  #output_file = "./vcg_data/vcg_output.RData"

  N = 98
  YY = array(c(0,0,0,0,0,0), c(3, 2, N))
  YY_M = YY

  for (i in 1:N){
    YY[,1,i] = t(X1[i,])
    YY[,2,i] = t(X2[i,])

    YY_M[,1,i] = t(X1_M[i,])
    YY_M[,2,i] = t(X2_M[i,])
  }

  data = YY;data_M = YY_M


  XX<-as.data.frame(model.matrix( ~ as.factor(AgeSex)-1, data=vcg_data ))

  Gender= XX[, 3]+XX[,4] # Female =1
  Age= XX[, 1]+XX[,3] # Age between 2 to 10 =1

  Y=t(YY[,1,] )
  X=as.matrix(cbind(Intercept = replicate(length(Age), 1), Gender=Gender, Age=Age, Gender_Age_Int=Gender*Age ))
  #X= cbind(replicate(dim(Y)[1], 1),as.integer(Maps$Maze)-1, as.integer(Maps$Trial.type)-1, Maps$L.c)



  n=dim(Y)[1] # NUmber of the samples
  p=dim(X)[2] # NUmber of the regression covariates
  d=dim(Y)[2] # Number of direcions in the direcional data
  #### bbeta is a matrix of dimension p\times d
  #bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
  sigma_square=1
  tau_square=1000

  beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)
  beta_init=beta_EM
  lst=MCMC_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =11000)

  i=4;j= 1
  p1=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  i=4;j= 2
  p2=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))

  library(cowplot)




pdf(file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/DirReg_WriteUpShared/fig/Plot_VCG_beta31_beta32_TripletPlot1.pdf", width = 10, height= 5)
  plot_grid(p1, p2, labels = c('A', 'B'))
  dev.off()



  ### All plots for appendix
  #Set1

  i=1;j= 1
  p1=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  j=2
  p2=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  j= 3
  p3=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  i=2;j=1
  p4=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  j=2
  p5=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  j= 3
  p6=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))



  pdf(file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/DirReg_WriteUpShared/fig/Plot_VCG_Set1_TripletPlot1.pdf", width = 10, height= 15)
  plot_grid(p1, p2,p3,p4,p5,p6, ncol = 2)
  dev.off()



  #Set2
  i=3;j= 1
  p1=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  j=2
  p2=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  j= 3
  p3=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  i=4;j=1
  p4=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  j=2
  p5=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))
  j= 3
  p6=Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i-1)][.(j)]))



  pdf(file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/DirReg_WriteUpShared/fig/Plot_VCG_Set2_TripletPlot1.pdf", width = 10, height= 15)
  plot_grid(p1, p2,p3,p4,p5,p6, ncol = 2)
  dev.off()





  Beta_est<-apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)

  Beta_sd<-apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = sd)
  Beta_est<- matrix(paste0(round(c(Beta_est),2),"(", round(c(Beta_sd),2),")"), nrow=4)
  write.csv(Beta_est, file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/VCG_Beta_est.csv")


  i=4; j=1
  xx<-(lst$MC$Mc_Beta[1001:11000, i, j])
  round(c(quantile(xx, .025), quantile(xx, .975)),2)



#############################
  ##### Example 1
  #########################

  library(foreign)
  #library(Cairo)
  germany <- read.dta("/Users/subhadippal/Dropbox/projects/Regression of Directional data/dataverse_files/germandata.dta")
  #attach(germany)
  summary(germany)
  # convert counterclockwise angles
  germany$direction.rad <- (360-germany$direction) / 360 * 2* pi
  #germany$direction.rad <- as.circular(germany$direction.rad, units = "radians",


 # Ideal direction = β0+β1 · Unemp + β2 · Outofwed + β3 · Reunification +
  #  β4 · SPD + β5 · CDU/CSU + β6 · Greens + β7 · PDS +
  #  β8 · Year + β9 · YearSQ + .

  #Constructing Response and Design Matrix
  Y= cbind(cos(germany$direction.rad), sin(germany$direction.rad))
  X_numerical= cbind(unemp=germany$unemp,
                     outofwed=germany$outofwed,
                     year=germany$year,
                     yearSQ=(germany$year)^2/100)

  X_categorical=cbind( Intercept=replicate(n = length(germany$direction.rad),1 ),
                      CduCsu=germany$cducsu,
                      Spd=germany$spd,
                      Green=germany$green,
                      Pds=germany$pds,
                      Reunification=germany$reunification)
  X_numerical =scale(X_numerical, scale = FALSE)
X=cbind(X_categorical, X_numerical)


  n=dim(Y)[1] # NUmber of the samples
  p=dim(X)[2] # NUmber of the regression covariates
  d=dim(Y)[2] # Number of direcions in the direcional data
  #### bbeta is a matrix of dimension p\times d
  #bbeta=matrix( rnorm(p*d), nrow=p, ncol=d)
  sigma_square=1
  tau_square=1000

  beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)
  beta_init=beta_EM
  lst_germany=MCMC_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =100)
  save(lst_germany, file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/workspaces/Germany_RUN_Example1_11000.rdata")


  lst_germany_LASSO=MCMC_BLASSO_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =200, lasso_lambda =.1, Sample_lasso_lambda = c(1, 100) )


  lst=lst_germany_LASSO
Calculate_DIC(lst_germany_LASSO)

  i=1;j= 1
  Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))
  Beta_est=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
  Beta_sd=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = sd)
  Beta_est1<- matrix(paste0(  round(c(Beta_est),2),"(", round(c(Beta_sd),2),")& "), nrow=10)
  Beta_est2= cbind(paste(colnames(X),"&"), Beta_est1, paste("\\\\"))
  paste0(Beta_est2, collapse="//")
  cat(Beta_est2)
  write.csv(Beta_est2, "/Users/subhadippal/Dropbox/projects/Regression of Directional data/workspaces/Germany_RUN_Example1_regCoef.csv")
 # Effect_of_var(lst$MC$Mc_Beta[1, , ],4)


  Beta_est_CI_L=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = function(x){quantile(x, .025)})
  Beta_est_CI_R=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = function(x){quantile(x, 1-.025)})
  CI_beta<- matrix(paste0("[",round(Beta_est_CI_L, 2)," , ", round(Beta_est_CI_R,2), "]"), ncol=10)



  ## CduCsu :Row= 2
  SPD_post<-apply(lst$MC$Mc_Beta, MARGIN = c(1), function(x){Effect_of_var(x,2)})
  p_1<-Plot_MCMC_Diag_Triplet(SPD_post,y_lab_text ="CduCsu")

  ## SPD :Row= 3
  SPD_post<-apply(lst$MC$Mc_Beta, MARGIN = c(1), function(x){Effect_of_var(x,3)})
  p_2<-Plot_MCMC_Diag_Triplet(SPD_post,y_lab_text ="SPD")
  ## Green :Row= 4
  Green_post<-apply(lst$MC$Mc_Beta, MARGIN = c(1), function(x){Effect_of_var(x,4)})
  p_3<-Plot_MCMC_Diag_Triplet(Green_post,y_lab_text ="Difference between `Green' and `FDP'")

  ## Green :Row= 5
  Psd_post<-apply(lst$MC$Mc_Beta, MARGIN = c(1), function(x){Effect_of_var(x,5)})
  p4<-Plot_MCMC_Diag_Triplet(Green_post,y_lab_text ="PSD")

  ###
  library(cowplot)

  pdf(file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/DirReg_WriteUpShared/fig/Plot_Ger_PartyEfffects_TripletPlot.pdf", width = 10, height= 10)
  plot_grid(p_1, p_2,p_3,p4, ncol = 2)
  dev.off()

  pdf(file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/DirReg_WriteUpShared/fig/Plot_Ger_GREEN_TripletPlot.pdf", width = 6, height= 5)
  #plot_grid(p_1, p_2,p_3,p4, ncol = 2)
  p_3
  dev.off()



  lst=lst_germany

  i=1;j= 1
  library(cowplot)
  i=1
  for(i in 1:9){
    i=10
             pdf(file=paste0("/Users/subhadippal/Dropbox/projects/Regression of Directional data/DirReg_WriteUpShared/fig/Plot_Ger_TripletPlot_Beta_",(i-1),".pdf"), width = 10, height= 5)
              #plot_grid(p_1, p_2,p_3,p4, ncol = 2)
             j=1; p1<-Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,1],y_lab_text = bquote(beta[.(i)][.(j)]))
             j=2; p2<-Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,2],y_lab_text = bquote(beta[.(i)][.(j)]))
              print(plot_grid(p1, p2, ncol = 2))
              #Sys.sleep(5)
              dev.off()
  }









  ## unemp :Row= 7
  unemp_post<-apply(lst$MC$Mc_Beta, MARGIN = c(1), function(x){Effect_of_var(x,7)})
  Plot_MCMC_Diag_Triplet(unemp_post,y_lab_text ="Effect of Unemployment Rate")

  ## Year :Row= 9
  year_post<-apply(lst$MC$Mc_Beta, MARGIN = c(1), function(x){Effect_of_var(x,9)})
  Plot_MCMC_Diag_Triplet(year_post,y_lab_text ="Year")

  Beta_est/Beta_sd





  ############ Example 1 German Political science Data plus


  sigma_square=1
  tau_square=1000

  beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)
  beta_init=beta_EM
  lst_germany=EM_BLASSO_Dir_regression_optimizer_V1(Y=Y, X=X,  beta_init = NULL,   lasso_lambda = .01)
  #save(lst_germany, file="/Users/subhadippal/Dropbox/projects/Regression of Directional data/workspaces/Germany_RUN_Example1_11000.rdata")

xx<-EM_BLASSO_Dir_regression_optimizer_V1.cv(Y=Y,
                                             X=X,
                                              beta_init = NULL,
                                             Max_EM_iter=1000,
                                             cv_k_fold = 10,
                                             cv_lambda_n = 50,
                                             epsilon_lambda_range_min = .0001
  )

  plot.cv.DirReg(xx)

  lst_germany_LASSO=MCMC_BLASSO_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =200, lasso_lambda =.1, Sample_lasso_lambda = c(1, 100) )











  Lasso_beta_EM=EM_BLASSO_Dir_regression_optimizer_V1(Y=Y, X=X,  beta_init = beta_EM,   EM_tolerence = .00001, lasso_lambda = .03, Max_EM_iter = 3000)

  Lasso_beta_EM

  beta_EM=EM_Dir_regression_optimizer_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,   EM_tolerence = .00001)
  beta_init=beta_EM



  ###### BLASSO MCMC and choosing Lambda ########
  lst_germany_LASSO=MCMC_BLASSO_Dir_regression_sampler_V1(Y=Y, X=X, prior=NULL, beta_init = NULL,  MCSamplerSize =11000, lasso_lambda = .5)
  lst=lst_germany_LASSO

  i=1;j= 1
  Plot_MCMC_Diag_Triplet(lst$MC$Mc_Beta[,i,j],y_lab_text = bquote(beta[.(i)][.(j)]))
  Beta_est=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)
  Beta_sd=apply(lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = sd)
  Beta_est1<- matrix(paste0(  round(c(Beta_est),2),"(", round(c(Beta_sd),2),")& "), nrow=10)
  Beta_est2= cbind(paste(colnames(X),"&"), Beta_est1, paste("\\\\"))
  paste0(Beta_est2, collapse="//")
  cat(Beta_est2)


  beta_EM_Lasso_LAM=EM_BLASSO_Dir_regression_optimizer_V1_EMLAMBDA(Y=Y, X=X, beta_init = NULL,   EM_tolerence = .0001, lasso_lambda = .01)


