StoreAllMuInMatrix<-function(MC_Object){
  # browser()
  NUm_of_cluster=length(MC_Object[[1]]$curr_param)-1;
  Mu_sample=array(0,dim=c(length(MC_Object), NUm_of_cluster, 3) )
  for(MCIndex in 1:length(MC_Object)){
    for(ClusterId in 1:NUm_of_cluster)
      Mu_sample[MCIndex,ClusterId,]<-MC_Object[[MCIndex]]$curr_param[[ClusterId]]$mu
  }

  return(Mu_sample)
}


StoreAllPiInMatrix<-function(MC_Object){
  #browser()
  NUm_of_cluster=length(MC_Object[[1]]$curr_param)-1;
  Pi_sample=array(0,dim=c(length(MC_Object), NUm_of_cluster) )
  for(MCIndex in 1:length(MC_Object)){
    for(ClusterId in 1:NUm_of_cluster)
      Pi_sample[MCIndex,]<-c(MC_Object[[MCIndex]]$pi_vec)
  }

  return(Pi_sample)
}


compute_point_estimate<-function(MCMC_sample,burn_in=1){
  MCMC_len = length(MCMC_sample)

  DIC_vec = rep(0.0,MCMC_len-burn_in+1)
  MCMC_sample_1 = MCMC_sample[burn_in:MCMC_len]
  kappa_all<-StoreAllkappaInMatrix(MCMC_sample_1)
  #l1 = lapply(MCMC_sample,function(x){(F_true,x$curr_param$F,true_nclust)})
  kappa_est<-apply(kappa_all, 2, 'mean')
  mu_all=StoreAllMuInMatrix(MCMC_sample_1)
  AA=apply(mu_all, c(2, 3), 'mean')
  Norm_AA<-apply(AA, MARGIN = 1, function(x){sqrt(sum(x^2))} )
  Mu_est<-AA/Norm_AA

  PiAll<-StoreAllPiInMatrix(MCMC_sample_1)
  Pi_est<-apply(PiAll, MARGIN = 2, 'mean')


return(list(kappa_est<-kappa_est, Mu_est<-Mu_est, Pi_est<-Pi_est))

}

StoreAllkappaInMatrix<-function(MC_Object){
  # browser()
  NUm_of_cluster=length(MC_Object[[1]]$curr_param)-1;
  kappa_sample=matrix(0,length(MC_Object), NUm_of_cluster )
  for(MCIndex in 1:length(MC_Object)){
    for(ClusterId in 1:NUm_of_cluster)
      kappa_sample[MCIndex,ClusterId]<-MC_Object[[MCIndex]]$curr_param[[ClusterId]]$kappa
  }

  return(kappa_sample)
}

calculate_DIC_spiegel <-function(MCMC_sample,data,burn_in){

  MCMC_len = length(MCMC_sample)

  DIC_vec = rep(0.0,MCMC_len-burn_in+1)

  #l1 = lapply(MCMC_sample,function(x){(F_true,x$curr_param$F,true_nclust)})

  MCMC_sample_1 = MCMC_sample[burn_in:MCMC_len]
  DIC_vec = lapply(MCMC_sample_1,function(x){-2*calculate_log_like_for_data(x$pi_vec,x$curr_param,data)} )


  DIC_vec = unlist(DIC_vec)

  kappa_mu_est=compute_point_estimate(MCMC_sample)
  kappa_est<-kappa_mu_est[[1]]
  Mu_est<-kappa_mu_est[[2]]
  Pi_est<-kappa_mu_est[[3]]

  browser()
  D_bar = mean(DIC_vec)
  D_var = var(DIC_vec)
  p_d= -2*calculate_log_like_for_data(Pi_est,x$curr_param,data)

    DIC = D_bar + D_var/2

  return(DIC)
}

