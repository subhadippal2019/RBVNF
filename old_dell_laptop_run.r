
#plot(density(DTI_real_data[,4]))
#VONF_3D_GEO_DA_POSTERIOR(Y =sel_DTI_real_DATA,  MCSamplerSize = 10, )



#Sel_Q1=Q[1,,which(Aniso>.4)]
#dim(Sel_Q1)



library(Directional)
Norm<-function(x){  sqrt(sum(x^2)) }


innerProd<-function(x,y){
  return(1-sum(x*y))
}

####### hierarchical clustering for initialization of parameters
create_hclust<-function(data,K=3){
  #browser()
  N = dim(data)[2]
  #dist_mat = matrix(rep(0.0,N*N),nrow=N,ncol=N)
  # for(i in 1:N){
  #   for(j in 1:N){
  #
  #     d = 1.0-sum(diag(t(data[,i])%*%data[,j]))
  #     dist_mat[i,j] = d
  #   }
  # }
  dist_mat=dist(t(data))


  input_hclust_dist = as.dist(dist_mat)
  h = hclust(input_hclust_dist)
  id_arr = cutree(h,K)



  cluster_initial_est = list()
  for(clust_id in 1:K){
    idx = which(id_arr == clust_id)
    cluster_initial_est[[clust_id]] = vmf(t(as.matrix(data[,idx])))

    #initial_parameter_est_each_cluster()
  }
  cluster_initial_est$id_arr = id_arr
  return(cluster_initial_est)

  ## 2-sum(diag(t(L$curr_param$M[,,1])%*%init_est$M[,,1]))

}



# part of the function: cluster_param_update.r
cluster_param_update<-function(data, curr_param, cluster_assign_vec, cluster_id_vec=c(1),hyper){
  curr_param_New=curr_param;
  ### cluster_id_vec contains all the id for those clusters that needed to be udpated (e.g: c(3,1,2))
  for(cluster_id in cluster_id_vec){
    ID_cluster = which(cluster_assign_vec == cluster_id)
    #################################################
    ################################################
    ################################################
    if(length(ID_cluster) == 0){
      Valid_cluster_flag=FALSE
    }else{
      Valid_cluster_flag=TRUE
    }
    ####################################
    ############# Cluster updates ############
    ####################################
    if(!Valid_cluster_flag){curr_param_New[[cluster_id]]=curr_param[[cluster_id]] }



    if(Valid_cluster_flag){
      curr_param_tmp_OneCluster= param_update_for_One_cluster(data[,cluster_assign_vec == cluster_id], curr_param[[cluster_id]], hyper )
      curr_param_New[[cluster_id]]=curr_param_tmp_OneCluster
    } ### if N > 0
  }

  return(curr_param_New)

}



param_update_for_One_cluster<-function(sel_data, curr_param_i_th , hyper=NULL ){
  if(is.null(hyper)){hyper$beta_prior_kappa=0;hyper$alpha_prior_kappa=1}
  #datasummary
  beta_prior_kappa=hyper$beta_prior_kappa;
  alpha_prior_kappa=hyper$alpha_prior_kappa

  Y=t(sel_data)
  Y_bar=apply(Y,2,'mean');Y_SUM= apply(Y,2,'sum');n=dim(Y)[1]

  #initialization
  Z_GEO=replicate(n,-1);  kappa=curr_param_i_th$kappa;
  #McSample_kappa=replicate(MCSamplerSize,-1);McSample_mu=replicate(MCSamplerSize,0*Y_bar)


  for(i in 1:n){    Z_GEO[i]=rgeom(n=1, prob=(1-exp(-2*kappa)))  }

  #### mu: mean direction sampler #####################
  Post_mean_dir=Y_bar/Norm(Y_bar);  Post_concentration=n*kappa*Norm(Y_bar);
  mu=rvmf(n=1, mu=Post_mean_dir, k=Post_concentration)
  #### kappa: concentration sampler #####################
  #beta_Kappa_post= calculate_beta_kappa_post_GEO(Z_GEO,mu,Y_SUM)
  beta_Kappa_post=sum(2*Z_GEO+1)-(mu)%*%(Y_SUM)+beta_prior_kappa
  alpha_kappa_post=n+alpha_prior_kappa
  kappa=rgamma(n = 1,shape = alpha_kappa_post,rate = beta_Kappa_post)
  curr_param_i_th$kappa=kappa
  curr_param_i_th$mu=mu
  return(curr_param_i_th)
}




#"C:\Users\subha\Dropbox\projects\Data Augmentation for Vonmises Distribution\Codes\RealDataCode\finiteMixtureMOdel"
### parametric mixture modeling of Matrix Langevin with fixed clusters
finiteMixtureML <- function(data,nc,max_iter=2,run_id=1){
  #browser()
  #library(gtools)
  #init_run()
  set.seed(43185);N = dim(data)[[2]]; vague_prior = 1;hyper=  hyper_param_setup(vague_prior)
  ########################

  MCMC_sample = vector("list", max_iter)

  MCMC_output_file = paste0("MCMC_sample_NC_",nc,"d_DATAID_",run_id,".RData")
  init_param_output_file = sprintf("init_param_MLE_NC_",nc,"_DATAID_",run_id,".RData")
  ########### hyper parameters #######
  ### need to select empirical prior
  # hyper=NULL
  #   hyper$G = (matrix(c(0,0,0,0,0,0),ncol = 2))
  #   hyper$H = (matrix(c(0,0,0,0),ncol = 2))
  #   hyper$alpha = 1.0
  #   hyper$beta = 0.0
  #hyper$dir_alpha = rep(0.0,nc)
  #   for (j in 1:nc)
  #     hyper$dir_alpha[j] = 1.0
  #

  ####################################

  #### initialization ####

  ###curr_param = init_param_true(nc)
  init_param = create_hclust(data,nc)
  save(init_param_MLE=init_param,file=init_param_output_file)
  ########################

  curr_param = init_param
  Z = curr_param$id_arr ## Z is the latent variable here for the class association.

  ### prior for mixture weights
  pi_vec = rep(1,nc)/nc
  cluster_id_cnt = sapply(1:nc, function(x) sum(curr_param$id_arr == x))
  pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)



  ###################################################################################
  ###################################################################################
  ######################## Sampling for the Z #############################
  ###################################################################################
  ###################################################################################

  d_tmp_val = matrix(rep(0,N*nc),ncol=nc)

  for(iter in 1:max_iter){
    ######################
    if(iter%%100 == 0){
      print(paste0("MCMC main iter = ",iter))
    }
    ##### gibbs step #####
    for (i in 1:N){
      clust_assign_prob_i = rep(0,nc)
      for (cluster_id in 1:nc){
        d_tmp_val[i,cluster_id] =dVonMIsesFishar(data[,i], curr_param[[cluster_id]]$mu,curr_param[[cluster_id]]$kappa  )
        clust_assign_prob_i[cluster_id] = pi_vec[cluster_id]*d_tmp_val[i,cluster_id]
      }
      Z[i] = sample(1:nc,1,FALSE,clust_assign_prob_i)
      #print(clust_assign_prob_i)
    }
    cluster_assign_vec = Z
    #write.table(d_tmp_val,file="d_tmp_val.txt");    #tmp_table = as.matrix(table(Z))
    #cluster_id_cnt = cbind(as.integer(row.names(tmp_table)),tmp_table);    #cluster_id_cnt = rep(0,length(hyper$dir_alpha))
    #for (cluster_id in 1:nc){    #  cluster_id_cnt[cluster_id] = sum(cluster_assign_vec == cluster_id)    #} # as an alternative to previous 3 statements

    cluster_id_cnt = sapply(1:nc, function(x) sum(cluster_assign_vec == x))
    #print(cluster_id_cnt)
    pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)

    #### for all the clusters ####
    #if(iter == 1){    #  load("ML_dataset.Rdata")     #  curr_param = L$curr_param     #  cluster_assign_vec = L$clust     #}
    curr_param_new = cluster_param_update(data, curr_param, cluster_assign_vec, cluster_id_vec=1:nc,hyper)
    curr_param_new$id_arr=   cluster_assign_vec;
    ### need to update with updated cluster_assign_vec
    ##################################



    MCMC_sample[[iter]]$curr_param = curr_param_new
    MCMC_sample[[iter]]$pi_vec = pi_vec

    if(iter%%5000==0){      save("MCMC_sample", file = MCMC_output_file) ; print(paste0("Mixture modeling: MCMC iteration update ",iter))  }
  }
  save("MCMC_sample", file = MCMC_output_file)
  return(MCMC_sample)
}



setwd("/Users/subha/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/RealDataCode/")
hyper_param_setup <-function(nc,vague_prior=1){
  hyper=NULL

  if(vague_prior == 1){
    hyper$alpha_prior_kappa = 1
    hyper$beta_prior_kappa = 0
    hyper$alpha = 1.0
    hyper$beta = 0.15 ###### note that, min 5 data made a cluster 5*0.01 = 0.05,
    ###### for N==1 if S>0.99 then distribution of d1 hv a heavy tail

    hyper$dir_alpha = rep(0.0,nc)
    for (j in 1:nc)
      hyper$dir_alpha[j] = 1.0
  } ## emperical


  hyper$debug = 0

  return(hyper)

}


dVonMIsesFishar<-function(x, mu,kappa , IfLog=FALSE ){
  #  if(sum(x^2)!=1)(print("make sure that the data is on sphere"))
  # if(sum(mu^2)!=1)(print("make sure that the mu is on sphere"))
  if(kappa<0)(print("concentration, k, should be non negative. "))
  #browser()
  p=length(mu); nu=p/2-1
  log_density=(kappa*sum(x*mu))+     (nu)*log(kappa)-( log(besselI( kappa,nu, expon.scaled = TRUE))+kappa)  - (p/2)*log(2*pi)

  if(IfLog){return(log_density)}
  if(!IfLog){return(exp(log_density))}
}

library(MCMCpack)
#


# Run real data
#C:\Users\subha\Dropbox\projects\Data Augmentation for Vonmises Distribution\BVNF
load(file = '/Users/subha/Dropbox/projects/Data Augmentation for Vonmises Distribution/BVNF/DTI_processed_data.RData')

DTI_real_data=DTI_list$Dti_data

selected_index<-(DTI_real_data[,4]>=.35)
sel_DTI_real_DATA<-DTI_real_data[selected_index, 1:3 ]
rm(DTI_list);rm(DTI_real_data)




for(nCluster in 10:15){
  MC_Object=finiteMixtureML(data = t(sel_DTI_real_DATA),nc =nCluster , max_iter = 10000, run_id = 26022021)
  rm(MC_Object)
}



