

# This code required the simulation Output files/MCMC samples for all the different simulation setting.
# Total number of MC from Dataset are 4*100*3=1200.  To gfenerate all the simulated replication,
# it takes time as there are total 1200 different datasets of different settings.
# Once all the MC are generated this code is used to summarize the results.

# Summary MSE
location="/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/"



all_d<-c(2,3,5, 10 )
all_DataId<-1:100
all_sampleSize<-c(750, 1000, 1500)
PredictionError_simNum_d_samSize=array(data = NA, dim = c(length(all_DataId),length(all_d), length(all_sampleSize)))


for(simData_id in all_DataId ){
  for(d_id in 1:length(all_d) ){
    for(sampleSize_id in 1:length(all_sampleSize)  ){

      sample_size_text<-ifelse(all_sampleSize[sampleSize_id]==1000, "",paste0("n_",all_sampleSize[sampleSize_id],"_") )

      fileName= paste0("MC_SIM_Reg_Dir_Data_",sample_size_text,"d_eq_",all_d[d_id],"_SimNUmber_",simData_id,".RData")
      file_with_path= paste0(location,fileName)
      assign('Mc_obj', get( load(file=file_with_path ) ))
      burnIN<-c(1:500)
      Mc_Beta_burnin<-((Mc_obj$MC_lst$MC$Mc_Beta[-burnIN, , ]))

      Beta_est_Post_Mean<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = mean)
      True_beta<-Mc_obj$Sim_Data_lst$beta

      Error_sq<- mean((Beta_est_Post_Mean-True_beta)^2)
      Error_simNum_d_samSize[simData_id,d_id, sampleSize_id ]= Error_sq

    }
  }
  print(simData_id)
}


