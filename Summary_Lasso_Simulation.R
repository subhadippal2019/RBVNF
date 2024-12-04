


# Summary MSE
#location="/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/"
location="/Users/subhadippal/Desktop/Lasso_Simulation_RBVNF/"



all_d<-c(3)
all_DataId<-c(1:97)
#all_sampleSize<-c(750, 1000, 1500)
Error_sim=array(data = NA, dim = c(length(all_DataId),length(all_d)))
classification_Table= array(data = NA, dim = c(length(all_DataId),length(all_d),2 , 2))

for(simData_id in all_DataId ){
  for(d_id in 1:length(all_d) ){
    #for(sampleSize_id in 1:length(all_sampleSize)  ){

      #sample_size_text<-ifelse(all_sampleSize[sampleSize_id]==1000, "",paste0("n_",all_sampleSize[sampleSize_id],"_") )

      fileName= paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",all_d[d_id],"_SimNUmber_",simData_id,".RData")
      file_with_path= paste0(location,fileName)
      assign('Mc_obj_lso', get( load(file=file_with_path ) ))
      burnIN<-c(1:500)
      Mc_Beta_burnin<-((Mc_obj_lso$MC_lst$MC$Mc_Beta[-burnIN, , ]))

      Beta_est_Post_Mean<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = mean)
      Beta_est_Post_Lower_q<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = function(xx){quantile(xx, 0.025)})
      Beta_est_Post_upper_q<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = function(xx){quantile(xx, 0.975)})

      estimated_zero_nonzero<-c((sign(Beta_est_Post_Lower_q*Beta_est_Post_upper_q)+1)/2)
      True_beta<-(Mc_obj_lso$Sim_Data_lst$beta)
      True_beta_zero_nonzero<-c((Mc_obj_lso$Sim_Data_lst$beta!=0)*1)
      classification_Table[ simData_id,d_id  , ,]<-table(estimated_zero_nonzero, True_beta_zero_nonzero)

      Error_sq<- (mean((Beta_est_Post_Mean-True_beta)^2)
      Error_sim[simData_id,d_id ]= Error_sq

   # }
  }
  print(all_DataId[simData_id])
}


matrix(apply(X = classification_Table, MARGIN = c(2, 3, 4), FUN = sum),  nrow=2)
matrix(apply(X = classification_Table, MARGIN = c(2, 3, 4), FUN = sd),  nrow=2)

misclassification_val_d_3=1:97*0-1
for(i in 1:97){

  misclassification_val_d_3[i]=classification_Table[i, 1, 2,1 ]/ sum(classification_Table[i, 1, ,1 ])
}

summary(misclassification_val_d_2)
summary(misclassification_val_d_3)

matrix(apply(X = classification_Table, MARGIN = c(2, 3, 4), FUN = function()),  nrow=2)
matrix(apply(X = classification_Table, MARGIN = c(2, 3, 4), FUN = sd),  nrow=2)

