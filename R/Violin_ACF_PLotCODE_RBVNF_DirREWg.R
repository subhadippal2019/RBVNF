load("/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/SimulationSummary_ACF_Lag_200.Rdata")




Lag=30
      file="/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/figures/Vioplot_ACF_DirReg_n_eq_750_lag_"
      file=paste0(file, Lag,".pdf")
  pdf(file=file, width = 6.5, height =3.5 )
  plot_Acf_violiin_RBVNF(ACF_simNum_d_samSize[,,,1],MAXLAG = Lag, size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()




#file="/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/figures/Vioplot_ACF_DirReg_n_eq_1000_lag25.pdf"
  file="/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/figures/Vioplot_ACF_DirReg_n_eq_1000_lag_"
  file=paste0(file, Lag,".pdf")
#load("/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/SimulationSummary_ACF.Rdata")

pdf(file=file, width = 6.5, height =3.5 )
  plot_Acf_violiin_RBVNF(ACF_simNum_d_samSize[,,,2],MAXLAG = Lag, size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()




#file="/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/figures/Vioplot_ACF_DirReg_n_eq_1500_lag25.pdf"
file="/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/figures/Vioplot_ACF_DirReg_n_eq_1500_lag_"
file=paste0(file, Lag,".pdf")
#load("/Users/subhadippal/Documents/Google Drive/MyRPackages/MC_SIM_Reg_Dir_Data/SimulationSummary_ACF.Rdata")

pdf(file=file, width = 6.5, height =3.5 )
plot_Acf_violiin_RBVNF(ACF_simNum_d_samSize[,,,3],MAXLAG = Lag, size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()



sel_dim=3# n=1500 ; 2== n=1000, 1==n=750
acf_mean<-(apply(ACF_simNum_d_samSize[,,,sel_dim], MARGIN = c(1,3), FUN = mean))[2:11,]
acf_sd<-(apply(ACF_simNum_d_samSize[,,,sel_dim], MARGIN = c(1,3), FUN = sd))[2:11, ]

aa<-matrix(paste0(round(acf_mean,3)," (", round(acf_sd, 3),") "), ncol=4)




for(ii in 1:dim(aa)[1]){
    cat(paste("Lag=", ii," & ", paste(aa[ii,], collapse=" & "), "\\\\ \\hline \n"))
}


