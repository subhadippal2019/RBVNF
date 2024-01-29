

#' @export
Sim_check_Plot_estimated_true_beta<-function(MCSimObj){
        library(gridExtra)
        True_beta<-c(MCSimObj$Sim_Data_lst$beta)
        Estimated_beta<-c(  apply(X = MCSimObj$MC_lst$MC$Mc_Beta, MARGIN = c(2,3), FUN = mean)  )
        Initial_EM_beta<- c(MCSimObj$MC_lst$RunDetails$beta_init)
        dd=data.frame(cbind(True_beta=True_beta, Estimated_beta=Estimated_beta, Initial_EM_beta=Initial_EM_beta ))
           #plot(True_beta, Estimated_beta)
        plt1<-ggplot(dd, aes(x=True_beta, y=Estimated_beta)) +  geom_point(size=2, shape=19)
        plt2<-ggplot(dd, aes(x=True_beta, y=Initial_EM_beta)) + geom_point(size=2, shape=19)
        plt<-grid.arrange(plt1, plt2, ncol=2)
        return(plt)
}

#' @export
Sim_Diagonistic_Plot<-function(MCSimObj,selInd=c(1,1) ){
        ii=selInd[1];jj=selInd[2]
        plt<-RBVNF::Plot_MCMC_Diag_Triplet(data = MCSimObj$MC_lst$MC$Mc_Beta[, ii,jj])
  return(plt)
}




