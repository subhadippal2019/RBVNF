


gen_post_mu_kappa_matrix_RealData<-function(MCMC_sample){
Num_of_mcmc_sample=length(MCMC_sample)
Num_of_cluster=length(MCMC_sample[[1]]$curr_param)-1
kappa=matrix(NA, ncol=Num_of_cluster, nrow =Num_of_mcmc_sample )
MuAll=array(NA, dim=c(Num_of_mcmc_sample, Num_of_cluster, 3 ))
#browser()
for( Cluster_id in 1:Num_of_cluster){
  for( i in 1:Num_of_mcmc_sample){

    aa=MCMC_sample[[i]]$curr_param
    kappa[i, Cluster_id]=aa[[Cluster_id]]$kappa
    MuAll[i, Cluster_id,]=aa[[Cluster_id]]$mu
  }
  print(Cluster_id)
}

par=(list(kappaAll=kappa, muAll=MuAll))

return(par)
}


#' Plotting ACF, Density, Trace and Cumulative Average plots.

#' @examples
#'   y_lab_text=expression(paste('Cumulative Average of ',Var_name  ))
#'   y_lab_text=expression(beta^T*beta)
#'
#'   Plot_post_sample_cumAvg_Trace_Acf(data = kappa_all[,12],horizontal = FALSE, ACF = TRUE, Density=TRUE, y_lab_text = expression(paste("Posterior Samples of ", kappa)))
#' Plot_post_sample_cumAvg_Trace_Acf
#' @export

Plot_MCMC_Diag_Triplet<- function(data, lag.max = 24, ci = 0.95, large.sample.size = TRUE, horizontal = TRUE,CumType="both", ThemeBoth=2,ACF=TRUE, Density=TRUE,y_lab_text=NULL,SelBw="nrd0",...) {

  require(ggplot2)
  require(dplyr)
  require(cowplot)

  #browser()

  if(horizontal == TRUE) {numofrow <- 1} else {numofrow <- 2}

  if(ACF){
    list.acf <- acf(data, lag.max = lag.max, type = "correlation", plot = FALSE)
    N <- as.numeric(list.acf$n.used)
    df1 <- data.frame(lag = list.acf$lag, acf = list.acf$acf)
    df1$lag.acf <- dplyr::lag(df1$acf, default = 0)
    df1$lag.acf[2] <- 0
    df1$lag.acf.cumsum <- cumsum((df1$lag.acf)^2)
    df1$acfstd <- sqrt(1/N * (1 + 2 * df1$lag.acf.cumsum))
    df1$acfstd[1] <- 0
    #df1 <- select(df1, lag, acf, acfstd)

    list.pacf <- acf(data, lag.max = lag.max, type = "partial", plot = FALSE)
    df2 <- data.frame(lag = list.pacf$lag,pacf = list.pacf$acf)
    df2$pacfstd <- sqrt(1/N)
    #rowser()
    #if(large.sample.size == TRUE) {
    plot.acf <- ggplot(data = df1, aes( x = lag, y = acf)) +
      geom_area(aes(x = lag, y = qnorm((1+ci)/2)*acfstd), fill = "#babdbf") +
      geom_area(aes(x = lag, y = -qnorm((1+ci)/2)*acfstd), fill = "#babdbf") +
      geom_col(fill = "#4f5254", width = 0.7) +
      #geom_area(aes(x = lag, y = qnorm((1+ci)/2)*acfstd), fill = "#B9CFE7") +
      #geom_area(aes(x = lag, y = -qnorm((1+ci)/2)*acfstd), fill = "#B9CFE7") +
      #geom_col(fill = "#4373B6", width = 0.7) +
      ##scale_x_continuous(breaks = seq(0,max(df1$lag),6)) +
      ##scale_y_continuous(name = element_blank(),
      ##limits = c(min(df1$acf,df2$pacf),1)) +
      ##ggtitle("ACF") +
      theme_bw()

  }
  #df_cumAvg<-data.frame(x=1:length(data),y=data)


  if(CumType=="trace"){

    y_trace=data
    df_cumAvg<-data.frame(x=1:length(data),y=y_trace)
    cum_avg_plot<-ggplot(data=df_cumAvg, aes(x=x, y=y)) +
      geom_line() + xlab("Iteration Number") + ylab(expression(beta^T*beta)) +
      theme_bw()
  }
  if(CumType=="cumavg"){ y_trace=cumsum(data)/1:length(data)
  df_cumAvg<-data.frame(x=1:length(data),y=y_trace)
  cum_avg_plot<-ggplot(data=df_cumAvg, aes(x=x, y=y)) +
    geom_line() + xlab("Iteration Number") + ylab() +
    ylim(165, 180)+
    theme_bw()
  }

#browser()
  if(CumType=="both"){
    y_trace=cumsum(data)/1:length(data)
    df_cumAvg<-data.frame(x=1:length(data),y=y_trace, y1=data)
    if(ThemeBoth==1){

    cum_avg_plot<-ggplot(data=df_cumAvg, aes(x=x)) +
      geom_line(aes(y=y1), color="gray", size=.35, alpha=1) +
      geom_line(aes(y=y_trace), color = "white", size=4.5, alpha=.1) +
      geom_line(aes(y=y_trace), color = "white", size=3.5, alpha=.1)+
      geom_line(aes(y=y_trace), color = "white", size=2.5, alpha=.2)+
      geom_line(aes(y=y_trace), color = "white", size=1.5, alpha=.2)+
      geom_line(aes(y=y_trace), color = "white", size=.5, alpha=1)+
      xlab("Iteration Number") + ylab(y_lab_text)
    }

    if(ThemeBoth==2){
    cum_avg_plot<-ggplot(data=df_cumAvg, aes(x=x)) +
      geom_line(aes(y=y1), color="black", size=.35, alpha=1) +
      geom_line(aes(y=y_trace), color = "white", size=4.5, alpha=.1) +
      geom_line(aes(y=y_trace), color = "white", size=3.5, alpha=.1)+
      geom_line(aes(y=y_trace), color = "white", size=2.5, alpha=.2)+
      geom_line(aes(y=y_trace), color = "white", size=1.5, alpha=.2)+
      geom_line(aes(y=y_trace), color = "black", size=.25, alpha=1)+
      xlab("Iteration Number") + ylab(y_lab_text)
    }

    #+   theme_bw()
  }

  # plot.pacf <- ggplot(data = df2, aes(x = lag, y = pacf)) +
  #   geom_area(aes(x = lag, y = qnorm((1+ci)/2)*pacfstd), fill = "#B9CFE7") +
  #   geom_area(aes(x = lag, y = -qnorm((1+ci)/2)*pacfstd), fill = "#B9CFE7") +
  #   geom_col(fill = "#4373B6", width = 0.7) +
  #   scale_x_continuous(breaks = seq(0,max(df2$lag, na.rm = TRUE),6)) +
  #   scale_y_continuous(name = element_blank(),
  #                      limits = c(min(df1$acf,df2$pacf),1)) +
  #   ggtitle("PACF") +
  #   theme_bw()
  # }
  #else{
  #
  plt=cum_avg_plot
  if(ACF){
    if(Density){
      #expression_text<-parse(text =paste0("kappa[",sel_ind[ii],"]"))
      df_density<-data.frame(x=data)
      plt_density=ggplot(df_density, aes(x=x))+geom_density(color="black", fill="gray", alpha=.8, size=.25, bw=SelBw)+
                  xlab(y_lab_text)
     # plt=cowplot::plot_grid(plot.acf,plt_density,  cum_avg_plot,rel_widths = c(1, 1, 2),nrow = 2)

      first_row = plot_grid(plot.acf,   plt_density, nrow=1)
      second_row = plot_grid(cum_avg_plot)
      plt = plot_grid(first_row, second_row, labels=c('', ''), ncol=1)
    }
    if(!Density){
    plt=cowplot::plot_grid(plot.acf, cum_avg_plot, nrow = numofrow)
    }
  }



  #plt
  return(plt)
}




Plot_MCMC_Diag_Triplet_multiple<-function(MCMC, NUM=4){
  kappa_mu<-gen_post_mu_kappa_matrix_RealData(MCMC_sample)
  kappaAll<-kappa_mu$kappaAll
   Num_of_cluster<-ncol(kappaAll)

   sel_ind<-sample(1:Num_of_cluster, NUM, replace = FALSE)

   plot_list<-list()
   for( ii in 1:NUM ){
     expression_text<-parse(text =paste0("kappa[",sel_ind[ii],"]"))
     #expression_text<-parse(text =paste0("kappa 'for Culster ' ",sel_ind[ii],""))
     plot_list[[ii]]<- Plot_MCMC_Diag_Triplet(kappaAll[,sel_ind[ii]], y_lab_text = expression_text)

   }
   first_row=list()
    if( NUM%%2==0){
      for( ii in 1:(NUM/2) ){
        first_row[[ii]] = plot_grid(plot_list[[2*ii-1]],  plot_list[[2*ii]], nrow=1)
      }
      code_str<-paste0("plt = plot_grid(", paste0("first_row[[",1:(NUM/2),"]],", collapse = " ")," ncol=1)")
      eval(parse(text=code_str))
      return(plt)
    }
   #browser()
   if( NUM%%2){
     if(!NUM%%3){
     for( ii in 1:(NUM/3) ){
       first_row[[ii]] = plot_grid(plot_list[[3*ii-2]], plot_list[[3*ii-1]],   plot_list[[3*ii]], nrow=1)
     }
     code_str<-paste0("plt = plot_grid(", paste0("first_row[[",1:(NUM/3),"]],", collapse = " ")," ncol=1)")
     eval(parse(text=code_str))
     return(plt)
     }
   }




}

#file_name="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/Diag_triplet.pdf"

#ggsave(plt, file=file_name, width = 10 , height = 10)



Plot_MCMC_Diag_Triplet_multiple_selected<-function(MCMC, selected_ind=c(1)){
  kappa_mu<-gen_post_mu_kappa_matrix_RealData(MCMC_sample)
  kappaAll<-kappa_mu$kappaAll
  Num_of_cluster<-ncol(kappaAll)

  #sel_ind<-sample(1:Num_of_cluster, NUM, replace = FALSE)
  sel_ind<-selected_ind
  NUM=length(sel_ind)
  plot_list<-list()
  for( ii in 1:NUM ){
    expression_text<-parse(text =paste0("kappa[",sel_ind[ii],"]"))
    #expression_text<-parse(text =paste0("kappa 'for Culster ' ",sel_ind[ii],""))
    plot_list[[ii]]<- Plot_MCMC_Diag_Triplet(kappaAll[,sel_ind[ii]], y_lab_text = expression_text)

  }
  first_row=list()
  if( NUM%%2==0){
    for( ii in 1:(NUM/2) ){
      first_row[[ii]] = plot_grid(plot_list[[2*ii-1]],  plot_list[[2*ii]], nrow=1)
    }
    code_str<-paste0("plt = plot_grid(", paste0("first_row[[",1:(NUM/2),"]],", collapse = " ")," ncol=1)")
    eval(parse(text=code_str))
    return(plt)
  }
  #browser()
  if( NUM%%2){
    if(!NUM%%3){
      for( ii in 1:(NUM/3) ){
        first_row[[ii]] = plot_grid(plot_list[[3*ii-2]], plot_list[[3*ii-1]],   plot_list[[3*ii]], nrow=1)
      }
      code_str<-paste0("plt = plot_grid(", paste0("first_row[[",1:(NUM/3),"]],", collapse = " ")," ncol=1)")
      eval(parse(text=code_str))
      return(plt)
    }
  }




}



plt1<-Plot_MCMC_Diag_Triplet_multiple_selected(MCMC = MCMC_sample, selected_ind = sel_kappa_ind)

ggsave(plt1, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MCMC_read_data_cluster_9_TripletPlot.pdf", width=10, height=10)


load(file = "/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/RealDataCode/MCMC_sample_NC_13d_DATAID_26022021.RData")



 point_est<-compute_point_estimate(MCMC_sample, burn_in = 1000)
sel_kappa_ind<-which(rank((point_est[[3]]))>5)
  Mat=point_est[[2]][sel_kappa_ind,]
cormat=round(1-abs(Mat%*%t(Mat)), 5)
# Reorder the correlation matrix


write.csv(point_est[[1]], file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MCMC_real_data_cluster_9_point_est_kappa.csv")


write.csv(point_est[[2]], file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MCMC_real_data_cluster_9_point_est_mu.csv")



write.csv(point_est[[3]], file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MCMC_real_data_cluster_9_point_est_pi.csv")

#cbind(pi=point_est[[3]], kappa=point_est[[1]], mu=point_est[[2]])
est<-round(cbind(pi=point_est[[3]], kappa=point_est[[1]], mu=point_est[[2]]), 3)
write.csv(est, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MCMC_real_data_cluster_9_point_est.csv")





get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(0,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)



load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_0.01.RData")

ii=13
MCMC_sample=AllList_3D[[ii]]$lst_GEO_3D$McSample_kappa
plot(MCMC_sample, type="l")

plt1<-Plot_MCMC_Diag_Triplet(data = MCMC_sample, SelBw=.012, y_lab_text = expression(kappa))
plt1

ggsave(plt1, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/GEO_DA_Kappa_01_3d_1.pdf", width=6, height=6)

ii=7
MCMC_sample=AllList_3D[[ii]]$lst_MPG_3D$McSample_kappa
plt1<-Plot_MCMC_Diag_Triplet(data = MCMC_sample, SelBw=.02, y_lab_text = expression(kappa))
plt1



ggsave(plt1, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MPG_DA_Kappa_01_3d_1.pdf", width=6, height=6)




#### comparison plots
alpha_seq=c(5, 10, 15, 20, 25, 30)
for (ii in 1: length(alpha_seq)){
alpha=alpha_seq[ii]; nu=10
plt<-Plot_MPG_GIG_CDF( alpha = alpha, nu = nu)
ggsave(plt, file=paste0("/Users/subhadippal/Box/DA for Directional data/WriteUp/6.6/figures/CDF plots/MPG(white)_GIG(black)_CDF_nu_", nu,"_alpha_", alpha,".pdf"), width=5, height=3)
}





#### comparison plots
alpha_seq=c(1, 5, 10, 15, 20, 25, 30)
for (ii in 1: length(alpha_seq)){
  alpha=alpha_seq[ii]; nu=40
  plt<-Plot_MPG_GIG_CDF_alt( alpha = alpha, nu = nu)
  ggsave(plt, file=paste0("/Users/subhadippal/Box/DA for Directional data/WriteUp/6.6/figures/CDF plots/MPG(white)_GIG_opt(black)_GIG(gray_dash)_CDF_nu_", nu,"_alpha_", alpha,".pdf"), width=5, height=3)
}


