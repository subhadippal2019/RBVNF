
#' Plotting ACF, Density, Trace and Cumulative Average plots.

#' @examples
#'   y_lab_text=expression(paste('Cumulative Average of ',Var_name  ))
#'   y_lab_text=expression(beta^T*beta)
#'
#'   Plot_post_sample_cumAvg_Trace_Acf(data = kappa_all[,12],horizontal = FALSE, ACF = TRUE, Density=TRUE, y_lab_text = expression(paste("Posterior Samples of ", kappa)))
#' Plot_post_sample_cumAvg_Trace_Acf
#' @export

Plot_MCMC_Diag_Triplet<- function(data,
                                  lag.max = 24,
                                  ci = 0.95,
                                  large.sample.size = TRUE,
                                  horizontal = TRUE,
                                  CumType="both",
                                  ThemeBoth=2,
                                  ACF=TRUE,
                                  Density=TRUE,
                                  y_lab_text=NULL,
                                  SelBw="nrd0",...){

  #require(ggplot2)
  #require(dplyr)
  #require(cowplot)

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
      geom_col(fill = "#4f5254", width = 0.7) +theme_bw()
      #geom_area(aes(x = lag, y = qnorm((1+ci)/2)*acfstd), fill = "#B9CFE7") +
      #geom_area(aes(x = lag, y = -qnorm((1+ci)/2)*acfstd), fill = "#B9CFE7") +
      #geom_col(fill = "#4373B6", width = 0.7) +
      ##scale_x_continuous(breaks = seq(0,max(df1$lag),6)) +
      ##scale_y_continuous(name = element_blank(),
      ##limits = c(min(df1$acf,df2$pacf),1)) +
      ##ggtitle("ACF") +
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
