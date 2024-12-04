



#Mc_Beta_burnin<-((Mc_obj_lso$MC_lst$MC$Mc_Beta[-(1:burnIN), , ]))

#Beta_est_Post_Mean<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = mean)
#Beta_est_Post_Lower_q<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = function(xx){quantile(xx, 0.25)})
#Beta_est_Post_upper_q<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = function(xx){quantile(xx, 0.75)})
#Beta_est_Post_median<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = function(xx){median(xx)})

location="/Users/subhadippal/Desktop/Lasso_Simulation_RBVNF/"
data_id=18
 fileName= paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",2,"_SimNUmber_",data_id,".RData")
      file_with_path= paste0(location,fileName)
      assign('Mc_obj_lso', get( load(file=file_with_path ) ))

plt<-plot_beta_p_20_d_2(Mc_obj_lso, PlotType = 'vio', zero_marking_linetype=6, zero_marking_linewidth=.5)
ggsave(plt, file=paste0(location,paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",2,"_SimNUmber_",data_id,"_beta_plot_bw.pdf")), width = 15, height = 7)

# library
library(ggplot2)

# create data

plot_beta_p_20_d_2<-function(Mc_obj=Mc_obj_lso, burnIN=1000, PlotType="vio", zero_marking_linetype=6, zero_marking_linewidth=.8){


  #location="/Users/subhadippal/Desktop/Lasso_Simulation_RBVNF/"
  #fileName= paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",2,"_SimNUmber_",5,".RData")
  #file_with_path= paste0(location,fileName)
  #assign('Mc_obj_lso', get( load(file=file_with_path ) ))

  #plot_beta_p_20_d_2(Mc_obj_lso)

  data_summary <- function(x) {
    m <- median(x)
    ymin <- as.numeric(quantile(x,.025))#m-sd(x)
    ymax <- as.numeric(quantile(x, 0.975))#m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }

  Mc_Beta_burnin<-((Mc_obj$MC_lst$MC$Mc_Beta[-(1:burnIN), , ]))
 mc_samples= Mc_Beta_burnin
 mc_len=dim(mc_samples)[1]
 value1 <- ( matrix(mc_samples, nrow=mc_len  ))

index_rearrange=0*1:40
      for(i in 1:20){
      index_rearrange[2*i-1]= i; index_rearrange[2*i]= 20+i
      }
value= value1[,index_rearrange]
#names <- c(rep("A", mc_len) , rep("B", 5) , rep("C", 30), rep("D", 100))
beta_names<-c(t(replicate(1, paste0("beta",paste0('[paste(',row(mc_samples[1, , ]), ",',',",col(mc_samples[1, , ]), ')]')))))[index_rearrange]
beta= c(t(replicate(mc_len, paste0("beta", 111:150))))
#beta<-c(t(replicate(mc_len, paste0("beta",paste0('[paste(',row(mc_samples[1, , ]), ",',',",col(mc_samples[1, , ]), ')]')))))

data <- data.frame(Names=c(beta),Value=c(value))

# prepare a special xlab with the number of obs for each group
my_xlab <- paste0("c(", ((paste("expression(",(beta_names),")",collapse=","))), ")")

# plot


p<-ggplot(data, aes(x=Names, y=Value,fill=Names ))
if(PlotType=='vio'){
  p<-p+with_shadow(geom_hline(yintercept=0, col="white",linewidth=zero_marking_linewidth, linetype=zero_marking_linetype ),  sigma = 1,
                   x_offset = 0,
                   y_offset = 0,
                   colour = "black" )
  p<-p+  with_shadow(geom_violin(alpha=0.9,
                                  linewidth = .05, scale = 'width'),
                     sigma = 3,
                     x_offset = 0,
                     y_offset = 0,
                     colour = "black" )
  #p<-p+stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")
 p<- p + with_shadow(stat_summary(fun.data=data_summary, size=.01, col="black"),  sigma = 2,
                     x_offset = 0,
                     y_offset = 0,
                     colour = "black" )

}
if(PlotType!='vio'){
p<-p+  with_shadow(geom_boxplot(alpha=0.85,
                                linewidth = .15,
                                outlier.size = .05, outlier.colour = "black"),
                                sigma = 2,
                                x_offset = 0,
                                y_offset = 0,
                                colour = "black" )
}

  p<-p+theme(legend.position="none") +scale_x_discrete(labels=eval(parse(text=my_xlab)))
p<-p+theme(panel.grid = element_line(color = "white",
                                     size = 0.15,
                                     linetype = 1))

p<-p+scale_fill_manual(values=c(replicate(4, "white"), replicate(36, "gray")))
p<-p+xlab("Regression Coefficients")+ ylab(" ")+ ggtitle(" ")
return(p)
}





##############



location="/Users/subhadippal/Desktop/Lasso_Simulation_RBVNF/"
data_id=45
fileName= paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",3,"_SimNUmber_",data_id,".RData")
file_with_path= paste0(location,fileName)
assign('Mc_obj_lso', get( load(file=file_with_path ) ))

plt<-plot_beta_p_20_d_3(Mc_obj_lso, PlotType = 'vio', zero_marking_linetype=6, zero_marking_linewidth=.5)
ggsave(plt, file=paste0(location,paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",3,"_SimNUmber_",data_id,"_beta_plot_bw.pdf")), width = 15, height = 7)


plot_beta_p_20_d_3<-function(Mc_obj=Mc_obj_lso, burnIN=1000, PlotType="vio", zero_marking_linetype=6, zero_marking_linewidth=.8){


  #location="/Users/subhadippal/Desktop/Lasso_Simulation_RBVNF/"
  #fileName= paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",2,"_SimNUmber_",5,".RData")
  #file_with_path= paste0(location,fileName)
  #assign('Mc_obj_lso', get( load(file=file_with_path ) ))

  #plot_beta_p_20_d_2(Mc_obj_lso)

  data_summary <- function(x) {
    m <- median(x)
    ymin <- as.numeric(quantile(x,.025))#m-sd(x)
    ymax <- as.numeric(quantile(x, 0.975))#m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }

  Mc_Beta_burnin<-((Mc_obj$MC_lst$MC$Mc_Beta[-(1:burnIN), , ]))
  mc_samples= Mc_Beta_burnin
  mc_len=dim(mc_samples)[1]
  value1 <- ( matrix(mc_samples, nrow=mc_len  ))

  index_rearrange=0*1:60
  for(i in 1:20){
    index_rearrange[3*i-2]= i; index_rearrange[3*i-1]= 20+i; index_rearrange[3*i]= 40+i
  }
  value= value1[,index_rearrange]
  #names <- c(rep("A", mc_len) , rep("B", 5) , rep("C", 30), rep("D", 100))
  beta_names<-c(t(replicate(1, paste0("beta",paste0('[paste(',row(mc_samples[1, , ]), ",',',",col(mc_samples[1, , ]), ')]')))))[index_rearrange]
  beta= c(t(replicate(mc_len, paste0("beta", 111:170))))
  #beta<-c(t(replicate(mc_len, paste0("beta",paste0('[paste(',row(mc_samples[1, , ]), ",',',",col(mc_samples[1, , ]), ')]')))))
  #browser()
  data <- data.frame(Names=c(beta),Value=c(value))

  # prepare a special xlab with the number of obs for each group
  my_xlab <- paste0("c(", ((paste("expression(",(beta_names),")",collapse=","))), ")")

  # plot


  p<-ggplot(data, aes(x=Names, y=Value,fill=Names ))
  if(PlotType=='vio'){
    p<-p+with_shadow(geom_hline(yintercept=0, col="white",linewidth=zero_marking_linewidth, linetype=zero_marking_linetype ),  sigma = 1,
                     x_offset = 0,
                     y_offset = 0,
                     colour = "black" )
    p<-p+  with_shadow(geom_violin(alpha=0.9,
                                   linewidth = .05, scale = 'width'),
                       sigma = 3,
                       x_offset = 0,
                       y_offset = 0,
                       colour = "black" )
    #p<-p+stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")
    p<- p + with_shadow(stat_summary(fun.data=data_summary, size=.01, col="black"),  sigma = 2,
                        x_offset = 0,
                        y_offset = 0,
                        colour = "black" )

  }
  if(PlotType!='vio'){
    p<-p+  with_shadow(geom_boxplot(alpha=0.85,
                                    linewidth = .15,
                                    outlier.size = .05, outlier.colour = "black"),
                       sigma = 2,
                       x_offset = 0,
                       y_offset = 0,
                       colour = "black" )
  }

  p<-p+theme(legend.position="none") +scale_x_discrete(labels=eval(parse(text=my_xlab)))
  p<-p+theme(panel.grid = element_line(color = "white",
                                       size = 0.15,
                                       linetype = 1))

  p<-p+scale_fill_manual(values=c(replicate(6, "white"), replicate(54, "gray")))
  p<-p+xlab("Regression Coefficients")+ ylab(" ")+ ggtitle(" ")
  return(p)
}







#old_theme <- theme_update(
 # plot.background = element_rect(fill = "lightblue3", colour = NA),
#  panel.background = element_rect(fill = "lightblue", colour = NA),
#  axis.text = element_text(colour = "linen"),
#  axis.title = element_text(colour = "linen")
#)

#plt+theme_set(old_theme)

