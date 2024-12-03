



#Mc_Beta_burnin<-((Mc_obj_lso$MC_lst$MC$Mc_Beta[-(1:burnIN), , ]))

#Beta_est_Post_Mean<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = mean)
#Beta_est_Post_Lower_q<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = function(xx){quantile(xx, 0.25)})
#Beta_est_Post_upper_q<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = function(xx){quantile(xx, 0.75)})
#Beta_est_Post_median<-apply(Mc_Beta_burnin, MARGIN = c(2,3), FUN = function(xx){median(xx)})

location="/Users/subhadippal/Desktop/Lasso_Simulation_RBVNF/"
 fileName= paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",2,"_SimNUmber_",30,".RData")
      file_with_path= paste0(location,fileName)
      assign('Mc_obj_lso', get( load(file=file_with_path ) ))

plot_beta_p_20_d_2(Mc_obj_lso)


# library
library(ggplot2)

# create data

plot_beta_p_20_d_2<-function(Mc_obj=Mc_obj_lso, burnIN=1000){


  #location="/Users/subhadippal/Desktop/Lasso_Simulation_RBVNF/"
  #fileName= paste0("MC_SIM_BLASSO_Reg_Dir_Data_d_eq_",2,"_SimNUmber_",5,".RData")
  #file_with_path= paste0(location,fileName)
  #assign('Mc_obj_lso', get( load(file=file_with_path ) ))

  #plot_beta_p_20_d_2(Mc_obj_lso)

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
p<-p+  with_shadow(geom_boxplot(varwidth = TRUE,
                                alpha=0.85,
                                linewidth = .15,
                                outlier.size = .05, outlier.colour = "black"),
                                sigma = 2,
                                x_offset = 0,
                                y_offset = 0,
                                colour = "black" )

  p<-p+theme(legend.position="none") +scale_x_discrete(labels=eval(parse(text=my_xlab)))
p<-p+theme(panel.grid = element_line(color = "white",
                                     size = 0.15,
                                     linetype = 1))

p<-p+scale_fill_manual(values=c(replicate(4, "darkolivegreen2"), replicate(36, "lightskyblue")))
p<-p+xlab("Regression Coefficients")+ ylab(" ")+ ggtitle(" ")
return(p)
}


ffff_delete-function(lst,
         bar_col="darksalmon",
         alpha_bar=.8,
         marked_region_col="darkseagreen1",
         alpha_marked_region=.8,
         Arrow_Mark_Lambda_Min=TRUE,
         color_theme=NULL){

  # Need library(ggplot)
  # Need library(ggfx)
  if(color_theme==1){
    bar_col="darksalmon";alpha_bar=.9;  marked_region_col="darkseagreen1";  alpha_marked_region=.7;Arrow_Mark_Lambda_Min=TRUE
  }
  if(color_theme==2){
    bar_col="darkorchid";alpha_bar=.9;  marked_region_col="darkolivegreen2";  alpha_marked_region=.7;Arrow_Mark_Lambda_Min=TRUE
  }
  if(color_theme==1){
    bar_col="darksalmon";alpha_bar=.9;  marked_region_col="darkseagreen1";  alpha_marked_region=.7;Arrow_Mark_Lambda_Min=TRUE
  }
  if(color_theme==3){
    bar_col="lightskyblue";alpha_bar=.9;  marked_region_col="lightseagreen";  alpha_marked_region=.7;Arrow_Mark_Lambda_Min=TRUE
  }
  xx=lst

  #browser()

  log_lambda=log(xx$lambda)
  diff<-(log_lambda[2]-log_lambda[1])/4

  x_min_all=log_lambda-diff
  x_max_all=log_lambda+diff

  y_min_all= xx$cvlo
  y_max_all=xx$cvup
  min_y_min= min(y_min_all, na.rm = TRUE)*.9
  max_y_max=max(y_max_all, na.rm = TRUE)*1.2
  x_marked_min<-log(xx$lambda.min)-2*diff
  x_marked_max<-log(max(xx$lambda.1se))+2*diff

  ########################## Arrow Specification #####  ##############################################
  arrow1_x_start<-log(xx$lambda.min);    arrow1_x_end<-log(xx$lambda.min)
  arrow1_y_start<-max_y_max*.98;    arrow1_y_end<-min(y_min_all)*.82

  arrow2_x_start<-log(max(xx$lambda.1se));    arrow2_x_end<-log(max(xx$lambda.1se))
  arrow2_y_start<-max_y_max*.98;    arrow2_y_end<-min(y_min_all)*.82

  ########################################   ####################   #################################
  #browser()
  df<-data.frame(x=log_lambda, y= xx$cvm)
  p <- ggplot(df, aes(x = x, y = y))
  p=p+with_outer_glow(annotate("rect", xmin =x_marked_min , xmax =x_marked_max , ymin = min_y_min, ymax = max_y_max,
                               alpha = alpha_marked_region,fill =marked_region_col, col="white",linewidth=.7 ),sigma = 2, x_offset = 0, y_offset = 0, colour = "black")
  #with_outer_glow
  for(i in 1:length(y_min_all)){
    p<-p+with_outer_glow(annotate("rect", xmin = x_min_all[i], xmax = x_max_all[i], ymin = y_min_all[i], ymax = y_max_all[i],
                                  alpha = alpha_bar,fill = bar_col, col=bar_col, linewidth=.2),sigma = 2, x_offset = 0, y_offset = 0, colour = "black")
  }

  p<-p+ with_shadow(geom_line(linewidth=.2, col="white"),sigma = 2, x_offset = 0, y_offset = 0, colour = "black" )
  p<-p+with_shadow(geom_point(col="black", size=1.3), sigma = 2, x_offset = 0, y_offset = 0, colour = "white")
  #p<-p+geom_point(col="black", size=.5)
  p<-p+labs(y = "Cross-Validation Error", x = expression(Log(lambda)))

  #browser()
  if(Arrow_Mark_Lambda_Min){
    p<-p+with_shadow(annotate("segment", x =arrow1_x_start , y = arrow1_y_start*.9, xend = arrow1_x_end, yend = arrow1_y_end,
                              size = .3, linejoin = "mitre",arrow = arrow(type = "closed", length = unit(0.01, "npc")), col="white")
                     ,sigma = 2, x_offset = 0, y_offset = 0, colour = "black")
    p<-p+with_shadow(annotate("segment", x =arrow2_x_start , y = arrow2_y_start*.9, xend = arrow2_x_start, yend = arrow2_y_end,
                              size = .3, linejoin = "mitre",arrow = arrow(type = "closed", length = unit(0.01, "npc")), col="white")
                     ,sigma = 2, x_offset = 0, y_offset = 0, colour = "black")

    p<-p+with_shadow(annotate("segment", x =arrow1_x_start , y = arrow1_y_start, xend = arrow1_x_end, yend = arrow1_y_start*.9,
                              size = 4.2, linejoin = "mitre", arrow = arrow(type = "closed", length = unit(0.01, "npc")), col="darkolivegreen2")
                     ,sigma = 2, x_offset = 0, y_offset = 0, colour = "black")
    p<-p+with_shadow(annotate("segment", x =arrow2_x_start , y = arrow2_y_start, xend = arrow2_x_end, yend = arrow2_y_start*.9,
                              size = 4.2, linejoin = "mitre",arrow = arrow(type = "closed", length = unit(0.01, "npc")), col="darkolivegreen2")
                     ,sigma = 2, x_offset = 0, y_offset = 0, colour = "black")


    p<-p+annotate("text", x=arrow1_x_start,y=arrow1_y_start, label =  expression(Log(lambda["min"])), color = "black", angle = 90, hjust = 1.5, size = 2.75, fontface = "bold")
    p<-p+annotate("text", x=arrow2_x_start,y=arrow1_y_start, label =  expression(Log(lambda["1.SE"])), color = "black", angle = 90, hjust = 1.5, size = 2.75, fontface = "bold")

    #p<-p+annotate("text", x=arrow1_x_start,y=arrow1_y_start, label =  expression(Log(lambda["min"])), color = "black", angle = 0, hjust = .1, size = 2.5, fontface = "bold")
    #p<-p+annotate("text", x=arrow2_x_start,y=arrow1_y_start, label =  expression(Log(lambda["1SE"])), color = "black", angle = 0, hjust = .1, size = 2.5, fontface = "bold")

  }
  return(p)
}
