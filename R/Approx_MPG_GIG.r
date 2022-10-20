f_mpg<-function(x){
  if(length(x)>1){
    return(apply(matrix(x, ncol=1), MARGIN = 1, function(y){log_f_nu(t=y, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = FALSE)}))
  }
  return(log_f_nu(t=x, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = FALSE))

}


f_gig<-function(x){
  val=dgig(x =x, lambda =-(nu+1), chi = .5, psi = 2*alpha^2 )
  return(val)
}




log_f_mpg<-function(x, alpha,nu,  j_0, J_1){
  if(length(x)>1){
    return(apply(matrix(x, ncol=1), MARGIN = 1, function(y){log_f_nu(t=y, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = TRUE)}))
  }
  return(log_f_nu(t=x, a=alpha, nu=nu,j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = TRUE))

}


log_f_gig<-function(x, alpha, nu){
  val=dgig(x =x, lambda =-(nu+1), chi = .5, psi = 2*alpha^2 , log = TRUE)
  return(val)
}




f_diff<-function(x, alpha, j_0, J_1 ){
  #ff1=f_mpg(x,alpha )
  #val=(ff1-f_gig(x,alpha))
  return(exp(log_f_mpg(x, alpha,j_0=j_0, J_1=J_1))-exp(log_f_gig(x=x, alpha=alpha)))
}





###################

#' Plot_Approximation_MPG_GIG_Relative_diff
#'
#' @examples
#' x=seq(from = .01,to=.15,by = .0001)
#' alpha=seq(from=10, to = 30, by = 1)
#' pp=Plot_Approximation_MPG_GIG_Relative_diff(x=seq(from = .01,to=.15,by = .001),alpha= seq(from=10, to = 20, by = 1))
#' pp$p
#' pp$p1
#' @export
Plot_Approximation_MPG_GIG_Relative_diff<-function(x, alpha, nu=0, bgCol="#686868"){

  #' x=seq(from = .01,to=.15,by = .0001)
  #' alpha=seq(from=10, to = 30, by = 1)
  #' pp=Plot_Approximation_MPG_GIG_Relative_diff(x=seq(from = .01,to=.15,by = .001),alpha= seq(from=10, to = 20, by = 1))
  #' pp$p
  #' pp$p1
#nu=0#### change nu for different results
j_0=bessel_zero_Jnu(nu = nu, s = 1:1000)
J_1=besselJ(j_0,nu=1+nu);


Alpha_all=alpha %x% replicate(length(x),1)
x_all=replicate(length(alpha),1)%x% x

data_mat=cbind(Alpha_all,x_all)
#diff_den1=apply(data_mat, MARGIN = 1, FUN = function(x_vec){ f_diff(x_vec[2], x_vec[1]) })
den1=apply(data_mat, MARGIN = 1, FUN = function(x_vec){ exp(log_f_mpg(x = x_vec[2],alpha =  x_vec[1],nu=nu, j_0 = j_0, J_1 = J_1)) })
den2=apply(data_mat, MARGIN = 1, FUN = function(x_vec){ exp(log_f_gig(x_vec[2], x_vec[1], nu) ) })

df=data.frame(x=x_all, alpha=Alpha_all, den1=den1, den2=den2, diff_den=den1-den2)

df$scaled_diff=df$diff_den*0
for(index in alpha){
  mx=max(df$den1[df$alpha==index])
  df$scaled_diff[df$alpha==index]=df$diff_den[df$alpha==index]/mx
}


p=ggplot(data=df, aes(x=x, y=scaled_diff, group=alpha)) +
  #geom_line(aes(color=alpha), size=.5)+scale_color_gradient(low="skyblue", high="#9999CC")+
  #geom_line(aes(color=alpha), size=.5)+scale_color_gradient(low="white", high="black")+
  geom_line(aes(color=alpha), size=.2)+scale_color_gradient(low="black", high="white")+
  ylim(-.012,.012)+
  xlab("Support of the Distributions")+
  ylab("Relative Difference Between the Densities ") +
  labs(color=expression(alpha))

# p+theme(
#   panel.background = element_rect(fill = bgCol,
#                                   colour = bgCol,
#                                   size = 0.5, linetype = "solid"),
#   panel.grid.major = element_line(size = 0.002, linetype = 'solid',
#                                   colour = "white"),
#   panel.grid.minor = element_line(size = 0.001, linetype = 'solid',
#                                   colour = "lightgray")
#   )
p+theme(
  panel.background = element_rect(fill =bgCol ,
                                  colour = bgCol,
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                  colour = "lightgray"),
  panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                  colour = "lightgray"))

#dev.off()
#p


p1=ggplot(data=df, aes(x=x, y=diff_den, group=alpha)) +
  geom_line(aes(color=alpha))+scale_color_gradient(low="skyblue", high="#9999CC")+
  ylim(-1,1)+
  xlab("Support of the Distributions")+
  ylab("Difference Between the Densities") +
  labs(color=expression(alpha))

#p1
return(list(p=p,p1=p1))

}











#' @export
plot_GIG_MPG_densities<-function(alpha=5, nu=0 , UpperLim=NULL, LowerLim=NULL){
  f_mpg<-function(x){
    if(length(x)>1){
      return(apply(matrix(x, ncol=1), MARGIN = 1, function(y){log_f_nu(t=y, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = FALSE)}))
    }
    return(log_f_nu(t=x, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = FALSE))

  }


  f_gig<-function(x){
    val=dgig(x =x, lambda =-(nu+1), chi = .5, psi = 2*alpha^2 )
    return(val)
  }
  #nu=.5
  #alpha=5;
  j_0=bessel_zero_Jnu(nu = nu, s = 1:1000)
  J_1=besselJ(j_0,nu=1+nu);
  if(is.null(UpperLim)){
    UpperLim=.05+.75/alpha
  }
  if(is.null(LowerLim)){LowerLim=.006}
  #UpperLim=.25*(nu==0)+ .15*(nu==.5)+ .1*(nu>.5)*(nu<4)+ .06*(nu>4)
  library(ggplot2)
  nu_name=nu
  if(nu==.5){nu_name="half"}
  #browser()
  print(alpha)


  p=ggplot(data.frame(x = rnorm(100)), aes(x)) +
    stat_function(fun = function(x){f_mpg(x = x )}, colour = "white", size=.3)+
    stat_function(fun = function(x){f_gig(x = x)}, colour = "black", size=.2)+
    xlim(LowerLim, UpperLim)+xlab("Support of the Distributions")+ylab("Density")

  p=p+theme(
    panel.background = element_rect(fill = "gray",
                                    colour = "gray",
                                    size = 0.05, linetype = "solid"),
    panel.grid.major = element_line(size = 0.02, linetype = 'solid',
                                    colour = "#f2f2f2"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                    colour = "#f2f2f2")
  )
  p
  #Fileloc="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/VonFCodes1.5/densityComparison plots/"
  #FileName=paste0(Fileloc,"plot_GIG_MPG_densities_nu_", nu_name, "alpha", alpha , ".pdf")
  #ggsave(file=FileName, height = 3,width=4)
  #dev.off()
  return(p)
}




#########################



#' Plot_MPG_GIG_CDF
#'
#' @examples
#' alpha=1
#' pp=Plot_MPG_GIG_CDF(alpha= alpha, nu=1, gigCompare=TRUE )
#' @export
Plot_MPG_GIG_CDF<-function(alpha,nu,gigCompare=TRUE, l_lim=NULL, r_lim=NULL ){

  mean_sd<-MPG_mean_sd(alpha = alpha , nu = nu)
  mean_pmg<-mean_sd$mean; sd_mpg<-mean_sd$sd

  if(is.null(l_lim)){
  l_lim<- max(mean_pmg-2.5*sd_mpg, .0001)
  }
  if(is.null(r_lim)){
  r_lim<- min(mean_pmg+6*sd_mpg, 1.5)
  }
  print(c(l_lim, r_lim))


  #library(ggplot2)

  plt<-ggplot()+xlim(l_lim, r_lim)
  plt<- plt +  stat_function(fun=pMPG, args = list(alpha=alpha,nu=nu,iflog=FALSE, Num_of_Terms=10000), colour = "white", size=.3, n=1000 )


  if(gigCompare){
    lambda_gig<- -nu-1; chi_gig<- .5; psi_gig<- 2*alpha^2
    #plt1<-ggplot(size=.25)+xlim(0.0001, 3)
    plt<- plt +  stat_function(fun=ghyp::pgig, args = list(lambda=lambda_gig,chi=chi_gig,psi=psi_gig ), colour = "black", size=.2, n=1000 )


  }




  plt=plt+theme(
    panel.background = element_rect(fill = "gray",
                                    colour = "gray",
                                    size = 0.03, linetype = "solid"),
    panel.grid.major = element_line(size = 0.01, linetype = 'solid',
                                    colour = "#d2d2d2"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                    colour = "#d2d2d2")
  )

 plt<-plt+  xlab("Support of the Distributions")+ ylab("Cumulative Distribution Function ")


  return(plt)

}


#detach_package("frmqa", TRUE)




pgig_new<-function(x, lambda,chi,psi ){
  val<-apply(X = as.matrix(x),MARGIN = 1, FUN = function(x) {ghyp::pgig(x, lambda=lambda,chi=chi,psi=psi)}  )
  return(val)
}


#
# detach_package <- function(pkg, character.only = FALSE)
# {
#   if(!character.only)
#   {
#     pkg <- deparse(substitute(pkg))
#   }
#   search_item <- paste("package", pkg, sep = ":")
#   while(search_item %in% search())
#   {
#     detach(search_item, unload = TRUE, character.only = TRUE)
#   }
# }






#' Plot_MPG_GIG_CDF
#' @param alpha: concentration parameter of the MPG distribution. A positive real number
#' @param nu: order of the MPG distribution. A positive real number
#' @examples
#' alpha=1
#' pp=Plot_MPG_GIG_CDF_alt(alpha= alpha, nu=1, gigCompare=TRUE )
#' @export
Plot_MPG_GIG_CDF_alt<-function(alpha,nu,gigCompare=TRUE, l_lim=NULL, r_lim=NULL ){

  mean_sd<-MPG_mean_sd(alpha = alpha , nu = nu)
  mean_pmg<-mean_sd$mean; sd_mpg<-mean_sd$sd

  if(is.null(l_lim)){
    l_lim<- max(mean_pmg-2.5*sd_mpg, .0001)
  }
  if(is.null(r_lim)){
    r_lim<- min(mean_pmg+6*sd_mpg, 1.5)
  }
  print(c(l_lim, r_lim))


  #library(ggplot2)

  plt<-ggplot()+xlim(l_lim, r_lim)
  plt<- plt +  stat_function(fun=pMPG, args = list(alpha=alpha,nu=nu,iflog=FALSE, Num_of_Terms=10000), colour = "white", size=.5, n=1000 )


  if(gigCompare){
    lambda_gig<- -nu-1; chi_gig<- .5; psi_gig<- 2*alpha^2
   # plt1<-ggplot(size=.25)
    plt<- plt +  stat_function(fun=ghyp::pgig, args = list(lambda=lambda_gig,chi=chi_gig,psi=psi_gig ), colour = "gray15", size=.4,linetype=  'dotdash', n=1000 )




      #lambda_gig<- -nu-1.5; chi_gig<- .5; psi_gig<- 2*alpha^2
      ##plt1<-ggplot(size=.25)+xlim(0.0001, 3)
      #plt<- plt +  stat_function(fun=ghyp::pgig, args = list(lambda=lambda_gig,chi=chi_gig,psi=psi_gig ), colour = "red", size=.2, n=1000 )


      #lambda_gig<- -nu-1-.5*(nu)/(nu+1); chi_gig<- .5; psi_gig<- 2*alpha^2
      #plt1<-ggplot(size=.25)+xlim(0.0001, 3)
      #plt<- plt +  stat_function(fun=ghyp::pgig, args = list(lambda=lambda_gig,chi=chi_gig,psi=psi_gig ), colour = "blue", size=.2, n=1000 )


      nu_calc<-matching_mean_GIG_nu_calculator(alpha = alpha, nu = nu)$root
      lambda_gig<- nu_calc; chi_gig<- .5; psi_gig<- 2*alpha^2
      #plt1<-ggplot(size=.25)+xlim(0.0001, 3)
      plt<- plt +  stat_function(fun=ghyp::pgig, args = list(lambda=lambda_gig,chi=chi_gig,psi=psi_gig ), colour = "black", size=.3, n=1000 )



  }




  plt=plt+theme(
    panel.background = element_rect(fill = "gray",
                                    colour = "gray",
                                    size = 0.03, linetype = "solid"),
    panel.grid.major = element_line(size = 0.01, linetype = 'solid',
                                    colour = "#d2d2d2"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                    colour = "#d2d2d2")
  )

  plt<-plt+  xlab("Support of the Distributions")+ ylab("Cumulative Distribution Function ")


  return(plt)

}


#detach_package("frmqa", TRUE)






#' matching_mean_GIG_nu_calculator
#' @param alpha: a positive real number
#' @param nu: A positive real number of the Ratio bessel \deqn{\frac{I_{\nu+1}(\alpha)}{I_{\nu}(\alpha)}}
#' @return opt_nu_for_GIG: An uniroot object for the optimum. opt_nu_for_GIG$root is the optimum nu_1 that satisfies the equation $\frac{K_{\nu_1+1}(\alpha)}{K_{\nu_1}(\alpha)}= \frac{I_{\nu+1}(\alpha)}{I_{\nu}(\alpha)}$
#' @examples
#' alpha=1
#' matching_nu_calculator(alpha=10, nu=1 )
#' @export
matching_mean_GIG_nu_calculator<-function(alpha, nu){
  x=alpha
  val= BesselI_Ratio(x, nu1 = nu)
  f_to_solve_K_ratio<-function( nu){
    #nu1=-nu-1
    val1= BesselK_Ratio(x = x, nu1 = nu)- val
    return(val1)
  }

  init_val_left=-nu-.8; init_val_right=-nu-3
 sol_root<- uniroot(f_to_solve_K_ratio, interval = c(init_val_left, init_val_right ), tol = 0.0000001)
opt_nu_for_GIG<-sol_root

  return(opt_nu_for_GIG)
}










#######



#' @export
plot_GIG_MPG_densities_alt<-function(alpha=5, nu=0 , UpperLim=NULL, LowerLim=NULL){
  f_mpg<-function(x){
    if(length(x)>1){
      return(apply(matrix(x, ncol=1), MARGIN = 1, function(y){log_f_nu(t=y, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = FALSE)}))
    }
    return(log_f_nu(t=x, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = FALSE))

  }


  f_gig_alt<-function(x){
    nu_calc<-matching_mean_GIG_nu_calculator(alpha = alpha, nu = nu)$root
    lambda_gig<- nu_calc; chi_gig<- .5; psi_gig<- 2*alpha^2
    val=dgig(x =x, lambda =lambda_gig, chi = chi_gig, psi = psi_gig )
    return(val)
  }


  f_gig<-function(x){
    val=dgig(x =x, lambda =-(nu+1), chi = .5, psi = 2*alpha^2 )
    return(val)
  }
  #nu=.5
  #alpha=5;
  j_0=bessel_zero_Jnu(nu = nu, s = 1:1000)
  J_1=besselJ(j_0,nu=1+nu);
  if(is.null(UpperLim)){
    UpperLim=.05+.75/alpha
  }
  if(is.null(LowerLim)){LowerLim=.006}
  #UpperLim=.25*(nu==0)+ .15*(nu==.5)+ .1*(nu>.5)*(nu<4)+ .06*(nu>4)
  library(ggplot2)
  nu_name=nu
  if(nu==.5){nu_name="half"}
  #browser()
  print(alpha)


  p=ggplot(data.frame(x = rnorm(100)), aes(x)) +
    stat_function(fun = function(x){f_mpg(x = x )}, colour = "white", size=.3)+
    stat_function(fun = function(x){f_gig(x = x)}, colour = "black", size=.2)+
    stat_function(fun = function(x){f_gig_alt(x = x)}, colour = "yellow", size=.2)+
    xlim(LowerLim, UpperLim)+xlab("Support of the Distributions")+ylab("Density")

  p=p+theme(
    panel.background = element_rect(fill = "gray",
                                    colour = "gray",
                                    size = 0.05, linetype = "solid"),
    panel.grid.major = element_line(size = 0.02, linetype = 'solid',
                                    colour = "#f2f2f2"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                    colour = "#f2f2f2")
  )
  p
  #Fileloc="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/VonFCodes1.5/densityComparison plots/"
  #FileName=paste0(Fileloc,"plot_GIG_MPG_densities_nu_", nu_name, "alpha", alpha , ".pdf")
  #ggsave(file=FileName, height = 3,width=4)
  #dev.off()
  return(p)
}




############################




#' CalculateTV_between_MPG_GIG_CDF
#' @param alpha: a positive real number
#' @param nu: A positive real number of the Ratio bessel \deqn{\frac{I_{\nu+1}(\alpha)}{I_{\nu}(\alpha)}}
#' @return Total Variation Distribution between MPG(alpha, nu) and the fowwowing distributions:
#' GIG(lambda=-(nu+1), .5, 2\alpha^2), GIG(lambda=\nu_{opt}, .5, 2\alpha^2)
#' @examples
#' alpha=1
#' CalculateTV_between_MPG_GIG_CDF(alpha=10, nu=1 )
#' @export
CalculateTV_between_MPG_GIG_CDF<-function(alpha, nu, search_n=100,l_lim=NULL, r_lim=NULL  ){
#browser()
  mean_sd<-MPG_mean_sd(alpha = alpha , nu = nu)
  mean_pmg<-mean_sd$mean; sd_mpg<-mean_sd$sd
  if(is.null(l_lim)){
    l_lim<- max(mean_pmg-1.8*sd_mpg, .0001)
  }
  if(is.null(r_lim)){
    r_lim<- min(mean_pmg+3*sd_mpg, 1.5)
  }
  #print(c(l_lim, r_lim))
  x_val=seq(from = l_lim, to = r_lim, length.out = search_n )



    MPG_CDF_val<- pMPG(x_val, alpha = alpha, nu = nu)


    lambda_gig<- -nu-1; chi_gig<- .5; psi_gig<- 2*alpha^2
    GIG_nu_plus_1_CDF_val<-ghyp::pgig(q = x_val,lambda =lambda_gig, chi =chi_gig, psi = psi_gig   )


    nu_calc<-matching_mean_GIG_nu_calculator(alpha = alpha, nu = nu)$root
    lambda_gig<- nu_calc; chi_gig<- .5; psi_gig<- 2*alpha^2
    GIG_nu_opt_CDF_val<-ghyp::pgig(q = x_val,lambda =lambda_gig, chi =chi_gig, psi = psi_gig   )

  TV_MPG_GIG_nu_plus_1 <-max(MPG_CDF_val-GIG_nu_plus_1_CDF_val)- min(0, min((MPG_CDF_val-GIG_nu_plus_1_CDF_val)))
  TV_MPG_GIG_nu_opt    <- max(MPG_CDF_val-GIG_nu_opt_CDF_val)- min(0, min((MPG_CDF_val-GIG_nu_opt_CDF_val)))

  #TV_MPG_GIG_nu_plus_1 <-max(MPG_CDF_val-GIG_nu_plus_1_CDF_val)
  #TV_MPG_GIG_nu_opt    <- max(MPG_CDF_val-GIG_nu_opt_CDF_val)


  Total_Var_Dist<-list(nu_opt= nu_calc, TV_MPG_GIG_nu_opt=TV_MPG_GIG_nu_opt, TV_MPG_GIG_nu_plus_1=TV_MPG_GIG_nu_plus_1)
  return(Total_Var_Dist)
}




#' CalculateTV_between_MPG_GIG_based_on_pdf
#' @param alpha: a positive real number
#' @param nu: A positive real number of the Ratio bessel \deqn{\frac{I_{\nu+1}(\alpha)}{I_{\nu}(\alpha)}}
#' @return Total Variation Distribution between MPG(alpha, nu) and the fowwowing distributions:
#' GIG(lambda=-(nu+1), .5, 2\alpha^2), GIG(lambda=\nu_{opt}, .5, 2\alpha^2)
#' @examples
#' alpha=1
#' CalculateTV_between_MPG_GIG_based_on_pdf(alpha=10, nu=1 )
#' @export
CalculateTV_between_MPG_GIG_based_on_pdf<-function(alpha, nu, search_n=100,l_lim=NULL, r_lim=NULL  ){

  mean_sd<-MPG_mean_sd(alpha = alpha , nu = nu)
  mean_pmg<-mean_sd$mean; sd_mpg<-mean_sd$sd
  if(is.null(l_lim)){
    l_lim<- max(mean_pmg-2.2*sd_mpg, .0001)
  }
  if(is.null(r_lim)){
    r_lim<- min(mean_pmg+5*sd_mpg, 1.5)
  }
  #print(c(l_lim, r_lim))
  #x_val=seq(from = l_lim, to = r_lim, length.out = search_n )
  lambda_gig<- -nu-1; chi_gig<- .5; psi_gig<- 2*alpha^2
  f_gig_mpg<-function(x){
    abs(dMPG(x, alpha = alpha, nu = nu,iflog = FALSE)- ghyp::dgig(x =  x,lambda =lambda_gig, chi =chi_gig, psi = psi_gig   ))
  }

  TV_mpg_gig<-.5*(integrate(f = f_gig_mpg, lower =l_lim, upper = r_lim, subdivisions = search_n )$value)


  nu_calc<-matching_mean_GIG_nu_calculator(alpha = alpha, nu = nu)$root
  lambda_gig<- nu_calc; chi_gig<- .5; psi_gig<- 2*alpha^2
  f_gig_opt_mpg<-function(x){
    abs(dMPG(x = x, alpha = alpha, nu = nu, iflog = FALSE)- ghyp::dgig(x =  x,lambda =lambda_gig, chi =chi_gig, psi = psi_gig   ))
  }
  TV_mpg_gig_opt<-.5*(integrate(f = f_gig_opt_mpg, lower =l_lim, upper = r_lim, subdivisions = search_n )$value)

  #MPG_CDF_val<- pMPG(x_val, alpha = alpha, nu = nu)



  #GIG_nu_plus_1_CDF_val<-ghyp::pgig(q = x_val,lambda =lambda_gig, chi =chi_gig, psi = psi_gig   )



  #GIG_nu_opt_CDF_val<-ghyp::pgig(q = x_val,lambda =lambda_gig, chi =chi_gig, psi = psi_gig   )







  TV_MPG_GIG_nu_plus_1 <-TV_mpg_gig
  TV_MPG_GIG_nu_opt    <- TV_mpg_gig_opt
  Total_Var_Dist<-list(nu_opt= nu_calc, TV_MPG_GIG_nu_opt=TV_MPG_GIG_nu_opt, TV_MPG_GIG_nu_plus_1=TV_MPG_GIG_nu_plus_1)
  return(Total_Var_Dist)
}





###################################
 # The table of TV generator
# alpha_seq= c(10, 15, 20)
# nu_seq= c(0, 5, 10)
#cat(generate_Table_TV_dist(alpha_seq = alpha_seq, nu_seq = nu_seq)$lat_string)


# alpha_seq= c(5, 10, 15, 20, 25, 30 )
# nu_seq= c(0, .5, 10, 15, 20, 25, 30  )

generate_Table_TV_dist<-function(alpha_seq, nu_seq, search_n=5000, ifLatex=TRUE){
  row_num= length(alpha_seq)* length( nu_seq)
  TV_table=matrix(NA, nrow=row_num, ncol=5)
  row_ind=1
  for(i in 1:length(nu_seq)){
    for(j in 1:length( alpha_seq)){
          list_val<-CalculateTV_between_MPG_GIG_CDF(alpha=alpha_seq[j], nu=nu_seq[i], search_n = search_n)
          TV_table[row_ind, 1]= alpha_seq[j]
          TV_table[row_ind, 2]= nu_seq[i]
          TV_table[row_ind, 3]= round( list_val$TV_MPG_GIG_nu_plus_1, 6)
          TV_table[row_ind, 4]= round(list_val$nu_opt, 3)
          TV_table[row_ind, 5]= round(list_val$TV_MPG_GIG_nu_opt, 6)
          row_ind=row_ind+1
    }
  }
  TV_table_lst=TV_table
  if(ifLatex){
    lat_string<-apply(X = TV_table,MARGIN = 1, FUN = function(x){paste0(paste(x, collapse = "& "), "\\\\\\hline\n")} )
    TV_table_lst<- list(Tab=TV_table, lat_string=lat_string )
  }

  return(TV_table_lst)


}




###################################
# The table of TV generator
# alpha_seq= c(10, 15, 20)
# nu_seq= c(0, 5, 10)
#cat(generate_Table_TV_dist(alpha_seq = alpha_seq, nu_seq = nu_seq)$lat_string)


# alpha_seq= c(5, 10, 15, 20, 25, 30 )
# nu_seq= c(0, .5, 10, 15, 20, 25, 30  )
generate_Table_TV_dist_pdf_based<-function(alpha_seq, nu_seq, search_n=1000, ifLatex=TRUE){
  row_num= length(alpha_seq)* length( nu_seq)
  TV_table=matrix(NA, nrow=row_num, ncol=5)
  row_ind=1
  for(i in 1:length(nu_seq)){
    for(j in 1:length( alpha_seq)){
      list_val<-CalculateTV_between_MPG_GIG_based_on_pdf(alpha=alpha_seq[j], nu=nu_seq[i], search_n = search_n)
      TV_table[row_ind, 1]= alpha_seq[j]
      TV_table[row_ind, 2]= nu_seq[i]
      TV_table[row_ind, 3]= round( list_val$TV_MPG_GIG_nu_plus_1, 6)
      TV_table[row_ind, 4]= round(list_val$nu_opt, 3)
      TV_table[row_ind, 5]= round(list_val$TV_MPG_GIG_nu_opt, 6)
      row_ind=row_ind+1
    }
  }
  TV_table_lst=TV_table
  if(ifLatex){
    lat_string<-apply(X = TV_table,MARGIN = 1, FUN = function(x){paste0(paste(x, collapse = "& "), "\\\\\\hline\n")} )
    TV_table_lst<- list(Tab=TV_table, lat_string=lat_string )
  }

  return(TV_table_lst)


}




######################################################################################################

#' Plot_MPG_GIG_CDF
#' @param alpha: concentration parameter of the MPG distribution. A positive real number
#' @param nu: order of the MPG distribution. A positive real number
#' @examples
#' alpha=1
#' pp=Plot_MPG_GIG_CDF_alt(alpha= alpha, nu=1, gigCompare=TRUE )
#' @export
Plot_MPG_CDF_discrete_approx<-function(alpha,nu,gigCompare=TRUE, l_lim=NULL, r_lim=NULL ){

  mean_sd<-MPG_mean_sd(alpha = alpha , nu = nu)
  mean_pmg<-mean_sd$mean; sd_mpg<-mean_sd$sd

  if(is.null(l_lim)){
    l_lim<- max(mean_pmg-2.5*sd_mpg, .0001)
  }
  if(is.null(r_lim)){
    r_lim<- min(mean_pmg+6*sd_mpg, 1.5)
  }
  print(c(l_lim, r_lim))


  #library(ggplot2)

  plt<-ggplot()+xlim(l_lim, r_lim)
  plt<- plt +  stat_function(fun=pMPG, args = list(alpha=alpha,nu=nu,iflog=FALSE, Num_of_Terms=10000), colour = "white", size=.5, n=1000 )


  plt<- plt +  stat_function(fun=pMPG, args = list(alpha=alpha,nu=nu,iflog=FALSE, Num_of_Terms=30), colour = "black", size=.5, n=1000 )

  plt<- plt +  stat_function(fun=pMPG, args = list(alpha=alpha,nu=nu,iflog=FALSE, Num_of_Terms=40), colour = "gray50", size=.5, n=1000 )


  plt=plt+theme(
    panel.background = element_rect(fill = "gray",
                                    colour = "gray",
                                    size = 0.03, linetype = "solid"),
    panel.grid.major = element_line(size = 0.01, linetype = 'solid',
                                    colour = "#d2d2d2"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                    colour = "#d2d2d2")
  )

  plt<-plt+  xlab("Support of the Distributions")+ ylab("Cumulative Distribution Function ")


  return(plt)

}






######################################################################################
######################################################################################
######################################################################################

Bin_Avg_CDF<-function(stat_x_bar,n, p){
  x=n*(stat_x_bar/sqrt(n)+p)
  prob<-pbinom(q = x, size = n, prob = p)
  return(prob)
}


#p=.5
#library(ggplot2)
#plt<-ggplot()+stat_function(fun =Bin_Avg_CDF, args = list(n=30, p=.5) )+xlim(-2,2)
#plt<-plt+stat_function(fun =pnorm, args = list(mean=0, sd= sqrt(p*(1-p))), col='red' )

TV_dist<-function(size, p, n_grid=1000 ){
  x_seq=seq(-2, 2, length.out =n_grid )
  bin_cdf     <-  Bin_Avg_CDF(stat_x_bar = x_seq,n = size, p = p )
  norm_cdf    <-  pnorm(q = x_seq, mean = 0, sd =sqrt(p*(1-p))  )
  return(max(-bin_cdf+norm_cdf))
}



