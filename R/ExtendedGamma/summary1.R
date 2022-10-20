load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_5.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

kappa_est_GEO=c()
for( i in 1:20){
  kappa_est_GEO[i]=mean(AllList_3D[[i]]$lst_GEO_3D$McSample_kappa[3000:5000])
}

MSE_GEO_5=sum((kappa_est_GEO-5)^2/20)


kappa_est_MPG=c()
for( i in 1:20){
  kappa_est_MPG[i]=mean(AllList_3D[[i]]$lst_MPG_3D$McSample_kappa[3000:5000])
}

MSE_MPG_5=sum((kappa_est_MPG-5)^2/20)
##################################################################################
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_10.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

kappa_est_GEO=c()
for( i in 1:20){
  kappa_est_GEO[i]=mean(AllList_3D[[i]]$lst_GEO_3D$McSample_kappa[3000:5000])
}

MSE_GEO_10=sum((kappa_est_GEO-10)^2/20)


kappa_est_MPG=c()
for( i in 1:20){
  kappa_est_MPG[i]=mean(AllList_3D[[i]]$lst_MPG_3D$McSample_kappa[3000:5000])
}

MSE_MPG_10=sum((kappa_est_MPG-10)^2/20)

#####################################################################################

##################################################################################
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_15.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

kappa_est_GEO=c()
for( i in 1:20){
  kappa_est_GEO[i]=mean(AllList_3D[[i]]$lst_GEO_3D$McSample_kappa[3000:5000])
}

MSE_GEO_15=sum((kappa_est_GEO-15)^2/20)


kappa_est_MPG=c()
for( i in 1:20){
  kappa_est_MPG[i]=mean(AllList_3D[[i]]$lst_MPG_3D$McSample_kappa[3000:5000])
}

MSE_MPG_15=sum((kappa_est_MPG-15)^2/20)

############################
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_1.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

kappa_est_GEO=c()
for( i in 1:20){
  kappa_est_GEO[i]=mean(AllList_3D[[i]]$lst_GEO_3D$McSample_kappa[3000:5000])
}

MSE_GEO_1=sum((kappa_est_GEO-1)^2/20)
############################
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_0.01.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

kappa_est_GEO=c()
for( i in 1:20){
  kappa_est_GEO[i]=mean(AllList_3D[[i]]$lst_GEO_3D$McSample_kappa[3000:5000])
}

MSE_GEO_01=sum((kappa_est_GEO-.1)^2/20)


kappa_est_MPG=c()
for( i in 1:20){
  kappa_est_MPG[i]=mean(AllList_3D[[i]]$lst_MPG_3D$McSample_kappa[3000:5000])
}

MSE_MPG_01=sum((kappa_est_MPG-.1)^2/20)


MSE_mat<-cbind(kappa_tru=c(0.1,1,5,10,15),GEO_MSE=c(MSE_GEO_01, MSE_GEO_1, MSE_GEO_5, MSE_GEO_10, MSE_GEO_15), MPG_MSE=c(MSE_GEO_01, MSE_MPG_1, MSE_MPG_5, MSE_MPG_10, MSE_MPG_15))



apply(round(t(MSE_mat), 3), 1, function(x){paste(x[1], " & ", x[2], "& ",x[3] ,"& ",x[4], "& ",x[5],  "\\\\hline")})




##################### mu estimates

load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_5.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

MSE_term_Mu_est_GEO=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_GEO_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  MSE_term_Mu_est_GEO[i ]=(xx[1])
}

MSE_GEO_5=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((MSE_term_Mu_est_GEO-1)^2/20)


Mu_est_MPG=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_MPG_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
Mu_est_MPG[i]=xx[1]
}

MSE_MPG_5=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((Mu_est_MPG-1)^2/20)
##################################################################################
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_10.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

MSE_term_Mu_est_GEO=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_GEO_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  MSE_term_Mu_est_GEO[i ]=(xx[1])
}

MSE_GEO_10=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((MSE_term_Mu_est_GEO-1)^2/20)


Mu_est_MPG=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_MPG_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  Mu_est_MPG[i]=xx[1]
}

MSE_MPG_10=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((Mu_est_MPG-1)^2/20)


#####################################################################################

##################################################################################
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_15.RData")

MSE_term_Mu_est_GEO=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_GEO_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  MSE_term_Mu_est_GEO[i ]=(xx[1])
}

MSE_GEO_15=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((MSE_term_Mu_est_GEO-1)^2/20)


Mu_est_MPG=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_MPG_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  Mu_est_MPG[i]=xx[1]
}

MSE_MPG_15=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((Mu_est_MPG-1)^2/20)



############################
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_1.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

MSE_term_Mu_est_GEO=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_GEO_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  MSE_term_Mu_est_GEO[i ]=(xx[1])
}

MSE_GEO_1=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((MSE_term_Mu_est_GEO-1)^2/20)


Mu_est_MPG=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_MPG_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  Mu_est_MPG[i]=xx[1]
}

MSE_MPG_1=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((Mu_est_MPG-1)^2/20)




load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_0.01.RData")


#plot(density(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa[1:100]),main='MPG Augmentation')

MSE_term_Mu_est_GEO=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_GEO_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  MSE_term_Mu_est_GEO[i ]=(xx[1])
}

MSE_GEO_01=sum((1-abs(MSE_term_Mu_est_GEO))/20)


Mu_est_MPG=c()
for( i in 1:20){
  xx=apply(AllList_3D[[i]]$lst_MPG_3D$McSample_mu[,3000:5000], 1, 'mean')
  xx=xx/norm_vec(xx)
  Mu_est_MPG[i]=xx[1]
}

MSE_MPG_01=sum((1-abs(MSE_term_Mu_est_GEO))/20)#sum((Mu_est_MPG-1)^2/20)






MSE_mat<-cbind(kappa_tru=c(0.01,1,5,10,15),GEO_MSE=c(MSE_GEO_01, MSE_GEO_1, MSE_GEO_5, MSE_GEO_10, MSE_GEO_15), MPG_MSE=c(MSE_GEO_01, MSE_MPG_1, MSE_MPG_5, MSE_MPG_10, MSE_MPG_15))



aa=apply(round(t(MSE_mat), 5), 1, function(x){paste(x[1], " & ", x[2], "& ",x[3] ,"& ",x[4], "& ",x[5],  "\\\\hline")})


cat(aa)




###################################################################################



acf_value=cbind( acf= acf(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa,lag=10)$acf,lag=as.factor(1:11))
for(i in 2:20){
acf_value= rbind(acf_value,cbind( acf= acf(AllList_3D[[i]]$lst_GEO_3D$McSample_kappa,lag=10)$acf,lag=as.factor(1:11)))

}



plot(density(acf_value[3,]))

acf_value=data.frame(acf_value)
library(vioplot)
x1 <- acf_value$acf[acf_value$lag==2]
x2 <- acf_value$acf[acf_value$lag==3]
x3 <- acf_value$acf[acf_value$lag==4]
x4 <- acf_value$acf[acf_value$lag==5]
x5 <- acf_value$acf[acf_value$lag==6]
x6 <- acf_value$acf[acf_value$lag==7]
x7 <- acf_value$acf[acf_value$lag==8]
x8 <- acf_value$acf[acf_value$lag==9]
x9 <- acf_value$acf[acf_value$lag==10]
x10 <- acf_value$acf[acf_value$lag==10]

vioplot(x1, x2, x3,x4,x5,x6,x7,x8,x9,x10,  col="magenta")
title("Violin Plots of Miles Per Gallon")

library(ggplot2)
# Basic violin plot
p <- ggplot(acf_df, aes(factor(acf_df$lag),acf_df$acf)) + geom_violin()
p
# Rotate the violin plot
p + coord_flip()
# Set trim argument to FALSE
ggplot(ToothGrowth, aes(x=dose, y=len)) +    geom_violin(trim=FALSE)

ToothGrowth
