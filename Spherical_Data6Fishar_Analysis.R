
#Datab6fisher Analysis Spherical data

library(CircNNTSR)
data(Datab6fisher)

phi<-Datab6fisher$V2*pi/180
theta<-Datab6fisher$V1*pi/180



calculate_radian<-function(dd){
  val<-atan2( dd[2], dd[1])
  if(val<0){ val=val+2*pi   }
  return(val)
}

calculate_radians<-function(MC_Sample_mu){
  x=MC_Sample_mu[1, ]
  y=MC_Sample_mu[2, ]
  z=MC_Sample_mu[3, ]

  #tan(phi)=y/x
  dd=cbind(x,y)
  val1<-apply(dd, MARGIN = 1, FUN = calculate_radian)
  val1_1=atan(y/x)
  val2_2= atan(sqrt(x^2+y^2)/z)
  #val2_2= ifelse( sign(x)==1, atan(sqrt(x^2+y^2)/z),  ifelse( sign(y)==1,    pi+atan(sqrt(x^2+y^2)/z ),-pi+atan(sqrt(x^2+y^2)/z )    ))
  #val3<-acos(x/sin(val1))
return()
}



x=cos(phi)*sin(theta)
y =sin(phi)*sin(theta)
z= cos(theta)







library(rgl)
spheres3d(x,y,z,col="red",radius=0.02)
data=cbind(x=x, y=y, z=z)

MC_OBject_geo1<-VONF_3D_GEO_DA_POSTERIOR(Y = data, MCSamplerSize = 5000)
plt<-Plot_MCMC_Diag_Triplet(MC_OBject_geo$McSample_kappa, y_lab_text = expression(kappa))


ggsave(plot = plt, filename ="/Users/subhadippal/Box/DA for Directional data/WriteUp/6.6/figures/Datab6fisher_kappa_MC_1.pdf", height = 4.5 , width = 6 )




MC_OBject=VONF_pD_MPG_INVCDF_DA_POSTERIOR(Y=data, MCSamplerSize=5000)
Plot_MCMC_Diag_Triplet(MC_OBject$McSample_kappa)
#data(Datab3fisher)
#phi<-Datab3fisher$V2
#theta<-Datab3fisher$V1



#post_rad<-apply(MC_OBject$McSample_mu, MARGIN = 2, calculate_radian)
#Plot_MCMC_Diag_Triplet(post_rad)
