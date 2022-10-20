# Ants Data Analysis


calculate_radian<-function(dd){
  val<-atan2( dd[2], dd[1])
  if(val<0){ val=val+2*pi   }
  return(val)
}



library(CircNNTSR)
data(Ants)
data=cbind(x=cos(Ants*pi/180), y=sin(Ants*pi/180))
MC_OBject=VONF_pD_MPG_INVCDF_DA_POSTERIOR(Y=data, MCSamplerSize=5000)
plt<-Plot_MCMC_Diag_Triplet(MC_OBject$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plot = plt, filename ="/Users/subhadippal/Box/DA for Directional data/WriteUp/6.6/figures/Ants_kappa_MC.pdf", height = 6 , width = 6 )

post_rad<-apply(MC_OBject$McSample_mu, MARGIN = 2, calculate_radian)
plt1<-Plot_MCMC_Diag_Triplet(post_rad, y_lab_text = 'modal direction (in radians)' )
ggsave(plot = plt1, filename ="/Users/subhadippal/Box/DA for Directional data/WriteUp/6.6/figures/Ants_mu_MC.pdf", height = 6 , width = 6 )


