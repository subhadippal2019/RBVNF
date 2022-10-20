

file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/ACF_Sim1_3D_MPG_DA_vioplot_cr.pdf"
pdf(file=file, width = 6.5, height =3.5 )

data("ACF_Sim1_3D_MPG_DA_vioplot_data")
plot_Acf_violiin(Vio_ACF_MPG, size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()




file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/ACF_Sim1_3D_GEO_DA_vioplot_cr.pdf"
pdf(file=file, width = 6.5, height =3.5 )

data(ACF_Sim1_3D_GEO_DA_vioplot_data)
plot_Acf_violiin(Vio_ACF_GEO, size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()









file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/RealDataCode/MCMC_sample_NC_12d_DATAID_26022021.RData"

load(file)


file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/RealDataCode/MCMC_sample_NC_12d_DATAID_26022021.RData"

load(file)

length(MCMC_sample)



# A moc run of DIC
calculate_DIC(MCMC_sample, data =t(sel_DTI_real_DATA) , burn_in = 1000)






file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/ACF_Sim2_pD_MPG_DA_vioplot.pdf"
pdf(file=file, width = 6.5, height =3.5 )
data(ACF_sim2_pD_MPG_DA_vioplot_data)
plot_Acf_violiin(Vio_ACF_MPG_pd,size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()




file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/SLICE_ACF_vioplot.pdf"
#pdf(file=file, width = 6.5, height =3.5 )
data("ACF_Slice")
plot_Acf_violiin(ACF_Slice)
plt=plot_Acf_violiin(ACF_Slice,size=.2, ifGray=TRUE, GridLines="Vertical", Ylimit = c(.9, 1))
Sys.sleep(5)
ggsave(filename = file, plot=plt,width = 6.5, height =3.5)
#Sys.sleep(5)
#dev.off()



file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/Circular_MPG_ACF_vioplot1.pdf"

dd=circular_data_MPG_acf(25)
plt=plot_Acf_violiin(dd,size=.2, ifGray=TRUE, GridLines="Vertical", Ylimit = c(0, 1))
Sys.sleep(5)
ggsave(filename = file, plot=plt,width = 6.5, height =3.5)

tab=CalculateTable_MPG_circular(dd)
kappa_str=paste0("$\\kappa= ", c(1,5,10,15), "$")
tab1=cbind(kappa_str, tab$TableFormat1)
aa=apply(tab1,MARGIN = 1 ,FUN = function(x){ paste0( paste0(x, collapse = "&"), "\\\\ \n \\hline")})
cat(aa)




#####################  Simulation Summary Codes #####################################################

load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_5.RData")


plt=Plot_MCMC_Diag_Triplet(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/GEO_DA_Kappa_5_3d.pdf", width=6, height=6)


load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_10.RData")

plt=Plot_MCMC_Diag_Triplet(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/GEO_DA_Kappa_10_3d.pdf", width=6, height=6)


plt=Plot_MCMC_Diag_Triplet(AllList_3D[[1]]$lst_MPG_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MPG_DA_Kappa_10_3d.pdf", width=6, height=6)




#########################################
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_1.RData")

plt=Plot_MCMC_Diag_Triplet(AllList_3D[[2]]$lst_MPG_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MPG_DA_Kappa_1_3d.pdf", width=6, height=6)


plt=Plot_MCMC_Diag_Triplet(AllList_3D[[2]]$lst_GEO_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/GEO_DA_Kappa_1_3d.pdf", width=6, height=6)

#####
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_1.RData")

plt=Plot_MCMC_Diag_Triplet(AllList_3D[[2]]$lst_MPG_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MPG_DA_Kappa_1_3d.pdf", width=6, height=6)


plt=Plot_MCMC_Diag_Triplet(AllList_3D[[2]]$lst_GEO_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/GEO_DA_Kappa_1_3d.pdf", width=6, height=6)

#####
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_0.01.RData")

plt=Plot_MCMC_Diag_Triplet(AllList_3D[[2]]$lst_MPG_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MPG_DA_Kappa_01_3d.pdf", width=6, height=6)


plt=Plot_MCMC_Diag_Triplet(AllList_3D[[2]]$lst_GEO_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/GEO_DA_Kappa_01_3d.pdf", width=6, height=6)

#####
load("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/AllList_3D_Kappa_eq_15.RData")

plt=Plot_MCMC_Diag_Triplet(AllList_3D[[4]]$lst_MPG_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/MPG_DA_Kappa_15_3d.pdf", width=6, height=6)


plt=Plot_MCMC_Diag_Triplet(AllList_3D[[5]]$lst_GEO_3D$McSample_kappa, y_lab_text = expression(kappa))

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/GEO_DA_Kappa_15_3d.pdf", width=6, height=6)







dic<-read.csv("/Volumes/NO NAME/DTIRealData/DIC_new_2021_5000.csv")
library(ggplot2)
# Basic line plot with points
plt<-ggplot(data=data.frame(dic[c(1:7, 9:10), ]), aes(x=lag_index, y=Store_DIC, group=1)) +
  geom_line()+
  geom_point(size=6, alpha=.3)+
  geom_point(size=2, alpha=1, col='white')+
  geom_point(size=1, alpha=.8, col='black')+
  xlab("Number of Clusters")+
  ylab("DIC")

ggsave(plt, file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/DIC.pdf", width=6, height=3.5)



load(file = "/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/RealDataCode/MCMC_sample_NC_9_DATAID_25022021.RData")
plt1=Plot_MCMC_Diag_Triplet_multiple(MCMC = MCMC_sample, NUM = 4)









