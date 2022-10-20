#library(fmri)

#file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Real Data/448347_txx.nii.gz"
#ddd<-read.NIFTI(filename = file)


library(oro.nifti)

#img_grid_sel1=img_grid_sel1[1:img_grid_sel, ]
Extract_FA_Dir<-function(xx){
  #browser()
  X_mat<-cbind( c(xx[1:3]) ,  c(xx[2], xx[4], xx[5] ) , c( xx[3], xx[5], xx[6]    )  )
  norm_X_sq=sum(xx^2)
  if(norm_X_sq ==0){
    FA=0; P_eigen_vec=c(0,0,0)
  }
  if(norm_X_sq >0){
    ev_decom=eigen(X_mat)
    ev=  ev_decom$values
    P_eigen_vec=ev_decom$vectors[,1]
    FA=sd(ev)/norm_vec(ev);
  }
  return( c(P_eigen_vec , FA )  )
  #xx=c(1, 2, 3, 4, 5, 6)
}




file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Real Data/448347_txx.nii.gz"
dxx<-readNIfTI(fname =  file)
file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Real Data/448347_txy.nii.gz"
dxy<-readNIfTI(fname =  file)
file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Real Data/448347_txz.nii.gz"
dxz<-readNIfTI(fname =  file)
file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Real Data/448347_tyy.nii.gz"
dyy<-readNIfTI(fname =  file)
file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Real Data/448347_tyz.nii.gz"
dyz<-readNIfTI(fname =  file)
file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Real Data/448347_tzz.nii.gz"
dzz<-readNIfTI(fname =  file)


#######################################################################
#######################################################################
img_size = dim(dzz)
img_grid = matrix( 0, nrow=img_size[1]*img_size[2]*img_size[3], ncol=6)

img_grid[,1]=c(dxx);
img_grid[,2]=c(dxy);
img_grid[,3]=c(dxz);
img_grid[,4]=c(dyy);
img_grid[,5]=c(dyz);
img_grid[,6]=c(dzz)


sel_id= which(img_grid[,1]!=0)
img_grid_sel=(img_grid[sel_id, ])


#######################################################################
#######################################################################


img_grid_sel1=0*(img_grid_sel[,1:4])

for(i in 1:nrow(img_grid_sel1) ){
  img_grid_sel1[i, ]=Extract_FA_Dir(img_grid_sel[i, ])
  if( !(i %% 10000)) print(i)
}


final_sel_id<-sel_id[(img_grid_sel1[,4]>.3)]

#########################################################################
#########################################################################

DTI_list<-list(final_sel_id=sel_id, Dti_data=img_grid_sel1, Original_data=img_grid, main_img_dim=dim(dzz)   )


save(DTI_list, file="DTI_processed_data.RData" )




#### cutoff point for FA is .3
