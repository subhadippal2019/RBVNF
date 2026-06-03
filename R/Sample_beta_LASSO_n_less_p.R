


r_sample_beta_n_less_p_LASSO<-function(X, Y, T_aug,Tau_ij_sq_all, d, n ){
#T_aug
#Diag_T

#d;n

iMat_d<-diag(rep(1,d))
#D_tau_inv=diag(c(1/Tau_ij_sq_all))
y_dagger=c(Y)
#Diag_T=diag(T_aug)

#D_tau= diag(c(Tau_ij_sq_all))


U1=rnorm(n = p*d, mean = 0, sd = sqrt(Tau_ij_sq_all))
U2=rnorm(n=n*d,mean = 0, sd = 1 )
#X_star= sqrt(2)*sqrt(Diag_T)%*%X # Next line does the same thing more efficeintly
X_star= sqrt(2)*sqrt(T_aug)*X ####

#M=kronecker(iMat_d,X_star ) %*% D_tau%*% kronecker(iMat_d,t(X_star) )+diag(d*n) ## Next line does the same thing more efficeintly
M=kronecker(iMat_d,X_star ) %*% (c(Tau_ij_sq_all)* kronecker(iMat_d,t(X_star) ))+diag(d*n) #####
#browser()
V1= kronecker(iMat_d,X_star )%*%U1+U2
V2= solve(M)%*%( kronecker(iMat_d, diag(1/sqrt(2*c(T_aug)))   ) %*%y_dagger  -V1)

#beta_dagger=U1+D_tau%*%kronecker(iMat_d,t(X_star) ) %*% V2## Next line does the same thing more efficeintly
beta_dagger=U1+c(Tau_ij_sq_all)*kronecker(iMat_d,t(X_star) ) %*% V2
#cbind(beta_dagger, beta_vec)
#browser()
return(beta_dagger)
}
