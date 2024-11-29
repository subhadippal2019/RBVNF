

soft_thresholding = function(x,a){
  sign(x) * pmax(abs(x)-a,0)
}




lasso_coord_desc = function(X,y,beta,lambda,tol=1e-6,maxiter=1000){
  #browser()
  beta = as.matrix(beta)
  X = as.matrix(X)
  omega = rep(1/length(y),length(y))
  obj = numeric(length=(maxiter+1))
  betalist = list(length(maxiter+1))
  betalist[[1]] = beta
  #beta0list = numeric(length(maxiter+1))
  #beta0 = sum(y-X%*%beta)/(length(y))
  #beta0list[1] = beta0
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      r = y - X[,-k]%*%beta[-k]
      beta[k] = (1/sum(omega*X[,k]^2))*
        soft_thresholding(t(omega*r)%*%X[,k] , lambda)
    }
    betalist[[j+1]] = beta
    obj[j] = (1/2)*(1/length(y))*norm(omega*(y - X%*%beta),'F')^2 + lambda*sum(abs(beta))
    if (norm(rbind(betalist[[j]]) - rbind(beta),'F') < tol) { break }
  }
  browser()
  return(list(obj=obj[1:j],beta=beta)) }

 chicago = read.table("http://freakonometrics.free.fr/chicago.txt",header=TRUE,sep=";")


      X = model.matrix(lm(Fire~.,data=chicago))[,2:4]
      for(j in 1:3) X[,j] = (X[,j]-mean(X[,j]))/sd(X[,j])
      y = chicago$Fire
      y = (y-mean(y))/sd(y)


      beta_init = lm(Fire~0+.,data=chicago)$coef
      lasso_coord_desc(X,y,beta_init,lambda=.001)


