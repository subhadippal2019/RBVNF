#' rMatrixNormal(M,V1,V2)
#' This function generates single sample from a matrix normal random variable 
#' @param $M$ the matrix which is the mean of the matrix Normal random variable
#' @param $V1$ and $V2$ are variance component such that the variance of Vec(X) is V2 cross V1 
#' @export
rMatrixNormal<-function(M,V1,V2){
  #M is the mean matrix, and the V1, V2 are symmetric positive definite matrices
  
  nrow=nrow(M)
  ncol=ncol(M)
  mu=as.vector(M)
  V=kronecker(V2 , V1)
  y=mvrnorm(n=1,mu,V)
  Y=matrix(y,nrow,ncol)
  return(Y)
}
