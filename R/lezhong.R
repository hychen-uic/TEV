#' Estimate the precision matrix when the network structure is known
#'
#' This function find the precision matrix estimate by the method of Le and Zhong (2022)
#'
#' @param dat  is a nxp data matrix
#' @param network is the fixed sparse struct of the precision matrix given
#'
#' @details This approach works for sparse precision matrix
#'
#' @return  precision matrix estimate
#'
#' @references  Le, T. M. and Zhong, P. S. (2022) High-dimensional precision matrix estimation with a known graphical structure. Stat, 11, e424.
#'

lezhong=function(dat,network){
  n=dim(dat)[1]
  p=dim(dat)[2]
  diag(network)=rep(1,p) #make sure diagonal elements in

  omega=matrix(0,nrow=p,ncol=p)
  dat=dat-rep(1,n)%*%t(apply(dat,c(2),mean))
  S=t(dat)%*%dat/(n-1)
  for(i in 1:p){
    if(sum(network[,i])!=0){
      B=matrix(0,nrow=p, ncol=sum(network[,i]))
      m=0
      for(j in 1:p){
        if(as.integer(network[i,j])!=0){
          m=m+1
          B[j,m]=1
        }
      } #end of assignment of B
      omega[,i]=B%*%solve(t(B)%*%S%*%B)%*%B[i,]
    }
  }
  omega=(omega+t(omega))/2 # symmetrize

  return(omega)
}
