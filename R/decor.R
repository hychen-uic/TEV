#' Decorrelate the data matrix to be used in the heritability estimation
#'
#' This function decorrelate the data matrix to make it approximately uncorrelated.
#'
#' @param dat  is a nxp data matrix
#' @param decorate the block size in term of a fraction of sample size n.
#' @param inter  0 means do block wise once, inter=1 means do another interlanced decorrelation
#'
#' @details It does not require to know the correlation structure and can handle $n<=p$.
#' The idea is to estimate the inverse of the blockwised covariance matrix rather than
#' the inverse of the covariance matrix which may not be feasible when n<=p.
#'
#' @return  precision matrix estimate
#'
#' @references To be added.
#'
decor=function(dat,decorate=0.5,inter=0){
  #drate<1, m=n*drate used for sub-block inversion

  n=dim(dat)[1]
  p=dim(dat)[2]
  m=as.integer(n*decorate)
  decorx=matrix(0, nrow=n,ncol=p)

  if(p<=m){
    decorx=dat%*%sqrtpdm(solve(cov(dat)))[[1]]
  }else{
    KN=as.integer(p/m)
    for(k in 1:KN){
      start=(k-1)*m+1
      end=k*m
      decorx[,start:end]=dat[,start:end]%*%sqrtpdm(solve(cov(dat[,start:end])))[[1]]
    }
    start=m*KN+1
    if(p>start){
      decorx[,start:p]=dat[,start:p]%*%sqrtpdm(solve(cov(dat[,start:p])))[[1]]
    }

    if(inter==1){ #interlance decorrelation
      hm=as.integer(m/2)
      if(KN>1){
        for(k in 1:(KN-1)){
          start=(k-1)*m+hm+1
          end=k*m+hm
          decorx[,start:end]=decorx[,start:end]%*%sqrtpdm(solve(cov(decorx[,start:end])))[[1]]
        }
        start=m*KN+hm+1
        if(p>start){
          decorx[,start:p]=decorx[,start:p]%*%sqrtpdm(solve(cov(decorx[,start:p])))[[1]]
        }
      }else{
        start=hm+1
        if(p>start){
          decorx[,start:p]=decorx[,start:p]%*%sqrtpdm(solve(cov(decorx[,start:p])))[[1]]
        }
      }
    } #end interlance decorrelation
  }

  return(decorx)
}
