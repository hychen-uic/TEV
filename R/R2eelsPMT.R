#' Permutation test for the least squares approach
#'
#' One paragraph describing this function
#'
#' @param y outcome
#' @param x covariates
#' @param npm permutation sample size. Default is 1000
#'
#' @details Details of this function...
#'
#' @return Ouput of this function...
#'
#' @references Reference 1 here...
#' @references Reference 2 here...
#'
#' @examples \dontrun{R2eelsPMT(y,x,npm=1000)}
#'
#' @export
R2eelsPMT=function(y,x,npm=1000){

  n=dim(x)[1]
  p=dim(x)[2]

  #1. Standardization

  for(j in 1:p){
    mu=mean(x[,j])
    sdx=sd(x[,j])
    x[,j]=(x[,j]-mu)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Compute the estimator
  xsvd=svd(x)
  r2=1-(n-1-sum((t(y)%*%xsvd$u)^2))/(n-p)
  r2=min(1,max(0,r2))

  #4. Compute permutation p-value for H0: r2=0.
  result=rep(0,npm)
  for(ii in 1:npm){
    samp=sample(n)
    yy=y[samp]

    result[ii]=1-(n-1-sum((t(yy)%*%xsvd$u)^2))/(n-p)
    result[ii]=min(1,max(result[ii],0))
  }

  pvalueEST=mean(1.0*(result>r2))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(pvalueEST,pvalueBOUND),r2,result)

}
