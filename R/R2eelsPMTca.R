#' Conditional permutation test for no explained variation using the least-square approach
#'
#' This function performs conditional tests of no variation explained by the additional
#' covariates when adjusted for the one set of covariates.
#'
#' @param y outcome
#' @param x all the covariates
#' @param pa dimension of the first part of covariates to be adjusted
#' @param npm permutation sample size. Default is 1000
#'
#' @details The conditional permutation test needs to perform the second part of the
#' covariates and can be computationally slow.
#'
#' @return Output includes the p-values (estimate and bound) of the test, estimate of
#' proportion of the extra explained variation, and simulation results.
#'
#' @references Chen, H.Y.; Li, H.; Argos, M.; Persky, V.; Turyk, M.
#' Statistical methods for assessing explained variations of a health outcome by mixtures of exposures.
#' Under review for Prep. Spec. Issue Int. J. Environ. Res. Public Health 2022.
#'
#' @examples \dontrun{R2eelsPMTca(y, x, pa, npm = 1000)}
#'
#' @export
R2eelsPMTca=function(y, x, pa, npm = 1000){

  n=dim(x)[1]
  p=dim(x)[2]

  #1. Standardization

  for(j in 1:p){
    mux=mean(x[,j])
    sdx=sd(x[,j])
    x[,j]=(x[,j]-mux)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Compute the estimato
  xasvd=svd(x[,1:pa])
  r2a=1-(n-1-sum((t(y)%*%xasvd$u)^2))/(n-pa)
  r2a=min(1,max(0,r2a))

  xsvd=svd(x)
  r2=1-(n-1-sum((t(y)%*%xsvd$u)^2))/(n-p)
  r2=min(1,max(0,r2))

  r2ba=max(0,r2-r2a)

  print(c(r2,r2a,r2ba))
  #4. Compute permutation p-value for H0: r2=0.
  result=rep(0,npm)
  for(jj in 1:npm){
    if(jj-as.integer(jj/50)*50==1){print(c(jj,jj))}

    samp=sample(n)
    xx=cbind(x[,1:pa],x[samp,(pa+1):p])

    svdxx=svd(xx)
    r2p=1-(n-1-sum((t(y)%*%svdxx$u)^2))/(n-p)

    r2p=min(1,max(r2p,0))

    result[jj]=max(0,r2p-r2a)
  }

  pvalueEST=mean(1.0*(result>r2ba))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(pvalueEST,pvalueBOUND),c(r2,r2a,r2ba),result)

}
