#' Permutation test for no explained variation
#'
#' This function performs test of no explained variation by
#' permuting the outcome using least square approach.
#'
#' @param y outcome
#' @param x covariates
#' @param npm permutation sample size. Default is 1000
#'
#' @details The computation permutes outcome and is computationally fast.
#'
#' @return Output is the estimate of the proportion of the explained variation
#' and the p-value of the test.
#'
#' @references Chen, H.Y.; Li, H.; Argos, M.; Persky, V.; Turyk, M.
#' Statistical methods for assessing explained variations of a health outcome by mixtures of exposures.
#' Prep. Spec. Issue Int. J. Environ. Res. Public Health 2022.
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
