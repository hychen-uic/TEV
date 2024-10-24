#' @import scalreg
#' @import lars
NULL
#'  Estimate the explained variation in high-dimensional linear model using the method CHIVE proposed
#'  by Cai and Guo (2020)
#'
#' This function estimates:
#'     (1). the proportion of the explained variation
#'     (2). the explained variation
#' by covariates in a linear model assuming the covariates have sparse effects.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param xext supplementary covariates of Nxp dimension.
#' @param alpha confidence level is 100*(1-alpha)%
#'
#' @details  Both point estimate and confidence intervals are computed.
#'
#' @return Estimate of the proportion of explained variation, variance estimates,
#' and confidence intervals.
#'
#' @references  Cai, T. T. Guo, Z. (2020).  Semisupervised inference for explained
#' variance in high dimensional linear regression and its applications,
#' Journal of Royal Statistical Society, Ser. B., 82, 391-419.
#'
#' @examples \dontrun{CHIVE(y,x,lam=0.05)}
#'
#'@export
CHIVE=function(y,x,xext=NULL,alpha=c(0.05)){

  n=dim(x)[1]
  p=dim(x)[2]
  if(!is.null(xext)){
    N=dim(xext)[1]
  }else{
    N=0
  }

  #1. Fit scaled lasso to estimate beta and sigma^2

  y=y-mean(y) # need to be centered.

  fit=scalreg(x,y)
  sigma2=fit[[1]]
  beta=fit[[2]]

  #2. Estimate \Sigma

  if(!is.null(xext)){
    xsig=(t(x)%*%x+t(xext)%*%xext)/(n+N)
  }else{
    xsig=t(x)%*%x/n
  }

  #3. Calibration
  pred=x%*%beta
  Q=as.numeric(t(beta)%*%xsig%*%beta)+2*sum((y-pred)*pred)/n
  vy=as.numeric(var(y))
  r2=min(1,max(0,Q/vy))

  #4. variance estimate
  evQ=4*sigma2*Q/n+sum((pred^2-Q)^2)/(n+N)^2
  evr2=evQ/vy^2+(Q/vy^2)^2*as.numeric(var(y^2))/n #This is not a consistent variance estimator

  #5 confidence interval
  len=length(alpha)
  cis2=Q+sqrt(evQ)*qnorm(c(alpha/2,1-alpha/2))
  cis2[c(1:len)]=cis2[c(1:len)]*(cis2[c(1:len)]>0)

  cir2=r2+sqrt(evr2)*qnorm(c(alpha/2,1-alpha/2))
  cir2[c(1:len)]=cir2[c(1:len)]*(cir2[c(1:len)]>0)
  cir2[len+c(1:len)]=cir2[len+c(1:len)]*(cir2[len+c(1:len)]<1)+1.0*(cir2[len+c(1:len)]>=1)

  #6. output result

  ind=rep(0,2*len)
  ind[2*c(1:len)-1]=c(1:len)
  ind[2*c(1:len)]=len+c(1:len)
  list(c(r2,evr2),cir2[ind],
       c(Q,evQ),cis2[ind])

  #[[1]]==> R2: estimate, estimated variance, estimated variance under normal error.
  #[[2]]==>confidence interval for R2
  #[[3]]==> sigma_s^2: estimate, estimated variance, estimated variance under normal error.
  #[[4]]==>confidence interval for sigma_s^2

}
