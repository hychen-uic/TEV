#' Estimation and inference on the explained variation
#' by the estimating equation approach with supplementary covariate data
#'
#' This function estimates the proportion of the explained variation using
#' the estimating equation approach incorporating supplementary covariate data.
#'
#' @param y outcome
#' @param x covariates
#' @param X supplementary covariate data
#' @param lam parameter adjusting the format of the weighting matrix. Default is 0.2
#' @param niter number of iterations for updating lambda. Default is 3
#' @param VA pre-calculated variance component
#' @param EB pre-calculated variance componenet
#' @param know if VA and EB are known, options include "yes" and "no". Default is "yes"
#' @param nrep Monte Carlo sample size for computing VA and EB
#'
#' @details The estimation approach does not assume independent covariates and can deal with the case n <= p.
#' This differs from \code{R2eesd} in the way the variance is estimated.
#'
#' @return Output include the estimator of the proportion of variation explained and
#' its variance, and the confidence intervals.
#'
#' @examples \dontrun{R2eesd0(y,x,X,lam=0.2,niter=3,VA=0,EB=0,know="no",nrep=1000)}
#'
#' @references Chen, H.Y.; Li, H.; Argos, M.; Persky, V.; Turyk, M.
#' Statistical methods for assessing explained variations of a health outcome by mixtures of exposures.
#' Prep. Spec. Issue Int. J. Environ. Res. Public Health 2022.
#'
#' @references An additional reference is to be added.
#'
#' @export
R2eesd0=function(y,x,X,lam=0.2,niter=3,VA=0,EB=0,know="yes",nrep=1000){

  n=dim(x)[1]
  p=dim(x)[2]
  N=dim(X)[1]
  if(dim(X)[2]!=p){
    # print("Stop: supplement data dimension does not match")
    # break
    stop("Supplement data dimension does not match!")
  }
  #1. Standardization

  XX=t(cbind(t(x),t(X))) #combine existing and supplement data on covariates
  for(j in 1:p){
    mu=mean(c(XX[,j]))
    sdx=sd(c(XX[,j]))
    X[,j]=(X[,j]-mu)/sdx
    x[,j]=(x[,j]-mu)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Sigular value decomposition
  M=x%*%solve(t(x)%*%x+t(X)%*%X)%*%t(x)*(n+N)/p

  r2=lam/(1+lam) #initial value
  for(ii in 1:niter){
    if(ii>1 & r2<1){lam=r2/(1-r2)}

    IM=solve(diag(rep(1,n))+lam*M)
    W=IM%*%(M-diag(rep(1,n)))%*%IM
    #3. Compute the estimators
    den=sum(diag(W%*%(M-diag(rep(1,n)))))
    num=t(y)%*%W%*%y-sum(diag(W))

    r2=as.numeric(num/den)
    r2=min(1,max(0,r2))
  }

  #4. Estimate the variance
  u=c(rep(1,p/2)/sqrt(p/2),rep(0,p/2))
  if(know!='yes'){ # calculating VA and EB by simulation
    SA=rep(0,nrep)
    SB=rep(0,nrep)
    for(j in 1:nrep){
      z=matrix(rnorm(n*p),ncol=p)
      zu=z%*%u
      Z=matrix(rnorm(N*p),ncol=p)
      SM=z%*%solve(t(z)%*%z+t(Z)%*%Z)%*%t(z)*(n+N)/p
      ISM=solve(diag(rep(1,n))+lam*SM)
      SW=ISM%*%(SM-diag(rep(1,n)))%*%ISM
      SWzu=SW%*%zu
      SA[j]=as.numeric(t(zu)%*%SW%*%zu)
      SB[j]=t(SWzu)%*%SWzu
    }
    VA=var(SA)
    EB=mean(SB)
  }

  # Variance under normal random error
  S=sum(diag(W%*%W))
  vest=(r2^2*VA+4*r2*(1-r2)*EB+2*(1-r2)^2*S)/den^2

  # Variance without normal random error assumption
  T=sum(diag(W)^2)
  veps1=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2*r2^2
  vest1=vest+T*(max(veps1,0)-2*(1-r2)^2)/den^2

  #veps1=mean(y^4)-6*r2+3*r2^2-(1-r2)^2
  #vest1=vest+(T*(max(0,veps1)-2*(1-r2)^2))/den^2

  ci=r2+sqrt(vest)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci[2*c(1:3)-1]=ci[2*c(1:3)-1]*(ci[2*c(1:3)-1]>0)
  ci[2*c(1:3)]=ci[2*c(1:3)]*(ci[2*c(1:3)]<1)+1.0*(ci[2*c(1:3)]>=1)

  ci1=r2+sqrt(vest1)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci1[2*c(1:3)-1]=ci1[2*c(1:3)-1]*(ci1[2*c(1:3)-1]>0)
  ci1[2*c(1:3)]=ci1[2*c(1:3)]*(ci1[2*c(1:3)]<1)+1.0*(ci1[2*c(1:3)]>=1)

  #5. output result

  list(c(r2,vest,vest1),ci,ci1)
  #[[1]]==> estimate, estimated variance (under normal), estimated variance.
  #[[2]]==>confidence interval under normal
  #[[3]]==>confidence interval without normal
}
