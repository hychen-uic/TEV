#' confidence interval for estimating equation approach with supplementary covariate data
#'
#' One paragraph describing this function...
#'
#' @param y outcome
#' @param x covariates
#' @param X supplementary covariate data
#' @param lam parameter adjusting the format of the weighting matrix. Default is 0.2
#' @param niter number of iterations for updating lambda. Default is 3
#' @param V V[1:3] == V[1]=VA, V[2]=VAB, V[3]=VB
#' @param E E[1:3] == E[1]=EC, E[2]=ED, E[3]=EF
#' @param know if VA and EB are known, options include "yes" and "no". Default is "yes"
#' @param nrep Monte Carlo sample size for updating VA and EB
#'
#' @details Details of this function...
#'
#' @return Output of this function...
#'
#' @references reference 1 here...
#' @references reference 2 here...
#'
#' @examples \dontrun{R2eesd(y,x,X,lam=0.2,niter=3,V=rep(0,3),E=rep(0,3),know="yes",nrep=1000)}
#'
#'
#' @export
R2eesd=function(y,x,X,lam=0.2,niter=3,V=rep(0,3),E=rep(0,3),know="yes",nrep=1000){

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
  M=x%*%chol2inv(chol(t(x)%*%x+t(X)%*%X))%*%t(x)*(n+N)/p

  r2=lam/(1+lam) #initial value
  for(ii in 1:niter){
    if(ii>1 & r2<1){lam=r2/(1-r2)}

    IM=chol2inv(chol(diag(rep(1,n))+lam*M))
    W=IM%*%(M-diag(rep(1,n)))%*%IM
    #3. Compute the estimators
    den=sum(diag(W%*%(M-diag(rep(1,n)))))
    num=t(y)%*%W%*%y-sum(diag(W))

    r2=as.numeric(num/den)
    r2=min(1,max(0,r2))
  }

  #4. Estimate the variance
  if(know!='yes'){ # calculating VA and EB by simulation
    u=rep(1,p)/sqrt(p)
    SUZWZU=rep(0,nrep)
    SUZWWZU=rep(0,nrep)
    SUZZU=rep(0,nrep)
    STRW=rep(0,nrep)
    STRWM=rep(0,nrep)
    for(j in 1:nrep){
      z=matrix(rnorm(n*p),ncol=p)
      z=zscale(z)[[1]]
      zu=z%*%u
      Z=matrix(rnorm(N*p),ncol=p)
      Z=zscale(Z)[[1]]
      SM=z%*%chol2inv(chol(t(z)%*%z+t(Z)%*%Z))%*%t(z)*(n+N)/p
      ISM=chol2inv(chol(diag(rep(1,n))+lam*SM))
      SW=ISM%*%(SM-diag(rep(1,n)))%*%ISM
      SWzu=SW%*%zu

      SUZZU[j]=sum(zu^2)
      SUZWZU[j]=t(zu)%*%SW%*%zu
      SUZWWZU[j]=t(SWzu)%*%SWzu
      STRW[j]=sum(diag(SW))/n
      STRWM[j]=sum(diag(SW%*%SM))/n
    }
    VA=var(SUZWZU-SUZZU*STRWM)
    VB=var(SUZWZU-SUZZU*STRW)
    VAB=var(SUZWZU-SUZZU*STRWM,SUZWZU-SUZZU*STRW)
    EC=mean(SUZWWZU)+mean(SUZZU)*mean(STRW)^2-mean(SUZWZU)*mean(STRW)
    ED=2*mean(STRW)*mean(STRWM-STRW)*mean(SUZZU)-mean(STRWM-STRW)*mean(SUZWZU)
    EF=mean(STRWM-STRW)^2*mean(SUZZU)
  }else{
    VA=V[1];VAB=V[2];VB=V[3]
    EC=E[1];ED=E[2];EF=E[3]
  }

  # Variance under normal random error
  #1# are original variance estimates
  #1# S=sum(diag(W%*%W))
  #1# vest=(r2^2*VA+4*r2*(1-r2)*EB+2*(1-r2)^2*S)

  S0=sum(diag(W))/n
  S1=sum(diag(W%*%M))/n
  S2=sum(diag(W%*%W))
  T0=sum(diag(W)^2)
  vest=r2^2*(r2^2*VA+2*r2*(1-r2)*VAB+(1-r2)^2*VB)
  vest=vest+4*r2*(1-r2)*(EC+ED*r2+EF*r2^2)
  vest=vest+2*(1-r2)^2*(S2+n*(r2*S1-(1-r2)*S0)^2+2*n*(r2*S1-(1-r2)*S0)*S0)

  # Variance without normal random error assumption
  #1# T=sum(diag(W)^2)
  #1# veps2=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2*r2^2
  #1# veps2=max(0,veps2)
  #1# vest1=vest+(T*(veps2-2*(1-r2)^2))
  veps2=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2*r2^2
  veps2=max(0,veps2)
  vest1=vest+(veps2-2*(1-r2)^2)*(T0+n*(r2*S1-(1-r2)*S0)^2+2*n*(r2*S1-(1-r2)*S0)*S0)

  vest=vest/den^2
  ci=r2+sqrt(vest)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci[2*c(1:3)-1]=ci[2*c(1:3)-1]*(ci[2*c(1:3)-1]>0)
  ci[2*c(1:3)]=ci[2*c(1:3)]*(ci[2*c(1:3)]<1)+1.0*(ci[2*c(1:3)]>=1)

  vest1=vest1/den^2
  ci1=r2+sqrt(vest1)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci1[2*c(1:3)-1]=ci1[2*c(1:3)-1]*(ci1[2*c(1:3)-1]>0)
  ci1[2*c(1:3)]=ci1[2*c(1:3)]*(ci1[2*c(1:3)]<1)+1.0*(ci1[2*c(1:3)]>=1)

  #5. output result

  list(c(r2,vest,vest1),ci,ci1)
  #[[1]]==> estimate, estimated variance (under normal), estimated variance.
  #[[2]]==>confidence interval under normal
  #[[3]]==>confidence interval without normal
}
