#' Estimating equation approach to the proportion of the explained variation with supplementary covariate data
#'
#' This function computes the proportion of the explained variation adjusting for the correlation in covariates.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param X supplementary covariates: a matrix of Nxp dimension.
#' @param lam parameter for altering the weighting matrix.
#' @param niter number of iterations for updating lam.
#' @param VA a 3-dimensional vector for V1,V2,V3
#' @param EB a 3-dimensional vector for D3,D4,D5
#' @param know if (VA, EB) are known, options include "yes" and "no". Default is "yes" (usually for simulation only).
#' @param nrep Monte Carlo sample size for computing VA and EB.
#'
#' @details The estimation approach does not assume independent covariates and can deal
#' with the case \eqn{n\le p}. But require the sample sizes of x and X combined be greater than p.
#'
#' @return Estimate of the proportion of explained variation, variance estimates under normality
#' and non-normality assumptions, and confidence intervals under normality and non-normality assumptions.
#'
#' @references Chen, H. Y., Li, H., Argos, M., Persky, V. W., and Turyk, M. (2022). Statistical Methods
#' for Assessing Explained Variation of a Health Outcome by Mixture of Exposures. International Journal
#' of Environmental Research and Public Health.
#' @references reference 2 to be added.
#'
#' @examples \dontrun{R2eesd(y,x,X,lam=1,niter=1,V=rep(0,3),B=rep(0,3),know="no",nrep=1000)}
#'
#' @export
R2eesd=function(y, x, X, lam =1, niter = 1, VA= rep(0, 3),EB=rep(0,3), know = "yes", nrep = 1000){

  n=dim(x)[1]
  p=dim(x)[2]
  N=dim(X)[1]
  if(dim(X)[2]!=p){
    # print("Stop: supplement data dimension does not match")
    # break
    stop("Supplement data dimension does not match!")
  }
  #1. Standardization

  XX=rbind(x,X) #combine existing and supplement data on covariates
  for(j in 1:p){
    XX[,j]=(XX[,j]-mean(c(XX[,j])))/sd(c(XX[,j]))
    x[,j]=XX[1:n,j] #(x[,j]-mean(c(x[,j])))/sd(c(x[,j]))
  }
  y=(y-mean(y))/sd(y)

  #2. compute the estimate
  M=x%*%chol2inv(chol(t(XX)%*%XX/(n+N)))%*%t(x)/p

  r2=lam/(1+lam) #initial value
  for(ii in 1:niter){
    if(ii>1 & r2<1){lam=r2/(1-r2)}

    IM=chol2inv(chol(diag(rep(1,n))+lam*M))
    W=IM%*%(M-diag(rep(1,n)))%*%IM
    #3. Compute the estimators
    D1=sum(diag(W%*%M))
    D2=sum(diag(W))
    den=D1-D2
    num=t(y)%*%W%*%y-D2

    r2=as.numeric(num/den)
    r2=min(1,max(0,r2))
  }

  #4. Estimate the variance
  if(know!='yes'){ # calculating VA and EB by simulation
    u=rep(1,p)/sqrt(p)
    SUZWZU=rep(0,nrep)
    SUZWWZU=rep(0,nrep)
    SUZZU=rep(0,nrep)
    STRWM=rep(0,nrep)
    for(j in 1:nrep){
      z=matrix(rnorm(n*p),ncol=p)
      Z=matrix(rnorm(N*p),ncol=p)
      ZZ=zscale(rbind(z,Z))[[1]]
      z=zscale(z)[[1]] #ZZ[1:n,]
      zu=z%*%u

      SM=z%*%chol2inv(chol(t(ZZ)%*%ZZ/(n+N)))%*%t(z)/p
      ISM=chol2inv(chol(diag(rep(1,n))+lam*SM))
      SW=ISM%*%(SM-diag(rep(1,n)))%*%ISM
      SWzu=SW%*%zu

      SUZZU[j]=sum(zu^2)
      SUZWZU[j]=t(zu)%*%SW%*%zu
      SUZWWZU[j]=t(SWzu)%*%SWzu
      STRWM[j]=sum(diag(SW%*%SM))
    }
    B1=mean(SUZWWZU)/n  # keep these with known limits to have positive var.
    B2=mean(SUZWZU)/n
    B3=mean(SUZZU)/n
    A1=mean((SUZWZU-STRWM)^2)/n
    A2=mean((SUZWZU-STRWM)*(SUZZU-n))/n
    A3=mean((SUZZU-n)^2)/n

    EB=c(B1,B2,B3)
    VA=c(A1,A2,A3)
  }else{
    B1=EB[1];B2=EB[2];B3=EB[3]
    A1=VA[1];A2=VA[2];A3=VA[3]
  }

  # Variance under normal random error
  #1# are original variance estimates
  D1=D1/n; D2=D2/n
  A=A1-2*(r2*D1+(1-r2)*D2)*A2+(r2*D1+(1-r2)*D2)^2*A3
  B=B1-2*(r2*D1+(1-r2)*D2)*B2+(r2*D1+(1-r2)*D2)^2*B3
  ST=r2^2*(D1-D2)^2-D2^2
  S=sum(diag(W%*%W))/n+ST
  T=sum(diag(W)^2)/n+ST

  vest=r2^2*A+4*r2*(1-r2)*B+2*(1-r2)^2*S

  # Variance without normal random error assumption
  #1# veps2=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2*r2^2
  #1# veps2=max(0,veps2)
  veps2=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2*r2^2
  veps2=max(0,veps2)
  vest1=vest+(veps2-2*(1-r2)^2)*T

  vest=n*vest/den^2
  ci=r2+sqrt(vest)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci[2*c(1:3)-1]=ci[2*c(1:3)-1]*(ci[2*c(1:3)-1]>0)
  ci[2*c(1:3)]=ci[2*c(1:3)]*(ci[2*c(1:3)]<1)+1.0*(ci[2*c(1:3)]>=1)

  vest1=n*vest1/den^2
  ci1=r2+sqrt(vest1)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci1[2*c(1:3)-1]=ci1[2*c(1:3)-1]*(ci1[2*c(1:3)-1]>0)
  ci1[2*c(1:3)]=ci1[2*c(1:3)]*(ci1[2*c(1:3)]<1)+1.0*(ci1[2*c(1:3)]>=1)

  #5. output result

  list(c(r2,vest,vest1),ci,ci1)
  #[[1]]==> estimate, estimated variance (under normal), estimated variance.
  #[[2]]==>confidence interval under normal
  #[[3]]==>confidence interval without normal
}
