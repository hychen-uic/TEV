#' Estimating equation approach to the proportion of the explained variation with supplementary covariate data
#'
#' This function computes:
#'      (1). the proportion of the explained variation
#'      (2). The explained variation
#' adjusting for the correlation in covariates with possible suplementary data
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param X supplementary covariates: a matrix of Nxp dimension.
#' @param lam parameter for altering the weighting matrix.
#' @param niter number of iterations for updating lam.
#' @param KV, the first component of the vector KV=kappa_1, the second=kappa_2, the third=kappa_3
#' @param know if KV is known, options include "yes" and "no". Default is "no".
#' @param nrep Monte Carlo sample size for computinging KV.
#' @param alpha a vector of type I errors used to generate (1-alpha)confidence intervals.
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
#' @references Chen H. Y, Zhang, B. and Pan, G.(2023).Estimation and inference on explained variation
#'with possible supplementary data, Manuscript.
#'
#' @examples \dontrun{R2eesd(y,x,X,lam=0.2,niter=3,alpha=c(0.1,0.05,0.01),KV=rep(0,3),know="no",nrep=1000)}
#'
#'
#' @export
RVeesd=function(y,x,X,lam=0.2,niter=1,alpha=c(0.1,0.05,0.01),KV=rep(0,3),know="no",nrep=1000){

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
    mu=mean(XX[,j])
    sdx=sd(XX[,j])
    X[,j]=(X[,j]-mu)/sdx
    x[,j]=(x[,j]-mu)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Sigular value decomposition
  x=as.matrix(x)
  X=as.matrix(X)
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
  vy=sdy*sdy
  s2=vy*r2  #explained variation

  #4. Estimate the variance
  if(know!='yes'){ # calculating VA and EB by simulation
    u=rep(1,p)/sqrt(p)
    SUZWZU=rep(0,nrep)
    SUZWWZU=rep(0,nrep)
    SUZZU=rep(0,nrep)
    STRWM=rep(0,nrep)
    #STRW=rep(0,nrep)

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
      #STRW[j]=sum(diag(SW))/n
      STRWM[j]=sum(diag(SW%*%SM))/n
      }
    K1=var(SUZWZU-STRWM)
    K2=cov(SUZWZU-STRWM,SUZZU-n)
    K3=var(SUZZU-n)
  }else{
    K1=KV[1];K2=KV[2];K3=KV[3]
  }

  # Variance under normal random error

  D1=sum(diag(W%*%M))/n
  D2=sum(diag(W))/n
  DELTA=r2*D1+(1-r2)*D2

  W2=W%*%W
  B=sum(diag(W2%*%M))
  S=sum(diag(W2))
  T=sum(diag(W)^2)

  vestr2n=r2^2*(K1-2*DELTA*K2+DELTA^2*K3)
  vestr2n=vestr2n+4*r2*(1-r2)*(B-2*n*D1*DELTA+n*DELTA^2)
  vestr2n=vestr2n+2*(1-r2)^2*(S-2*DELTA*D2*n+n*DELTA^2)

  vests2n=r2^2*(K1-2*D2*K2+D2^2*K3)
  vests2n=vests2n+4*r2*(1-r2)*(B-2*n*D1*D2+n*D2^2)
  vests2n=vests2n+2*(1-r2)^2*(S-D2^2*n)
  vests2n=vests2n*vy^2

  # Variance without normality assumption

  veps2=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2*r2^2
  veps2=max(0,veps2)

  vestr2=vestr2n+(T-2*DELTA*D2*n+n*DELTA^2)*(veps2-2*(1-r2)^2)
  vests2=vests2n+vy^2*(T-D2^2*n)*(max(veps2,0)-2*(1-r2)^2)

## proportion of explained variation
  len=length(alpha)
  vestr2=vestr2/den^2
  cir2=r2+sqrt(vestr2)*qnorm(c(alpha/2,1-alpha/2))
    cir2[c(1:len)]=cir2[c(1:len)]*(cir2[c(1:len)]>0)
    cir2[len+c(1:len)]=cir2[len+c(1:len)]*(cir2[len+c(1:len)]<1)+1.0*(cir2[len+c(1:len)]>=1)
  vestr2n=vestr2n/den^2  #under normal error
  cir2n=r2+sqrt(vestr2n)*qnorm(c(alpha/2,1-alpha/2))
    cir2n[c(1:len)]=cir2n[c(1:len)]*(cir2n[c(1:len)]>0)
    cir2n[len+c(1:len)]=cir2n[len+c(1:len)]*(cir2n[len+c(1:len)]<1)+1.0*(cir2n[len+c(1:len)]>=1)

## explained variation
  vests2=vests2/den^2
  cis2=s2+sqrt(vests2)*qnorm(c(alpha/2,1-alpha/2))
    cis2[c(1:len)]=cis2[c(1:len)]*(cis2[c(1:len)]>0)
  vests2n=vests2n/den^2  #under normal error
  cis2n=s2+sqrt(vests2n)*qnorm(c(alpha/2,1-alpha/2))
    cis2n[c(1:len)]=cis2n[c(1:len)]*(cis2n[c(1:len)]>0)

  #5. output result

  len=length(alpha)
  ind=rep(0,2*len)
  ind[2*c(1:len)-1]=c(1:len)
  ind[2*c(1:len)]=len+c(1:len)
  list(c(r2,vestr2,vestr2n),cir2[ind],cir2n[ind],
       c(s2,vests2,vests2n),cis2[ind],cis2n[ind])


  #[[1]]==> R2: estimate, estimated variance, estimated variance under normal error.
  #[[2]]==>confidence interval for R2
  #[[3]]==>confidence interval for R2 under normal error.
  #[[4]]==> sigma_s^2: estimate, estimated variance, estimated variance under normal error.
  #[[5]]==>confidence interval for sigma_s^2
  #[[6]]==>confidence interval for sigma_s^2 under normal error.

}
