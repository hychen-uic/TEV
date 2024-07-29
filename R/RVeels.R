#' Least-square approach to proportion of the explained variation
#'
#' This function uses least-square estimates in computing the proportion of the explained variation.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param lam parameter for altering the weighting matrix.
#' @param alpha a vector of type I errors used to generate (1-alpha)confidence intervals.
#'
#'
#' @details This method works only for the case n>p. It uses the least-square approach
#' for the estimation. Covariates are allowed to be correlated.
#'
#' @return Estimate of proportion of the explained variation, variance estimates under normality
#'         and non-normality error, and the corresponding confidence intervals.
#'
#' @references Chen, H.Y. (2022). Statistical inference on explained variation in high-dimensional
#'  linear model with dense effects. arXiv:2201.08723
#' @references Chen, H. Y., Li, H., Argos, M., Persky, V. W., and Turyk, M. (2022). Statistical Methods
#' for Assessing Explained Variation of a Health Outcome by Mixture of Exposures. International Journal
#' of Environmental Research and Public Health.
#'
#' @examples \dontrun{R2eels(y,x)}
#'
#' @export
RVeels=function(y,x,lam=0.1, alpha=0.05 ){

  # y==outcome
  # x==covariates

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
  Xsvd=svd(x,nv=0)
  r2=1-(n-1-sum((t(y)%*%Xsvd$u)^2))/(n-p)
  r2=min(1,max(0,r2))

  vy=sdy*sdy
  s2=vy*r2   # the explained variation

  #3. estimate variance

  Mev=(Xsvd$d>0)*n/p
  Wev=(Mev-1)/(1+lam*Mev)^2
  Dev=Wev*Mev
  D1=sum(Dev)/n
  if(n>=p){
    D2=(sum(Wev)-n+p)/n
    A0=2*p*var(Dev)
    A01=2*p*cov(Dev,Mev)
    A1=2*p*var(Mev)
    B=sum(Wev^2*Mev)
    S=sum(Wev^2)+n-p

    W=Xsvd$u%*%diag(Wev+1)%*%t(Xsvd$u)-diag(rep(1,n))
    T=sum(diag(W^2))

    den=sum(Wev*(Mev-1))+n-p
  }else{
    D2=sum(Wev)/n
    A0=2*n*(mean(Dev^2)-(mean(Dev))^2*n/p)
    A01=2*n*(mean(Dev*Mev)-mean(Dev)*mean(Mev)*n/p)
    A1=2*n*(mean(Mev^2)-(mean(Mev))^2*n/p)
    B=sum(Wev^2*Mev)
    S=sum(Wev^2)

    W=Xsvd$u%*%diag(Wev)%*%t(Xsvd$u)
    T=sum(diag(W^2))
    den=sum(Wev*(Mev-1))

  }
  DELTA=r2*D1+(1-r2)*D2
  vestr2n=r2^2*(A0-2*DELTA*A01+DELTA^2*A1)
  vestr2n=vestr2n+4*r2*(1-r2)*(B-2*DELTA*D1*n+n*DELTA^2)
  vestr2n=vestr2n+2*(1-r2)^2*(S-2*DELTA*D2*n+n*DELTA^2)

  vests2n=r2^2*(A0-2*D2*A01+D2^2*A1)
  vests2n=vests2n+4*r2*(1-r2)*(B-2*n*D1*D2+n*D2^2)
  vests2n=vests2n+2*(1-r2)^2*(S-D2^2*n)
  vests2n=vests2n*vy^2

  #4(b). Without normal random error assumption

  M=Xsvd$u%*%diag(Mev)%*%t(Xsvd$u)
  veps2=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2*r2^2

  vestr2=vestr2n+(T-2*DELTA*D2*n+n*DELTA^2)*(max(veps2,0)-2*(1-r2)^2)
  vests2=vests2n+vy^2*(T-D2^2*n)*(max(veps2,0)-2*(1-r2)^2)

# proportion of explained variation
  len=length(alpha)
  vestr2=vestr2/den^2

  cir2=r2+sqrt(vestr2)*qnorm(c(alpha/2,1-alpha/2))
    cir2[c(1:len)]=cir2[c(1:len)]*(cir2[c(1:len)]>0)
    cir2[len+c(1:len)]=cir2[len+c(1:len)]*(cir2[len+c(1:len)]<1)+1.0*(cir2[len+c(1:len)]>=1)
  vestr2n=vestr2n/den^2  #under normal error
  cir2n=r2+sqrt(vestr2n)*qnorm(c(alpha/2,1-alpha/2))
    cir2n[c(1:len)]=cir2n[c(1:len)]*(cir2n[c(1:len)]>0)
    cir2n[len+c(1:len)]=cir2n[len+c(1:len)]*(cir2n[len+c(1:len)]<1)+1.0*(cir2n[len+c(1:len)]>=1)

# explained variation
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

}
