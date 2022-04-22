#' @import stats
NULL
#' Estimating equation approach to the proportion of the explained variation
#'
#' This approach estimates the proportion of the explained variation in a linear model
#' assuming the covariates are independent.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param lam parameter for altering the weighting matrix.
#' @param niter number of iterations for updating the lam parameter.
#'
#' @details  Both point estimate and confidence intervals
#' are computed. Two set of confidence intervals under normal or non-normal error are computed.
#'
#' @return Estimate of the proportion of explained variation, variance estimates under normality
#' and non-normality assumptions, and confidence intervals under normality and non-normality assumptions.
#'
#' @references Chen, H.Y. (2022). Statistical inference on explained variation in high-dimensional
#'  linear model with dense effects. arXiv:2201.08723.
#' @references Chen, H. Y., Li, H., Argos, M., Persky, V. W., and Turyk, M. (2022). Statistical Methods
#' for Assessing Explained Variation of a Health Outcome by Mixture of Exposures. International Journal
#' of Environmental Research and Public Health.
#'
#' @examples \dontrun{R2ee(y,x,lam=1.0,niter=1)}
#'
#'@export
R2ee=function(y, x, lam = 1.0, niter = 1){

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

  #2. Sigular value decomposition
  Xsvd=svd(x,nv=0) #nv=0 means not computing v matrix
  # singular value decomposition
  # $u%*%diag($d)%*%t($v)=X, t($u)%*%$u=I, t($v)%*%$v=I
  Mev=Xsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.

  r2=lam/(1+lam) # initial value
  for(ii in 1:niter){ #iteration to update lambda
    if(ii>1 & r2<1){lam=r2/(1-r2)}
    Wev=(Mev-1)/(1+lam*Mev)^2  #vector of eigenvalues of weight matrix
    #3. Compute the estimators
    uy=t(Xsvd$u)%*%y
    u1=t(Xsvd$u)%*%rep(1,n)
    if(n>=p){
      com=sum(u1^2*(Wev+1))/n-1
      num=sum(uy^2*(Wev+1))-sum(y^2)-sum(Wev)+n-p
      den=sum(Wev*(Mev-1))+n-p
    }else{
      com=sum(u1^2*Wev)/n
      num=sum(uy^2*Wev)-sum(Wev)
      den=sum(Wev*(Mev-1))
    }

    r2=(num+com)/(den+com)
    r2=min(1,max(0,r2))
  }
  #4. estimate variance
  #4(a). With normal random error assumption

  WMev=Wev*Mev
  D1=sum(WMev)/n
  if(n>=p){
    W=Xsvd$u%*%diag(Wev+1)%*%t(Xsvd$u)-diag(rep(1,n))
    D2=sum(Wev)/n-1+p/n
    S=sum(Wev^2)+n-p
  }else{
    W=Xsvd$u%*%diag(Wev)%*%t(Xsvd$u)
    D2=sum(Wev)/n
    S=sum(Wev^2)
  }
  tau1=sum(WMev^2)/p-(sum(WMev))^2/p^2
  tau01=sum(WMev*Mev)/p-sum(WMev)*sum(Mev)/p^2
  tau0=sum(Mev^2)/p-(sum(Mev))^2/p^2
  A=2*p*(tau1-2*(r2*D1+(1-r2)*D2)*tau01+(r2*D1+(1-r2)*D2)^2*tau0)

  B=sum((Wev-(r2*D1+(1-r2)*D2))^2*Mev)
  S=S-n*D2^2+n*r2^2*(D1-D2)^2

  vest=r2^2*A+4*r2*(1-r2)*B+2*(1-r2)^2*S

  #4(b). Without normal random error assumption

  M=Xsvd$u%*%diag(Mev)%*%t(Xsvd$u)
  T=sum(diag(W)^2)-n*D2^2+n*r2^2*(D1-D2)^2
  veps2=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2*r2^2
  vest1=vest+(max(veps2,0)-2*(1-r2)^2)*T

  #print(c((max(veps2,0)-2*(1-r2)^2)*T/den^2,den,T,veps2))

  vest1=vest1/den^2 #(den+com)^2
  vest=vest/den^2 #(den+com)^2

  ci=r2+sqrt(vest)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci1=r2+sqrt(vest1)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci[2*c(1:3)-1]=ci[2*c(1:3)-1]*(ci[2*c(1:3)-1]>0)
  ci[2*c(1:3)]=ci[2*c(1:3)]*(ci[2*c(1:3)]<1)+1.0*(ci[2*c(1:3)]>=1)
  ci1[2*c(1:3)-1]=ci1[2*c(1:3)-1]*(ci1[2*c(1:3)-1]>0)
  ci1[2*c(1:3)]=ci1[2*c(1:3)]*(ci1[2*c(1:3)]<1)+1.0*(ci1[2*c(1:3)]>=1)

  #5. output result

  list(c(r2,vest,vest1),ci,ci1)
  #[[1]]==> estimate, estimated variance (under normal), estimated variance.
  #[[2]]==>confidence interval under normal
  #[[3]]==>confidence interval without normal
}
