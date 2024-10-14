#' This function implement the maximum likelihood approach with variance estimated
#' by inverse of information matrix.
#'
#' This the reml for the high dimensional linear model
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param alpha a vector of type I errors used to generate (1-alpha)confidence intervals.
#' @param lam initial value
#' @param niter the number of iterations for finding the signal noise ratio
#' @param eps the convergence criterion for the iteration
#'
#' @details  This method assume the independent covariates with fixed effects but can be
#'         equivalently treated as random effects
#'
#' @return Estimate of proportion of the explained variation, variance estimates,
#'        and the corresponding confidence intervals.
#'
#' @references Restricted maximum likelihood estimator of the random effects model
#'
#' @examples \dontrun{REML(y,x)}
#'
#' @export
#'

REML=function(y,x, alpha=c(0.05),lam=1.0, niter=100,eps=1e-6){

  n = dim(x)[1]
  p = dim(x)[2]
  for (j in 1:p) {
    mu = mean(x[, j])
    sdx = sd(x[, j])
    x[, j] = (x[, j] - mu)/sdx
  }
  sdy = sd(y)
  y = (y - mean(y))/sdy

  xsvd=svd(x,nv=0) #nv=0 means not computing v matrix
  # singular value decomposition
  # $u%*%diag($d)%*%t($v)=X, t($u)%*%$u=I, t($v)%*%$v=I
  uy=t(xsvd$u)%*%y
  tau=xsvd$d^2/p-1
  dif=sum(y^2)-sum(uy^2)
  Wev=tau/(1+lam*(tau+1))^2  #vector of eigenvalues of weight matrix

  u1=t(xsvd$u)%*%rep(1,n)
  if(n>p){
    com=sum(u1^2*(Wev+1))/n-1 #negligible
    num=sum(uy^2*(Wev+1))-sum(y^2)-sum(Wev)+n-p
    den=sum(Wev*tau)+n-p
  }else{
    com=sum(u1^2*Wev)/n #negligible
    num=sum(uy^2*Wev)-sum(Wev)
    den=sum(Wev*tau)
  }
  r2=min(1,max(0,(num+com)/(den+com))) # initial value

  for(iter in 1: niter){
    if(n>=p){
      num=sum(tau*uy^2/(1+r2*tau)^2)-dif/(1-r2)^2
      num=num-((sum(uy^2/(1+r2*tau))+dif/(1-r2))/n)*(sum(tau/(1+r2*tau))-(n-p)/(1-r2))
      den=-2*(sum(tau^2*uy^2/(1+r2*tau)^3)+dif/(1-r2)^3)
      den=den+((sum(tau^2*uy^2/(1+r2*tau)^2)+dif/(1-r2)^2)/n)*(sum(tau/(1+r2*tau))-(n-p)/(1-r2))
      den=den+((sum(uy^2/(1+r2*tau))+dif/(1-r2))/n)*(sum(tau^2/(1+r2*tau)^2)+(n-p)/(1-r2))
    }else{ #No additional terms if n<=p
      num=sum(tau*uy^2/(1+r2*tau)^2)
      num=num-(sum(uy^2/(1+r2*tau))/n)*sum(tau/(1+r2*tau))
      den=-2*sum(tau^2*uy^2/(1+r2*tau)^3)
      den=den+(sum(tau^2*uy^2/(1+r2*tau)^2)/n)*sum(tau/(1+r2*tau))
      den=den+(sum(uy^2/(1+r2*tau))/n)*sum(tau^2/(1+r2*tau)^2)
    }

    factor=1
    while(r2-factor*num/den<0 | r2-factor*num/den>=1 ){
      factor=factor/2
    }
    r2=r2-factor*num/den

    if(abs(num/den)<eps){break}
  }
  #print("MLEa")
  #print(c(iter,factor,abs(num/den),r2))
  # REML estimate of the variance
  lam=r2/(1-r2)
  if(n>p){
    sigmay2=(sum(y^2)+sum(uy^2*xsvd$d^2/(p+lam*xsvd$d^2)))/n
    Iab=(-sum(y^2)+sum(uy^2)+sum(uy^2*(xsvd$d^2/p-1)/(1+lam*xsvd$d^2/p)^2))/(sigmay2^(3/2)*(1-r2)^2)
    Ib=-0.5*(sum((xsvd$d^2/p-1)^2/(1+lam*xsvd$d^2/p)^2)+n-p)/(1-r2)^2
    Ib=Ib+(sum(y^2)-sum(uy^2)+sum(uy^2*(xsvd$d^2/p-1)^2/(1+lam*xsvd$d^2/p)^3))/(sigmay2*(1-r2)^3)
  }else{
    sigmay2=sum(uy^2*p/(p+lam*xsvd$d^2))/n/(1-r2)
    Iab=sum(uy^2*(xsvd$d^2/p-1)/(1+lam*xsvd$d^2/p)^2)/(sigmay2^(3/2)*(1-r2)^2)
    Ib=-0.5*sum((xsvd$d^2/p-1)^2/(1+lam*xsvd$d^2/p)^2)/(1-r2)^2
    Ib=Ib+sum(uy^2*(xsvd$d^2/p-1)^2/(1+lam*xsvd$d^2/p)^3)/(sigmay2*(1-r2)^3)
  }
  Ia=2*n/sigmay2
print(c(Ia,Iab,Ib))
  evr2=max(0,Ia/(Ia*Ib-Iab^2))

  s2=r2*sigmay2
  evs2=max(0,evr2*sigmay2^2)

  len=length(alpha)
  cir2=r2+sqrt(evr2)*qnorm(c(alpha/2,1-alpha/2))
  cir2[c(1:len)]=cir2[c(1:len)]*(cir2[c(1:len)]>0)
  cir2[len+c(1:len)]=cir2[len+c(1:len)]*(cir2[len+c(1:len)]<1)+1.0*(cir2[len+c(1:len)]>=1)

  cis2=s2+sqrt(evs2)*qnorm(c(alpha/2,1-alpha/2))
  cis2[c(1:len)]=cis2[c(1:len)]*(cis2[c(1:len)]>0)

  #6. output result

  ind=rep(0,2*len)
  ind[2*c(1:len)-1]=c(1:len)
  ind[2*c(1:len)]=len+c(1:len)
  list(c(r2,evr2),cir2[ind],
       c(s2,evs2),cis2[ind])
}
