#' This function implement the maximum likelihood approach of Dicker and Erdogdu (2016).
#'
#' This the mle for the high dimensional linear model
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param alpha a vector of type I errors used to generate (1-alpha)confidence intervals.
#' @param niter the number of iterations for finding the signal noise ratio
#' @param eps the convergence criterion for the iteration
#'
#' @details  This method assume the independent covariates with fixed effects but can be
#'         equivalently treated as random effects
#'
#' @return Estimate of proportion of the explained variation, variance estimates,
#'        and the corresponding confidence intervals.
#'
#' @references Dicker, L. H. and Erdogdu, M. A. (2016). Maximum likelihood for
#'           variance estimation in high-dimensional linear models.
#'  Proceedings of the 19th International Conference on Articial Intelligence and Statistics
#'
#' @examples \dontrun{RVMLE(y,x)}
#'
#' @export
#'

RVmlea=function(y,x, alpha=c(0.05),niter=100,eps=1e-6){

  n = dim(x)[1]
  p = dim(x)[2]
  for (j in 1:p) {
    mu = mean(x[, j])
    sdx = sd(x[, j])
    x[, j] = (x[, j] - mu)/sdx
  }
  sdy = sd(y)
  y = (y - mean(y))/sdy

  r2=0.1
  xsvd=svd(x,nv=0) #nv=0 means not computing v matrix
  # singular value decomposition
  # $u%*%diag($d)%*%t($v)=X, t($u)%*%$u=I, t($v)%*%$v=I
  uy=t(xsvd$u)%*%y
  tau=xsvd$d^2/p-1
  dif=sum(y^2)-sum(uy^2)
  Wev=tau/(1+(r2/(1-r2))*(tau+1))^2  #vector of eigenvalues of weight matrix

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
    print("MLEa")
    print(c(iter,r2,factor,abs(num/den)))
    if(abs(num/den)<eps){break}
  }

  rho=p/n
  if(r2>0){
    z=(1-r2)/r2
    A=(1-rho*z-rho+sqrt((1-rho*z-rho)^2+4*rho*z))/(2*z)
    B=-A/z-(1+(1+rho+rho*z)/sqrt((1-rho*z-rho)^2+4*rho*z))*rho/z
    evr2=r2^4*(1/(A^2+B)+z/rho)
  }else{
    evr2=0
  }

  s2=sdy^2*r2
  evs2=sdy^4*evr2+r2^2*sdy^2/n # This variance estimator is not consistent

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
