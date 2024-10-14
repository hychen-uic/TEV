#' This function implement the maximum likelihood approach of Dicker and Erdogdu (2016).
#'
#' This the mle for the high dimensional linear model
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param alpha a vector of type I errors used to generate (1-alpha)confidence intervals.
#' @param eta2 initial value
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

RVmle=function(y,x, alpha=c(0.05),eta2=1.0,niter=100,eps=1e-6){

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
  Mev=xsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.
  Wev=(Mev-1)/(1+eta2*Mev)^2  #vector of eigenvalues of weight matrix

  u1=t(xsvd$u)%*%rep(1,n)
  if(n>p){
    com=sum(u1^2*(Wev+1))/n-1 #negligible
    num=sum(uy^2*(Wev+1))-sum(y^2)-sum(Wev)+n-p
    den=sum(Wev*(Mev-1))+n-p
  }else{
    com=sum(u1^2*Wev)/n #negligible
    num=sum(uy^2*Wev)-sum(Wev)
    den=sum(Wev*(Mev-1))
  }
  r2=min(1,max(0,(num+com)/(den+com))) # initial value

  eta2=r2/(1-r2) # initial value
  if(n>=p){ #when n>p, y^t[(I+eta2*XX^t/p)^(-1)-I]y=sum_{k=1}^p (U_k^t y)^2((1+eta2*xsvd$d^2/p)^{-1}-1)
            # This means y^t(I+eta2*XX^t/p)^(-1)y=sum_{k=1}^p (U_k^t y)^2(1+eta2*xsvd$d^2/p)^{-1}
            #                                      +y^ty-sum_{k=1}^p (U_k^t y)^2
    add=sum(y^2)-sum(uy^2)
  }else{ #No additional terms if n<=p
    add=0
  }

  for(iter in 1: niter){
    fact=1/(1+eta2*xsvd$d^2/p)

    num=sum(xsvd$d^2*fact^2*uy^2)/p-(sum(xsvd$d^2*fact)/p)*(sum(uy^2*fact)+add)/n
    den=-2*sum(xsvd$d^4*fact^3*uy^2)/p^2+(sum(xsvd$d^2*fact^2*uy^2)/p/n)*sum(xsvd$d^2*fact)/p
        +((sum(fact*uy^2)+add)/n)*sum(uy^4*fact^2)/p^2

    factor=1
    while(eta2-factor*num/den<0){
      factor=factor/2
    }
    eta2=eta2-factor*num/den

    if(abs(num/den)<eps){break}
  }
  #print("MLE")
  #print(c(iter,factor,abs(num/den),eta2))

  eta2=max(0,eta2)
  r2=eta2/(1+eta2)
  rho=p/n
  if(eta2>0){
    z=1/eta2
    A=(1-rho*z-rho+sqrt((1-rho*z-rho)^2+4*rho*z))/(2*z)
    B=-A/z-(1-(1+rho+rho*z)/sqrt((1-rho*z-rho)^2+4*rho*z))*rho/z/2
    evr2=r2^4*(1/(A^2+B)+z/rho)
  }else{
    evr2=0
  }
  evr2=evr2/n

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
