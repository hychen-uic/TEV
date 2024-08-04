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

RVmle=function(y,x, alpha=c(0.05),niter=100,eps=1e-5){

  n = dim(x)[1]
  p = dim(x)[2]
  for (j in 1:p) {
    mu = mean(x[, j])
    sdx = sd(x[, j])
    x[, j] = (x[, j] - mu)/sdx
  }
  sdy = sd(y)
  y = (y - mean(y))/sdy

  lam=0.1
  Xsvd=svd(x,nv=0) #nv=0 means not computing v matrix
  # singular value decomposition
  # $u%*%diag($d)%*%t($v)=X, t($u)%*%$u=I, t($v)%*%$v=I
  uy=t(Xsvd$u)%*%y
  Mev=Xsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.
  if(n>p){
    num=sum(uy^2*(Mev-1)/(1+lam*Mev)^2)-sum(y*y)+sum(uy^2)
    num=num-sum((Mev-1)/(1+lam*Mev)^2)+n-p
    den=sum((Mev-1)^2/(1+lam*Mev)^2)+n-p
  }else{
    num=sum(uy^2*(Mev-1)/(1+lam*Mev)^2)
    num=num-sum((Mev-1)/(1+lam*Mev)^2)
    den=sum((Mev-1)^2/(1+lam*Mev)^2)
  }
  r2=num/den # initial value


  xsvd=svd(x,nv=0)
  temp=t(xsvd$u)%*%y
  eta2=1.0 # initial value
  if(n>=p){ #when n>p, y^t[(I+eta2 XX^t/p)^(-1)-I]y=sum_{k=1}^p (U_k^t y)^2((1+eta2*xsvd$d^2/p)^{-1}-1)
            # This means y^t(I+eta2 XX^t/p)^(-1)y=sum_{k=1}^p (U_k^t y)^2(1+eta2*xsvd$d^2/p)^{-1}
            #                                      +y^ty-sum_{k=1}^p (U_k^t y)^2
    add=sum(y^2)-sum(temp^2)
  }else{ #No additional terms if n<=p
    add=0
  }
  for(iter in 1: niter){
    eta2old=eta2
    fact=1/(1+eta2*xsvd$d^2/p)

    num=sum(xsvd$d^2*fact^2*temp^2)/p-(sum(xsvd$d^2*fact)/p)*(sum(temp^2*fact)+add)/n
    den=-2*sum(xsvd$d^4*fact^3*temp^2)/p^2+(sum(xsvd$d^2*fact^2*temp^2)/p/n)*sum(xsvd$d^2*fact)/p
        +((sum(fact*temp^2)+add)/n)*sum(temp^4*fact^2)/p^2
    eta2=eta2-num/den
    #print(c(iter,eta2,abs(num/den)))
    if(abs(num/den)<eps){break}
  }

  eta2=max(0,eta2)
  r2=eta2/(1+eta2)
  rho=p/n
  if(eta2>0){
    z=1/eta2
    A=(1-rho*z-rho+sqrt((1-rho*z-rho)^2+4*rho*z))/(2*z)
    B=-A/z-(1+(1+rho+rho*z)/sqrt((1-rho*z-rho)^2+4*rho*z))*rho/z
    evr2=r2^4*(1/(A^2+B)+z/rho)
  }else{
    evr2=0
  }

  vy=as.numeric(var(y))
  s2=vy*r2
  evs2=vy^2*evr2+r2^2*as.numeric(var(y^2))/n # This variance estimator is not consistent

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
