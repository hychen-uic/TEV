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
#' @examples \dontrun{DEMLE(y,x)}
#'
#' @export
#'

DEMLE=function(y,x, alpha=0.05,niter=100,eps=1e-5){
  n=length(y)
  n = dim(x)[1]
  p = dim(x)[2]
  for (j in 1:p) {
    mu = mean(x[, j])
    sdx = sd(x[, j])
    x[, j] = (x[, j] - mu)/sdx
  }
  sdy = sd(y)
  y = (y - mean(y))/sdy

  xsvd=svd(x,nv=0)
  temp=t(xsvd$u)%*%y
  eta2=1 # initial value
  for(iter in 1: niter){
    eta2old=eta2
    fact=1/(1+eta2*xsvd$d^2/p)
    num=sum(xsvd$d^2*fact^2*temp^2)/p-(sum(xsvd$d^2*fact)/p)*sum(temp^2*fact)/n
    den=-2*sum(xsvd$d^4*fact^3*temp^2)/p^2+(sum(xsvd$d^2*fact^2*temp^2)/p/n)*sum(xsvd$d^2*fact)/p
        +(sum(fact*temp^2)/n)*sum(temp^4*fact^2)/p^2
    eta2=eta2-num/den
    #print(c(iter,eta2,abs(num/den)))
    if(abs(num/den)<eps){break}
  }

  r2=eta2/(1+eta2)
  rho=p/n
  z=1/eta2
  A=(1-rho*z-rho+sqrt((1-rho*z-rho)^2+4*rho*z))/(2*z)
  B=-A/z-(1+(1+rho+rho*z)/sqrt((1-rho*z-rho)^2+4*rho*z))*rho/z
  ev=r2^4*(1/(A^2+B)+z/rho)

  lowci=r2-qnorm(1-alpha/2)*sqrt(ev)
  uppci=r2+qnorm(1-alpha/2)*sqrt(ev)

  list(r2,ev,c(lowci,uppci))
}
