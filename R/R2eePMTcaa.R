#' Conditional permutation test adjusting for covariates for estimating equation approach
#'
#' This function performs test of no extra variation explained by the additional covariates
#' by permuting the additional covariates using the estimating equation approach.
#'
#' @param y outcome
#' @param x covariates
#' @param pa dimension of covariates to be adjusted
#' @param lam parameter adjusting the format of the weighting matrix. Default is 0.2
#' @param niter number of iterations for updating lambda. Default is 3
#' @param npm Monte Carlo sample size for permutation. Default is 1000
#'
#' @details The computation permutes the second part of the covariates and is computationally slow.
#' It differs from \code{R2eePMTca} in the way of matrix computation.
#'
#' @return Output includes the p-values (estimate and bound) of the test,
#' estimate of proportion of the extra explained variation, and simulation results.
#'
#' @references Chen, H.Y.; Li, H.; Argos, M.; Persky, V.; Turyk, M.
#'  Statistical methods for assessing explained variations of a health outcome by mixtures of exposures.
#'  Under review for Prep. Spec. Issue Int. J. Environ. Res. Public Health 2022.
#'
#' @examples \dontrun{R2eePMTcaa(y, x, pa, lam = 0.2, niter = 3, npm = 1000)}
#'
#' @export
R2eePMTcaa=function(y,x,pa,lam=0.2,niter=3,npm=1000){

  #1. Standardization

  n=dim(x)[1]
  p=dim(x)[2]
  for(j in 1:p){
    mux=mean(c(x[,j]))
    sdx=sd(c(x[,j]))
    x[,j]=(x[,j]-mux)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Compute initial estimator
  r2a=lam/(1+lam)
  for(ii in 1:niter){ #iteration to update lambda
    if(ii>1 & r2<1){lam=r2/(1-r2)}
    delta=chol2inv(chol(diag(rep(1,n))+lam*x[,1:pa]%*%t(x[,1:pa])/pa))
    sd1=sum(diag(delta))
    delta2=delta%*%delta
    sd2=sum(diag(delta2))
    den=n-2*(1+lam)*sd1+(1+lam)^2*sd2
    num=lam*(t(y)%*%delta%*%y-sd1)-lam*(1+lam)*(t(y)%*%delta2%*%y-sd2)
    r2a=num/den
    r2a=min(1,max(0,r2a))
  }

  r2=lam/(1+lam)
  for(ii in 1:niter){ #iteration to update lambda
    if(ii>1 & r2<1){lam=r2/(1-r2)}
    delta=chol2inv(chol(diag(rep(1,n))+lam*x%*%t(x)/p))
    sd1=sum(diag(delta))
    delta2=delta%*%delta
    sd2=sum(diag(delta2))
    den=n-2*(1+lam)*sd1+(1+lam)^2*sd2
    num=lam*(t(y)%*%delta%*%y-sd1)-lam*(1+lam)*(t(y)%*%delta2%*%y-sd2)
    r2=num/den
    r2=min(1,max(0,r2))
  }

  r2ba=r2-r2a
  r2ba=max(0,r2ba)

  print(c(r2,r2a,r2ba))

  result=rep(0,npm)
  for(jj in 1:npm){
    if(jj-as.integer(jj/50)*50==1){print(c(jj,jj))}

    samp=sample(n)
    xp=cbind(x[,1:pa],x[samp,(pa+1):p])

    #r2p=lam/(1+lam)
    #for(ii in 1:niter){ #iteration to update lambda
    #  if(ii>1 & r2p<1){lam=r2p/(1-r2p)}
    delta=chol2inv(chol(diag(rep(1,n))+lam*xp%*%t(xp)/p))
    sd1=sum(diag(delta))
    delta2=delta%*%delta
    sd2=sum(diag(delta2))
    den=n-2*(1+lam)*sd1+(1+lam)^2*sd2
    num=lam*(t(y)%*%delta%*%y-sd1)-lam*(1+lam)*(t(y)%*%delta2%*%y-sd2)
    r2p=num/den
    r2p=min(1,max(0,r2p))
    #  }

    result[jj]=max(0,r2p-r2a)
  }

  pvalueEST=mean(1.0*(result>r2ba))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(pvalueEST,pvalueBOUND),c(r2,r2a,r2ba),result)

}
