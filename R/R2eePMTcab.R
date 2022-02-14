#' Conditional permutation test adjusting for covariates for estimating equation approach
#'
#' This function performs test of no extra variation explained by the additional covariates by
#' permuting the additional covariates using the estimating equation approach.
#'
#' @param y outcome
#' @param x covariates
#' @param pa dimension of covariates to be adjusted
#' @param lam parameter adjusting the format of the weighting matrix. Default is 0.2
#' @param niter number of iterations for updating lambda. Default is 3
#' @param npm Monte Carlo sample size for permutation. Default is 1000
#'
#' @details The computation permutes the second part of the covariates and can be computationally slow.
#' It differs from \code{R2eePMTca} in the way of matrix computation.
#'
#' @return Output includes the p-values (estimate and bound) of the test,
#' estimate of proportion of the extra explained variation, and simulation results.
#'
#' @references Chen, H.Y.; Li, H.; Argos, M.; Persky, V.; Turyk, M.
#'  Statistical methods for assessing explained variations of a health outcome by mixtures of exposures.
#'  Under review for Prep. Spec. Issue Int. J. Environ. Res. Public Health 2022.
#'
#' @examples \dontrun{R2eePMTcab(y, x, pa, lam = 0.2, niter = 3, npm = 1000)}
#'
#' @export
R2eePMTcab=function(y, x, pa, lam = 0.2, niter = 3, npm = 1000){

  #1. Standardization

  n=dim(x)[1]
  p=dim(x)[2]
  for(j in 1:p){
    mu=mean(c(x[,j]))
    sdx=sd(c(x[,j]))
    x[,j]=(x[,j]-mu)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy
  xa=x[,1:pa]

  #2. Sigular value decomposition
  r2a=0.1#lam/(1+lam)
  for(ii in 1:niter){
    if(ii>1 & r2a<1){lam=r2a/(1-r2a)}

    Ma=xa%*%t(xa)/pa
    IMa=chol2inv(chol(diag(rep(1,n))+lam*Ma))
    Wa=IMa%*%(Ma-diag(rep(1,n)))%*%IMa

    dena=sum(diag(Wa%*%(Ma-diag(rep(1,n)))))
    numa=t(y)%*%Wa%*%y-sum(diag(Wa))
    r2a=as.numeric(numa/dena)
    r2a=min(1,max(0,r2a))
    #print(c(ii,r2a))
  }

  r2=0.1#lam/(1+lam) #initial value
  for(ii in 1:niter){
    if(ii>1 & r2<1){lam=r2/(1-r2)}

    M=x%*%t(x)/p
    IM=chol2inv(chol(diag(rep(1,n))+lam*M))
    W=IM%*%(M-diag(rep(1,n)))%*%IM

    den=sum(diag(W%*%(M-diag(rep(1,n)))))
    num=t(y)%*%W%*%y-sum(diag(W))
    r2=as.numeric(num/den)
    r2=min(1,max(0,r2))
    #print(r2)
  }

  r2ba=max(0,r2-r2a)
  print(c(r2,r2a,r2ba))

  #4. Compute permutation p-value for H0: r2=0.

  result=rep(0,npm)
  for(ii in 1:npm){
    if(ii-as.integer(ii/50)*50==1){print(c(ii,ii))}

    samp=sample(n)
    xp=cbind(x[,1:pa],x[samp,(pa+1):p])

    M=xp%*%t(xp)/p
    IM=chol2inv(chol(diag(rep(1,n))+lam*M))
    W=IM%*%(M-diag(rep(1,n)))%*%IM

    #3. Compute the estimators
    den=sum(diag(W%*%(M-diag(rep(1,n)))))
    num=t(y)%*%W%*%y-sum(diag(W))
    r2p=as.numeric(num/den)
    r2p=min(1,max(0,r2p))

    result[ii]=max(0,r2p-r2a)
  }

  pvalueEST=mean(1.0*(result>r2ba))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(pvalueEST,pvalueBOUND),c(r2,r2a,r2ba),result)

}
