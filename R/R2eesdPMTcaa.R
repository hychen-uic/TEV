#' Conditional permutation test for no extra explained variation by estimating equation approach
#' with supplementary covariate data
#'
#' This function permutes the outcome and uses the estimator for the proportion of variation
#' explained for performing test of no explained variation.
#'
#' @param y outcome
#' @param x covariates
#' @param X supplementary covariate data
#' @param pa dimension of covariates to be adjusted
#' @param lam parameter adjusting the format of the weighting matrix. Default is 0.2
#' @param niter number of iteration for updating lambda. Default is 3
#' @param npm Monte Carlo sample size for permutation. Default is 1000
#'
#' @details The algorithm can be slow because the second part of the covariates need to be permuted.
#' The approach differs from the \code{R2eesdPMTca} in the way the matrix operations are handled.
#'
#' @return Output includes the p-values (estimate and bound) of the test,
#' estimate of proportion of the extra explained variation, and simulation results.
#'
#' @references Chen, H.Y.; Li, H.; Argos, M.; Persky, V.; Turyk, M.
#' Statistical methods for assessing explained variations of a health outcome by mixtures of exposures.
#' Under review for Prep. Spec. Issue Int. J. Environ. Res. Public Health 2022.
#' @references An additional reference is to be added.
#'
#' @examples \dontrun{R2eesdPMTcaa(y, x, X, pa, lam = 0.1, niter = 3, npm = 1000)}
#'
#' @export
R2eesdPMTcaa=function(y,x,X,pa,lam=0.2,niter=3,npm=1000){

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
    muxx=mean(c(XX[,j]))
    sdxx=sd(c(XX[,j]))
    XX[,j]=(XX[,j]-muxx)/sdxx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Singular value decomposition
  r2a=lam/(1+lam)
  for(ii in 1:niter){
    if(ii>1 & r2a<1){lam=r2a/(1-r2a)}

    Ma=XX[1:n,1:pa]%*%chol2inv(chol(t(XX[,1:pa])%*%XX[,1:pa]))%*%t(XX[1:n,1:pa])*(n+N)/pa
    IMa=chol2inv(chol(diag(rep(1,n))+lam*Ma))
    Wa=IMa%*%(Ma-diag(rep(1,n)))%*%IMa

    dena=sum(diag(Wa%*%(Ma-diag(rep(1,n)))))
    numa=t(y)%*%Wa%*%y-sum(diag(Wa))
    r2a=as.numeric(numa/dena)
    r2a=min(1,max(0,r2a))
  }

  r2=lam/(1+lam) #initial value
  for(ii in 1:niter){
    if(ii>1 & r2<1){lam=r2/(1-r2)}

    M=XX[1:n,]%*%chol2inv(chol(t(XX)%*%XX))%*%t(XX[1:n,])*(n+N)/p
    IM=chol2inv(chol(diag(rep(1,n))+lam*M))
    W=IM%*%(M-diag(rep(1,n)))%*%IM

    den=sum(diag(W%*%(M-diag(rep(1,n)))))
    num=t(y)%*%W%*%y-sum(diag(W))
    r2=as.numeric(num/den)
    r2=min(1,max(0,r2))
  }

  r2ba=max(0,r2-r2a)
  print(c(r2,r2a,r2ba))

  #4. Compute permutation p-value for H0: r2=0.

  result=rep(0,npm)
  for(ii in 1:npm){
    if(ii-as.integer(ii/50)*50==1){print(c(ii,ii))}

    samp=sample(n+N)
    XXp=cbind(XX[,1:pa],XX[samp,(pa+1):p])

    M=XXp[1:n,]%*%chol2inv(chol(t(XXp)%*%XXp))%*%t(XXp[1:n,])*(n+N)/p
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
