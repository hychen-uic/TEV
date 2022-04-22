#' Conditional permutation test for no extra variation explained using the estimating equation approach
#'
#' This method performs permutation test by permuting the second part of the covariates and can be computationally slow.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param pa length of the covariates to be adjusted.
#' @param lam parameter adjusting the formation of the weighting matrix. Default is 0.12
#' @param niter number of iterations for updating lambda. Default is 3
#' @param npm permutation sample size for simulation computation of p-value.
#'
#' @details This method tests no extra variation explained by the second part of covariates
#' given that the first part of covariates in the model by permuting the second part of
#' covariates and estimating using the estimating equation approach. P-value is computed
#' using simulation approach. This function differs from \code{R2eePMTca} in the way the matrix operation is handled.
#'
#' @return The p-values (estimate and bound) of the test, estimate of proportion of
#' the explained variation for both parts together, the first part alone, and the second part
#' given the first part, and simulation results.
#'
#' @references Chen, H. Y., Li, H., Argos, M., Persky, V. W., and Turyk, M. (2022). Statistical Methods
#' for Assessing Explained Variation of a Health Outcome by Mixture of Exposures. International Journal
#' of Environmental Research and Public Health.
#' @references Reference 2 to be added.
#'
#' @examples \dontrun{R2eePMTcaa(y,x,pa,lam=0.1,niter=1,npm=1000)}
#'
#' @export
R2eePMTcaa=function(y,x,pa,lam=0.1,niter=1,npm=100){

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
    if(ii>1 & r2a<1){lam=r2a/(1-r2a)}
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
    result[jj]=r2p-r2a

    #r2p=min(1,max(0,r2p))
    #result[jj]=max(0,r2p-r2a)
  }
  predpvalue=1-pnorm((r2ba-mean(result))/sd(result))

  pvalueEST=mean(1.0*(result>r2ba))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(predpvalue,pvalueEST,pvalueBOUND),c(r2,r2a,r2ba),result)

}
