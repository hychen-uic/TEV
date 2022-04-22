#' Permutation test for no explained variation using the estimating equation approach with supplementary covariate data
#'
#' This method performs permutation test by permuting the outcome values.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param X supplementary covariates: a matrix of Nxp dimension.
#' @param lam parameter adjusting the formation of the weighting matrix.
#' @param niter number of iterations for updating lambda.
#' @param npm permutation sample size for simulation computation of p-value.
#'
#' @details This method tests no explained variation by permuting the outcome and estimating
#' using the estimating equation approach with supplementary data. P-value is computed using simulation approach.
#'
#' @return The p-values (normal approximation, empirical estimate, and bound) of the test, estimate of proportion of
#' the explained variation, and simulation results.
#'
#' @references Chen, H. Y., Li, H., Argos, M., Persky, V. W., and Turyk, M. (2022). Statistical Methods
#' for Assessing Explained Variation of a Health Outcome by Mixture of Exposures. International Journal
#' of Environmental Research and Public Health.
#' @references Reference 2 to be added.
#'
#' @examples \dontrun{R2eesdPMT(y,x,X,lam=0.1,niter=1,npm=1000)}
#'
#' @export
R2eesdPMT=function(y,x,X,lam=0.1,niter=1,npm=1000){

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
    x[,j]=(x[,j]-muxx)/sdxx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Sigular value decomposition
  M=x%*%chol2inv(chol(t(XX)%*%XX))%*%t(x)*(n+N)/p

  r2=lam/(1+lam) #initial value
  for(ii in 1:niter){
    if(ii>1 & r2<1){lam=r2/(1-r2)}

    IM=chol2inv(chol(diag(rep(1,n))+lam*M))
    W=IM%*%(M-diag(rep(1,n)))%*%IM
    #3. Compute the estimators
    den=sum(diag(W%*%(M-diag(rep(1,n)))))
    num=t(y)%*%W%*%y-sum(diag(W))

    r2=as.numeric(num/den)
    r2=min(1,max(0,r2))
  }

  #4. Compute permutation p-value for H0: r2=0.

  result=rep(0,npm)
  for(ii in 1:npm){
    samp=sample(n)
    yy=y[samp]
    num=t(yy)%*%W%*%yy-sum(diag(W))

    result[ii]=as.numeric(num/den)
    result[ii]=min(1,max(result[ii],0))
  }
  predpvalue=1-pnorm((r2-mean(result))/sd(result))

  pvalueEST=mean(1.0*(result>r2))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(predpvalue,pvalueEST,pvalueBOUND),r2,result)
}
