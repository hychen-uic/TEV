#' Permutation test for no explained variation using the estimating equation approach
#'
#' This method performs permutation test by permuting the outcome values.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param lam parameter adjusting the formation of the weighting matrix. Default is 0.12.
#' @param niter number of iterations for updating lambda. Default is 3.
#' @param npm permutation sample size for simulation computation of p-value.
#'
#' @details This method tests no explained variation by permuting the outcome and estimating
#' using the estimating equation approach. P-value is computed using simulation approach.
#'
#' @return The p-values (estimate and bound) of the test, estimate of proportion of
#' the explained variation, and simulation results.
#'
#' @references Chen, H. Y., Li, H., Argos, M., Persky, V. W., and Turyk, M. (2022). Statistical Methods
#' for Assessing Explained Variation of a Health Outcome by Mixture of Exposures. International Journal
#' of Environmental Research and Public Health.
#' @references Reference 2 to be added.
#'
#' @examples \dontrun{R2eePMT(y,x,lam=0.1,niter=1,npm=1000)}
#'
#' @export
R2eePMT=function(y, x, lam =0.1, niter = 1, npm = 1000){

  n=dim(x)[1]
  p=dim(x)[2]

  #1. Standardization

  for(j in 1:p){
    mu=mean(x[,j])
    sdx=sd(x[,j])
    x[,j]=(x[,j]-mu)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Sigular value decomposition
  Xsvd=svd(x,nv=0) #nv=0 means not computing v matrix
  # singular value decomposition
  # $u%*%diag($d)%*%t($v)=X, t($u)%*%$u=I, t($v)%*%$v=I
  Mev=Xsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.

  r2=lam/(1+lam) # initial value
  for(ii in 1:niter){ #iteration to update lambda
    if(ii>1 & r2<1){lam=r2/(1-r2)}
    Dev=(Mev-1)/(1+lam*Mev)^2  #vector of eigenvalues of weight matrix
    #3. Compute the estimators
    uy=t(Xsvd$u)%*%y
    u1=t(Xsvd$u)%*%rep(1,n)
    if(n>=p){
      com=sum(u1^2*(Dev+1))/n-1
      num=sum(uy^2*(Dev+1))-sum(y^2)-sum(Dev)+n-p
      den=sum(Dev*(Mev-1))+n-p
    }else{
      com=sum(u1^2*Dev)/n
      num=sum(uy^2*Dev)-sum(Dev)
      den=sum(Dev*(Mev-1))
    }

    r2=(num+com)/(den+com)
    r2=min(1,max(0,r2))
  }
  #4. Compute permutation p-value for H0: r2=0.
  result=rep(0,npm)
  for(ii in 1:npm){
    samp=sample(n)
    yy=y[samp]
    uy=t(Xsvd$u)%*%yy
    if(n>=p){
      num=sum(uy^2*(Dev+1))-sum(yy^2)-sum(Dev)+n-p
    }else{
      num=sum(uy^2*Dev)-sum(Dev)
    }

    result[ii]=(num+com)/(den+com)
    result[ii]=min(1,max(result[ii],0))
  }
  predpvalue=1-pnorm((r2-mean(result))/sd(result))

  pvalueEST=mean(1.0*(result>r2))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(predpvalue,pvalueEST,pvalueBOUND),r2,result)
}
