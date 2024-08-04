#' Estimating equation approach to the proportion of the explained variation with supplementary covariate data
#'
#' This function computes:
#'      (1). the proportion of the explained variation
#'      (2). The explained variation
#' adjusting for the correlation in covariates with possible suplementary data
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param xsup supplementary covariates: a matrix of Nxp dimension.
#'             when xsup=NULL and n>p, it performs least-square with the simulation variance estimate
#' @param lam parameter for altering the weighting matrix.
#'
#' @details The estimation approach does not assume independent covariates and can deal
#' with the case \eqn{n\le p}. But require the sample sizes of x and X combined be greater than p.
#'
#' @return Estimate of the proportion of explained variation, variance estimates under normality
#' and non-normality assumptions, and confidence intervals under normality and non-normality assumptions.
#'
#' @references Chen, H. Y., Li, H., Argos, M., Persky, V. W., and Turyk, M. (2022). Statistical Methods
#' for Assessing Explained Variation of a Health Outcome by Mixture of Exposures. International Journal
#' of Environmental Research and Public Health.
#' @references Chen H. Y, Zhang, B. and Pan, G.(2023).Estimation and inference on explained variation
#'with possible supplementary data, Manuscript.
#'
#' @examples \dontrun{RVsda(y,x,xsup,lam=0.2)}
#'
#'
#' @export
RVsda=function(y,x,xsup,lam=0.2){

  n=dim(x)[1]
  p=dim(x)[2]
  N=dim(xsup)[1]

  #1. Standardization


  for(j in 1:p){
    mu=mean(x[,j])
    sdx=sd(x[,j])
    x[,j]=(x[,j]-mu)/sdx
    xsup[,j]=(xsup[,j]-mu)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Sigular value decomposition

  xsvd=svd(x,nv=0) #nv=0 means not computing v matrix
  # singular value decomposition
  # $u%*%diag($d)%*%t($v)=X, t($u)%*%$u=I, t($v)%*%$v=I
  uy=t(xsvd$u)%*%y

  Mev=xsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.
  atemp=xsup%*%t(xsvd$v)%*%diag(1/xsvd$d)
  num=sum(uy^2*(Mev-1)/(1+lam*Mev)^2)-sum((Mev-1)/(1+lam*Mev)^2)
  den=sum(xsvd$d*(xsvd$d^2/p-1)/(1+lam*xsvd$d^2/p)^2
          *diag(chol2inv(chol(diag(rep(1,p))+t(atemp)%*%atemp))))*(n+N)/p
  if(n>p){
    num=num-sum(y*y)+sum(uy^2)
    num=num+n-p
    den=den-sum((Mev-1)/(1+lam*Mev)^2)+n-p
  }

  r2=num/den # initial value

  list(r2)

}
