#' This function implement the generalized random effects model approach.
#'
#' This function uses least-square estimates in computing the proportion of the explained variation.
#' It is simply just the least-square approach to estimate heritability.The general GRE assumes
#' covariates are blockwise independent. This function can be applied to one block. Variance estimate
#' is good for normal outcome only.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param alpha a vector of type I errors used to generate (1-alpha)confidence intervals.
#'
#' @details This method works only for the case n>p. It uses the least-square approach
#' for the estimation. Covariates are allowed to be correlated.
#'
#' @return Estimate of proportion of the explained variation, variance estimates,
#'        and the corresponding confidence intervals.
#'
#' @references Hou et al (2019).  Accurate estimation of SNP-heritability from
#' biobank-scale data irrespective of genetic architecture. NAture GeNeticS, 51, 1244–1251.
#'
#' @examples \dontrun{GRE(y,x)}
#'
#' @export
#'

GRE=function(y,x, alpha=0.05){
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

  beta=t(x)%*%y/n
  xsvd=svd(x)
  q=sum(1.0*(xsvd$d!=0))
  temp=t(xsvd$v)%*%beta
  #xxginv=xsvd$v%*%diag(1/xsvd$d^2)%*%t(xsvd$v)
  #r2=(n*t(beta)%*%xxginv%*%beta-q)/(n-q)
  r2=(n^2*sum(temp^2/xsvd$d^2)-q)/(n-q)
  ev=(1-r2)^2*2*q/(n-q)^2+4*r2*(1-r2)*n/(n-q)^2
  lowci=r2-qnorm(1-alpha/2)*sqrt(ev)
  uppci=r2+qnorm(1-alpha/2)*sqrt(ev)

  list(r2,ev,c(lowci,uppci))
}
