#' Estimating equation approach to the proportion of the explained variation using least square with
#'           simulated variance estimate for inference (more accurate)
#'
#' This function computes:
#'      (1). the proportion of the explained variation
#'      (2). The explained variation
#' adjusting for the correlation in covariates for n>p
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param lam parameter for altering the weighting matrix.
#' @param niter number of iterations for updating lam.
#' @param KV, the first component of the vector KV=kappa_1, the second=kappa_2, the third=kappa_3
#' @param know if KV is known, options include "yes" and "no". Default is "no".
#' @param nrep Monte Carlo sample size for computinging KV.
#' @param alpha a vector of type I errors used to generate (1-alpha)confidence intervals.
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
#' @examples \dontrun{R2eelsa(y,x,lam=0.2,niter=3,alpha=c(0.1,0.05,0.01),KV=rep(0,3),know="no",nrep=1000)}
#'
#'
#' @export
RVeelsa=function(y,x,lam=0.2,niter=1,alpha=c(0.1,0.05,0.01),KV=rep(0,3),know="no",nrep=1000){
  n=dim(x)[1]
  p=dim(x)[2]
  if(n<=p){print("n<=p, least-square approach is not applicable.")}
  else{
    if(know!='yes'){
      u = rep(1, p)/sqrt(p)
      SUZWZU = rep(0, nrep)
      SUZWWZU = rep(0, nrep)
      SUZZU = rep(0, nrep)
      STRWM = rep(0, nrep)
      for (j in 1:nrep) {
         z = matrix(rnorm(n * p), ncol = p)
         z = zscale(z)[[1]]
         zu = z %*% u
         Z=z
         SM = z %*% chol2inv(chol(t(z) %*% z + t(Z) %*% Z)) %*%
            t(z) * (n + n)/p
         ISM = chol2inv(chol(diag(rep(1, n)) + lam * SM))
         SW = ISM %*% (SM - diag(rep(1, n))) %*% ISM
         SWzu = SW %*% zu
         SUZZU[j] = sum(zu^2)
         SUZWZU[j] = t(zu) %*% SW %*% zu
         SUZWWZU[j] = t(SWzu) %*% SWzu
         STRWM[j] = sum(diag(SW %*% SM))/n
      }
      K1 = var(SUZWZU - STRWM)
      K2 = cov(SUZWZU - STRWM, SUZZU - n)
      K3 = var(SUZZU - n)
      KV=c(K1,K2,K3)
     }

    aa=RVeesd(y, x, x, lam =lam, alpha=alpha, niter = niter, know="yes", KV=KV )

    list(aa)
   }
    #list(c(r2,vestr2,vestr2n),cir2[ind],cir2n[ind],
    #   c(s2,vests2,vests2n),cis2[ind],cis2n[ind])
  #[[1]]==> R2: estimate, estimated variance, estimated variance under normal error.
  #[[2]]==>confidence interval for R2
  #[[3]]==>confidence interval for R2 under normal error.
  #[[4]]==> sigma_s^2: estimate, estimated variance, estimated variance under normal error.
  #[[5]]==>confidence interval for sigma_s^2
  #[[6]]==>confidence interval for sigma_s^2 under normal error.

}
