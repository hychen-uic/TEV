#' @import cccp
NULL
#' Estimating proportion of explained variation using the least-square approach or the EigenPrism approach
#'
#' EigenPrismFull procedure integrating the \eqn{n\le p} and \eqn{n>p} cases
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param alpha significance level: a vector.
#'
#' @details For \eqn{n\le p}, the estimate and confidence interval are obtained by EigenPrism approach.
#' For \eqn{n>p}, the estimate is obtained by least-square approach, and the confidence intervals
#' are obtained by inverting the chisquare test.
#'
#' @return Estimate of the proportion of the explained variation and confidence intervals for the proportion.
#'
#' @references Chen, H.Y. (2022). Statistical inference on explained variation in high-dimensional linear model with dense effects. arXiv:2201.08723
#' @references Janson, L., Barber, R. F., Candes, E. (2017). EigenPrism: inference for high-dimensional signal-to-noise ratios. Journal of Royal Statistical Society, Ser. B., 79, 1037-1065.
#' @references Lucas Janson. \url{http://lucasjanson.fas.harvard.edu/code/EigenPrism.R}.
#'
#' @examples \dontrun{RVep(y,x)}
#'
#' @export
#'
RVep=function(y,x,alpha=c(0.05)){

  n=dim(x)[1]
  p=dim(x)[2]
  sigy=var(y)
  ci=rep(0,2*length(alpha))

  if(n<=p){ #use EigenPrism

    aa=EigenPrism(y,x,alpha=alpha,target="ev",diagnostics=F)
    r2=aa[[1]]
    if(r2<0){r2=0}
    if(r2>1){r2=1}
    ci=aa[[2]]
    for(i in 1:length(alpha)){
      if(ci[2*i-1]<0){ci[2*i-1]=0}
      if(ci[2*i]>1){ci[2*i]=1}
    }
  }else{
    if(p==n-1){ #cannot estimate residual variance
      print("p=n-1, not computable")
      r2=0.5 # guess value
      for(i in 1:length(alpha)){
        ci[(2*i-1):(2*i)]=c(0,1)
      }
    }else{ #linear model fit
      aa=lm(y~x)
      resid=summary(aa)$sigma   # estimated residual variance
      r2=1-resid^2/sigy
      for(i in 1:length(alpha)){
        ci[(2*i-1):(2*i)]=1-(1-r2)*(n-p-1)/c(qchisq(alpha[i]/2,df=n-p-1),qchisq(1-alpha[i]/2,df=n-p-1))
        if(ci[2*i-1]<0){ci[2*i-1]=0}
        if(ci[2*i]>1){ci[2*i]=1}
      }
    }
  }
  list(r2,ci)
  #[[1]]==> Estimate
  #[[2]]==> confidence interval
}
