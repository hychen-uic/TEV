#' @import graphics
#' @import cccp
NULL
#' EigenPrism approach
#'
#' This function implements the EigenPrism procedure for estimating and generating
#' confidence intervals for variance components in high-dimensional linear model.
#'
#' @param y response vector of length n
#' @param X n by p design matrix. Columns are automatically centered and scaled to variance 1, and they cannot include
#'          one for the intercept terms.
#' @param invsqrtSig if columns of X are not independent, p by p positive definite matrix which is the square-root
#'                   of the inverse of Sig, where Sig is the *correlation* matrix of the X. Default is identity
#' @param alpha significance level for confidence interval. Default is 0.05
#' @param target target of estimation/inference: options include "beta2", "sigma2", or "heritability"
#' @param zero.ind vector of which indices of the weight vector w to constrain to zero. Default is none
#' @param diagnostics boolean variable indicating whether to generate the diagnostic plots for V_i,
#' lambda_i, and w_i. Default is \code{TRUE}.
#'
#' @details This function is a copy of Jansen's R program.
#'
#' @return estimate: unbiased estimate of the target (for heritability, only approximately unbiased);
#'         CI: 100*(1-alpha)% confidence interval for target
#'
#' @references Janson, L., Barber, R. F., Candes, E. (2017). EigenPrism: inference for high-dimensional signal-to-noise ratios.
#' *Journal of Royal Statistical Society, Ser. B.*, **79**, 1037-1065.
#' @references Lucas Janson. \url{http://lucasjanson.fas.harvard.edu/code/EigenPrism.R}.
#'
#' @examples \dontrun{EigenPrism(y, x)}
#'
#' @export
EigenPrism <- function(y, X, invsqrtSig = NULL,alpha = c(0.01, 0.05, 0.10),target = 'beta2', zero.ind = c(), diagnostics = T){
  # Author: Lucas Janson (statweb.stanford.edu/~ljanson)
  # Runs EigenPrism procedure for estimating and generating confidence
  #  intervals for variance components in high-dimensional linear model:
  #       y = X%*%beta + e,   rows of X iid~ N(0,Sig),   e iid~ N(0,sigma^2)
  #  Requires cccp package for solving second order cone optimization.
  #  Note confidence interval endpoints may lie outside parameter domain, so it may be appropriate
  #   to clip them after the fact.
  #
  # Inputs:
  #  y: response vector of length n (will automatically be centered)
  #  X: n by p design matrix; columns will automatically be centered and scaled to variance 1;
  #      should not contain intercept column, since both y and X will be centered
  #  invsqrtSig: if columns of X not independent, p by p positive definite matrix which is the square-root
  #               of the inverse of Sig, where Sig is the *correlation* matrix of the X (default is identity)
  #  alpha: significance level for confidence interval (default = 0.05)
  #  target: target of estimation/inference
  #		  'beta2' (default) is the squared 2-norm of the coefficient vector: sum(beta^2)
  #           'sigma2' is the noise variance sigma^2
  #           'heritability' is the fraction of variance of y explained by X%*%beta: t(beta)%*%Sig%*%beta/var(y)
  #  zero.ind: vector of which indices of the weight vector w to constrain to zero (default is none)
  #  diagnostics: boolean (default = T) for whether to generate diagnostic plots for the V_i, lambda_i, and w_i
  #
  # Outputs:
  #  estimate: unbiased estimate of the target (for heritability, only approximately unbiased)
  #  CI: 100*(1-alpha)% confidence interval for target

  # Get dimensionality of problem
  n = nrow(X)
  p = ncol(X)

  # Transform y and X to proper form
  y = y-mean(y)
  #X = scale(X,T,T)*n/(n-1)
  X = scale(X)[,]*n/(n-1)

  if(!is.null(invsqrtSig)) X = X%*%invsqrtSig

  # Take singular value decomposition and rescale singular values
  svd = svd(X)
  lambda = svd$d^2/p

  # Defined cone-constrained linear problem to optimize weights; [v; w] is vector of optimization variables
  q = c(1,rep(0,n)) #coefficient vector in objective function
  A = rbind(c(0,rep(1,n)),c(0,lambda)) #matrix for linear constraints
  b = c(0,1) #vector for linear constraints
  if(target=='sigma2') b = c(1,0) #switch constraints if target is sigma^2
  # Constrain some weights to be zero if desired
  if(!is.null(zero.ind)){
    A = rbind(A,cbind(rep(0,length(zero.ind)),diag(rep(1,n))[zero.ind,]))
    b = c(b,rep(0,length(zero.ind)))
  }
  # Define second-order cone constraints
  soc1 = socc(diag(c(1/4,rep(1,n))),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
  soc2 = socc(diag(c(1/4,lambda)),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
  prob = dlp(as.vector(q),A,as.vector(b),list(soc1,soc2))

  # Solve optimization problem and extract variables
  opt = cps(prob,ctrl(trace=F))
  v = getx(opt)[1]
  w = getx(opt)[-1]

  # Compute estimate and y's variance
  est = sum(w*(t(svd$u)%*%y)^2)
  yvar = sum(y^2)/n

  # Compute confidence interval
  CI=rep(0,length(alpha)*2)
  for(i in 1:length(alpha)){
    CI[(2*i-1):(2*i)]= est + yvar*sqrt(v)*qnorm(1-alpha[i]/2)*c(-1,1)
  }
  if(target=='ev'){
    est = est/yvar
    CI = CI/yvar
  }

  # Generate list with results
  result=list()
  result$estimate = est
  result$CI = CI

  # Generate diagnostic plots
  if(diagnostics){
    par(mfrow=c(1,3))

    # Check that eigenvectors are approximately Gaussian
    nV = floor(log10(n))
    srtV = svd$v[,10^(0:nV)]
    labs = c()
    for(i in 1:(nV+1)){
      srtV[,i] = sort(srtV[,i])
      ind = 10^(i-1)
      labs = c(labs,bquote(V[.(ind)]))
    }
    matplot(qnorm((1:p)/(p+1)),srtV,type="l",lwd=2,
            ylab="Quantiles of Eigenvectors",xlab="Gaussian Quantiles",
            main=expression(paste("Check Gaussianity of Eigenvectors ",V[i])))
    legend("topleft",as.expression(labs),col=1:(nV+1),lty=1:(nV+1),lwd=2)

    # Check that there are no outliers in the eigenvalues
    hist(lambda,main=expression(paste("Histogram of Normalized Eigenvalues ",lambda[i])),
         xlab=expression(lambda[i]))

    # Check that the weights are not dominated by just a few values
    srtw = sort(abs(w),T)
    plot(1:n,cumsum(srtw)/sum(srtw),type="l",lwd=2,
         main=expression(paste("Fraction of Total Weight in Largest k ",w[i])),
         xlab="k",ylab="Fraction of Total Weight")
  }

  return(result)
}
