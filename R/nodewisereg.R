#' Approximate MLE with given precision matrix structure
#'
##' This function find the MLE of the precision matrix approximately by individual regression
##' models when the sparse structure of the precision matrix is given.
#'
#' @param dat  is a nxp data matrix
#' @param network is the fixed sparse struct of the precision matrix given
#'
#' @details This approach works for very sparse precision matrix
#'
#' @return  precision matrix estimate
#'
#' @references This program is nodewise regression for precision matrix.
#'

nodewisereg=function(dat, network)
{
  network = network + t(network)
  network = ifelse(network == 0, 0, 1)
  diag(network) = 0
  n = dim(dat)[1]
  p = dim(dat)[2]
  precision = matrix(0, nrow = p, ncol = p)
  for (j in 1:p) {
    y = dat[, j]
    index = which(network[j, ] == 1)
    if (length(index) == 0) {
      precision[j, ] = 0
      precision[j, j] = 1/var(y)
      next
    }
    x = dat[, index]
    result = lm(y ~ 0 + x)
    alpha = result$coefficients
    sigma = t(result$residuals) %*% (result$residuals)/result$df.residual
    precision[j, index] = -alpha/sigma[1]
    precision[j, j] = 1/sigma
  }
  precision = (precision + t(precision))/2

  return(precision)
}
