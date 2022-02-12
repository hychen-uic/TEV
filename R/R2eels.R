#' Estimate the proportion of explained variation by the least squares approach
#'
#' This function works only for n > p. But covariates are not required to be independent.
#'
#' @param y outcome
#' @param x covariates
#'
#' @details This method uses the least square approach.
#'
#' @return Output includes estimates, variance, and confidence intervals.
#'
#' @references Chen, H.Y. (2022). Statistical inference on explained variation in high-dimensional
#' linear model with dense effects. arXiv:2201.08723
#'
#' @examples \dontrun{R2eels(y,x)}
#'
#' @export
R2eels=function(y,x){

  # y==outcome
  # x==covariates

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

  #2. Compute the estimator
  xsvd=svd(x)
  r2=1-(n-1-sum((t(y)%*%xsvd$u)^2))/(n-p)
  r2=min(1,max(0,r2))

  #3. estimate variance
  # Assuming normal random error
  #vest=2*n*(1-r2)^2/(n-p)
  vest=2*(1-r2)^2*(2*r2^2-1)+2*(1-r2)^2*n/(n-p)


  # Not assuming normal random error

  zext=cbind(x,matrix(0,nrow=n,ncol=n-p))  # create an extended matrix
  svdzext=svd(zext)
  #Find the remaining orthogonal matrix in svd(x)$u
  U=svdzext$u[1:n,(p+1):n]
  #residual error terms
  eps=t(U)%*%y
  # variance estimate without normality for random error
  veps2=(sum(eps^4)-3*(n-p)*(1-r2)^2)/sum(U^4)+2*(1-r2)^2
  #vest1=2*p*(1-r2)^2/(n-p)+max(veps2,0)
  Delta=diag(xsvd$u%*%t(xsvd$u))
  vest1=vest+(max(veps2,0)-2)*sum(((1-r2)^2-(1-r2)*(1-Delta)*n/(n-p))^2)/n

  vest=vest/n
  vest1=vest1/n

  #variance for logit transformed r^2

  ci=r2+sqrt(vest)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95)) #normal
  ci[2*c(1:3)-1]=ci[2*c(1:3)-1]*(ci[2*c(1:3)-1]>0)
  ci[2*c(1:3)]=ci[2*c(1:3)]*(ci[2*c(1:3)]<1)+1.0*(ci[2*c(1:3)]>=1)

  ci1=r2+sqrt(vest1)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95)) #non-normal
  ci1[2*c(1:3)-1]=ci1[2*c(1:3)-1]*(ci1[2*c(1:3)-1]>0)
  ci1[2*c(1:3)]=ci1[2*c(1:3)]*(ci1[2*c(1:3)]<1)+1.0*(ci1[2*c(1:3)]>=1)

  #5. output result

  list(c(r2,vest,vest1),ci,ci1)
  #[[1]]==> estimate, estimated variance (under normal), estimated variance.
  #[[2]]==>confidence interval under normal
  #[[3]]==>confidence interval without normal
}
