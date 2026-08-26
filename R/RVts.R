#' @import SILGGM
NULL
#' Estimating equation approach to the proportion of the explained variation,
#'    and the explained variation with estimated structured covariance matrix
#'
#' This function estimates:
#'     (1). the proportion of the explained variation
#'     (2). the explained variation
#' by covariates in a linear model assuming the covariates are correlated and the
#' covariance are estimated using the latent variable followed by
#' sparse precision matrix assumption.
#'
#' @param y outcome: a vector of length n.
#' @param x covariates: a matrix of nxp dimension.
#' @param xsup supplementary covariates (Nxp) if any.
#' @param lam parameter for altering the weighting matrix.
#' @param niter number of iterations for updating the lam parameter
#' @param alpha a vector of type I error for create the confidence intervals.
#' @param maxpc a number denotes the maximum number principal components allowed.
#'
#' @details  Both point estimate and confidence intervals
#' are computed. Two set of confidence intervals under normal or non-normal error are computed.
#'
#' @return Estimate of the proportion of explained variation, variance estimates, and confidence intervals,
#'  under normality and non-normality assumptions.
#'
#' @return Estimate of the explained variation, variance estimates, and confidence intervals
#'  under normality and non-normality assumptions.
#'
#' @references Chen, H.Y. (2026). Estimation of Proportion of Explained Variation by High-dimensional
#'  covariates with Structured covariance Matrix.
#'
#' @examples \dontrun{RVts(y,x)}
#'
#'@export
RVts=function(y,x,xsup=NULL,lam=1,niter=0,alpha=c(0.05), maxpc=10){
  n=dim(x)[1]
  p=dim(x)[2]
  if(is.null(xsup)){
    N=0
    XX=x
  }else{
    N=dim(xsup)[1]
    xx=rbind(x,xsup) #combine existing and supplement data on covariates

    if(dim(xsup)[2]!=p){
      # print("Stop: supplement data dimension does not match")
      # break
      stop("Supplement data dimension does not match!")
    }
  }
  #1. Standardization


  for(j in 1:p){
    mu=mean(xx[,j])
    sdx=sd(xx[,j])
    xx[,j]=(xx[,j]-mu)/sdx
    x[,j]=(x[,j]-mu)/sdx
  }

  sdy=sd(y)
  y=(y-mean(y))/sdy

  if(N>0){xx=rbind(x,xsup)}else{xx=x}

  # 1. identify latent variables (Principal components)

  svdx=svd(xx)
  npc=0
  for(i in 1:maxpc){
    if(svdx$d[i]>2*svdx$d[i+1]-svdx$d[i+10]){npc=i}
  }

  if(npc>=1){
    pc=xx%*%svdx$v[,1:npc]
    fit=lm(y~pc)
    yres=as.numeric(fit$residual)
    exvar=var(y-yres)/var(y) #proportion of variation explained by all pcs.
    coef=as.vector(fit$coefficients)[-1]
    residvar=(summary(fit)$sigma)^2/n
    vexvar=2*t(coef)%*%(t(pc)%*%pc/n)%*%coef*residvar/(var(y))^2/n # need to be updated
    # since diag(a number)=identity matrix with that dimension, it has to be cautious
    xx=xx-as.matrix(svdx$u[,1:npc])%*%diag(svdx$d)[1:npc,1:npc]%*%t(as.matrix(svdx$v[,1:npc]))
    # remove the PC components from x
  }else{
    yres=y
    exvar=0
    vexvar=0
  }


  # 2. identify sparse precision matrix for the residual data

  fit=SILGGM(xx[1:n,],alpha=alpha,global=TRUE)
  network=fit$global_decision[[1]]
  Omega=nodewisereg(dat=xx[1:n,], network=network) #refit for sparse covariance matrix

  sqrtOmega=sqrtpdm(Omega)[[1]]
  #print(sum(abs(eigen(x%*%Omega%*%t(x)/p)$values-eigen(x%*%solve(sqrtsig%*%sqrtsig)%*%t(x)/p)$values)))

  z=xx[1:n,]%*%sqrtOmega

  aa=TEV::RVee(yres,z,alpha=alpha,lam=lam,niter=1)
  r2=aa[[1]][1]*(1-exvar)+exvar   # estimator of R2, variance estimate under normal, variance estimate
  vestr2=aa[[1]][2]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2
  vestr2n=aa[[1]][3]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2


  s2=r2*var(y)
  vests2=vestr2*var(y)^2
  vests2n=vestr2n*var(y)^2

  # proportion of explained variation
  len=length(alpha)

  cir2=r2+sqrt(vestr2)*qnorm(c(alpha/2,1-alpha/2))
  # cut off below 0 and exceeding 1
  cir2[c(1:len)]=cir2[c(1:len)]*(cir2[c(1:len)]>0)
  cir2[len+c(1:len)]=cir2[len+c(1:len)]*(cir2[len+c(1:len)]<1)+1.0*(cir2[len+c(1:len)]>=1)

  cir2n=r2+sqrt(vestr2n)*qnorm(c(alpha/2,1-alpha/2)) # under normal
  # cut off below 0 and exceeding 1
  cir2n[c(1:len)]=cir2n[c(1:len)]*(cir2n[c(1:len)]>0)
  cir2n[len+c(1:len)]=cir2n[len+c(1:len)]*(cir2n[len+c(1:len)]<1)+1.0*(cir2n[len+c(1:len)]>=1)

  # explained variation

  cis2=s2+sqrt(vests2)*qnorm(c(alpha/2,1-alpha/2))
  # cut off below 0
  cis2[c(1:len)]=cis2[c(1:len)]*(cis2[c(1:len)]>0)

  cis2n=s2+sqrt(vests2n)*qnorm(c(alpha/2,1-alpha/2)) # under normal
  # cut off below 0
  cis2n[c(1:len)]=cis2n[c(1:len)]*(cis2n[c(1:len)]>0)


  #5. output result

  ind=rep(0,2*len)
  ind[2*c(1:len)-1]=c(1:len)
  ind[2*c(1:len)]=len+c(1:len)
  list(c(r2,vestr2,vestr2n),cir2[ind],cir2n[ind],
       c(s2,vests2,vests2n),cis2[ind],cis2n[ind])
}
