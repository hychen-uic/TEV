#' Generate data for simulation studies.
#'
#' This function generate data follows linear model based on parameter input.
#'
#' @param n sample size
#' @param p number of covariates
#' @param beta regression coefficients
#' @param xsig standard deviation of x variables
#' @param errsig residual standard deviation
#' @param powx power for x variable transformation
#' @param powy power of y variable transformation
#' @param sqrtsig square root decomposition of correlation matrix
#'
#' @return  x covariates.
#'          y outcome variable
#' @export
#'
makedat=function(n,p,beta,xsig,errsig,powx,powy,sqrtsig){
  #generate data from linear model
  # powx, powy parameters making non-normal x and random error
  # sd, square-root of correlation matrix, making correlated x.

  x=array(0,c(n,p))
  for(j in 1:p){
    x[,j]=rnorm(n,mean=0,sd=xsig)
    }
   x=sign(x)*abs(x)^powx # simulate correlated covariates
   if(powx==2){
     x=abs(x)
     }
   for(j in 1:p){ #standardization
     mu=mean(x[,j])
     sdx=sd(x[,j])
     x[,j]=(x[,j]-mu)/sdx
     }

   x=x%*%sqrtsig  # create correlated covariates
   x=x*(x<200)+200*(x>=200)

   err=rnorm(n,mean=0,sd=errsig)
   err=sign(err)*abs(err)^powy
   if(powy==2){
     err=abs(err)
     }
   err=err*(err<200)+200*(err>=200)

   y=x%*%beta+err
   y=as.numeric(y)

  list(x,y)
}

#' Generate data following high-dimensional linear model.
#'
#' This function generates data following a linear model where covariates follow
#'  an auto regressive model.
#'
#' @param n sample size
#' @param p number of covariates
#' @param beta regression coefficients
#' @param xsig standard deviation of x variables
#' @param errsig residual standard deviation
#' @param powx power for x variable transformation
#' @param powy power of y variable transformation
#' @param k the length of dependence in the auto-regression
#' @param a autoregressive coefficients of length k
#'
#' @return  x covariates.
#'          y outcome variable
#' @export
#'
makedata=function(n,p,beta,xsig,errsig,powx,powy,k,a=rep(0.5,k)){
  #generate data from linear model
  # powx, powy parameters making non-normal x and random error
  # sd, square-root of correlation matrix, making correlated x.

  x=matrix(rnorm(n*p),ncol=p)
  x=sign(x)*abs(x)^powx # simulate correlated covariates
  if(powx==2){
    x=abs(x)
  }
  x=x*(x<200)+200*(x>=200)

  for(j in (k+1):p){ # use auto regression to create correlation
    for(i in 1:k){
      x[,j]=x[,j]+a[k]*x[,j-k]
    }
  }
  for(j in 1:p){ #standardization
    mu=mean(x[,j])
    sdx=sd(x[,j])
    x[,j]=(x[,j]-mu)/sdx
  }

  err=rnorm(n,mean=0,sd=errsig)
  err=sign(err)*abs(err)^powy
  if(powy==2){
    err=abs(err)
  }
  err=err*(err<200)+200*(err>=200)

  y=x%*%beta+err
  y=as.numeric(y)

  list(x,y)
}


#' generate a square root of a correlation matrix
#'
#' This function generate data follows linear model based on parameter input.
#'
#' @param amply parameter controlling the correlation
#' @param tilt parameter controlling the correlation
#' @param shrink parameter(>=0) controlling the correlation
#' @param p dimensional of the matrix
#'
#' @return  corr correlaation matrix
#'
#' @export
#'
makecora=function(amply,tilt,shrink,p){

  x=matrix(rnorm(p*p,mean=amply,sd=1),ncol=p)
  y=matrix(runif(p*p),ncol=p)-0.5+tilt
  xy=x%*%y
  covar=t(xy)%*%xy
  covar=covar+diag(shrink*diag(covar))
  covar=abs(covar)
  corr=diag(1/sqrt(diag(covar)))%*%covar%*%diag(1/sqrt(diag(covar)))

  list(corr)
}

#' generate a square root of a correlation matrix by specifying correlation matrix
#'
#' This function generate data follows linear model based on parameter input.
#'
#' @param rho parameter controlling the correlation
#' @param fix the correlation is a constant if fix=TRUE, Otherwise, it is rho^|i-j|
#' @param p dimensional of the matrix
#'
#' @return  corr correlation matrix
#'
#' @export
#'
makecorb=function(rho,p,fix=TRUE){

  corr=matrix(0,ncol=p,nrow=p)

  corr[1,1]=1
  for(i in 2:p){
    corr[i,i]=1
    for(j in 1:(i-1)){
      if(fix==TRUE){
        corr[i,j]=rho
      }else{
        corr[i,j]=rho^(i-j)
      }
      corr[j,i]=corr[i,j]
    }
  }
  list(corr)
}

#' Find the square root of a positive definite matrix
#'
#' This function find a square root of the positive definite matrix by single value decomposition.
#'
#' @param corr correlation or covariance matrix
#'
#' @return  sqrtsig square root of a covariance matrix
#'
#' @export
#'
sqrtpdm=function(corr){

  svdcorr=svd(corr) # singular value decomposition
  sqrtsig=svdcorr$u%*%diag(sqrt(svdcorr$d))%*%t(svdcorr$v)# square-root a matrix

  list(sqrtsig)
}

#' Calculation of the true r2 and explained variance v2 by simulation
#'
#' @param nrep simulation sample size
#' @param p number of covariates
#' @param beta regression coefficients
#' @param xsig standard deviation of x variables
#' @param errsig residual standard deviation
#' @param powx power for x variable transformation
#' @param powy power of y variable transformation
#' @param sd square root decomposition of correlation matrix
#'
#' @return  x covariates.
#'          y outcome variable
#' @export
#'
trueR2=function(nrep,p,beta,xsig,errsig,powx,powy,sd){

  xy=makedat(nrep,p,beta,xsig,errsig,powx,powy,sd)
  s2=as.numeric(var(xy[[1]]%*%beta))
  r2=as.numeric(s2/var(xy[[2]]))

  list(r2,s2)
}



