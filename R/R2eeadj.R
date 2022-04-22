#' Estimating equation approach to the proportion of the explained variation with two subgroups
#'
#' This approach estimates the proportion of the explained variation in a linear model
#' assuming the covariates are independent.
#'
#' @param y outcome: a vector of length n.
#' @param xa covariates to be adjusted: a matrix of nxpa dimension.
#' @param xb covariate effects to be computed: a matrix of nxpb dimension.
#' @param niter number of iterations for updating lambda. Default is 1.
#'
#' @details This method uses two estimating equations to estimate the variations of the outcome explained
#' by the two parts of covariates simultaneously.
#'
#' @return Estimates of the proportions of explained variations for xa, xb, and x=(xa,xb) respectively.
#'        No variance estimate is available currently.
#' @references To be added.
#'
#' @examples \dontrun{R2eeadj(y,xa,xb,niter=1)}
#'
#' @export
R2eeadj=function(y,xa,xb,niter=1){

  # y==outcome
  # xa==covariates to be adjusted
  # xb==covariate effects to be computed
  # niter==number of iterations for updating lambda

  n=dim(xa)[1]
  pa=dim(xa)[2]
  if(dim(xb)[1]!=n){
    # print("Stop: covariates to be adjusted does not match with covariates to be computed")
    # break
    stop("Covariates to be adjusted does not match with covariates to be computed!")
  }
  pb=dim(xb)[2]

  #1. Standardization

  for(j in 1:pa){
    mua=mean(c(xa[,j]))
    sdxa=sd(c(xa[,j]))
    xa[,j]=(xa[,j]-mua)/sdxa
  }
  for(j in 1:pb){
    mub=mean(c(xb[,j]))
    sdxb=sd(c(xb[,j]))
    xb[,j]=(xb[,j]-mub)/sdxb
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Compute the estimators
  Ma=xa%*%t(xa)/pa
  Mb=xb%*%t(xb)/pb

  r2a=0.25 #initial value
  r2b=0.25
  r2b0=0.25
  for(ii in 1:niter){
    if(r2b0<1){
      lamb0=r2b0/(1-r2b0)
    }
    if(r2a<1){
      lama=r2a/(1-r2a)
    }
    if(r2b<1){
      lamb=r2b/(1-r2b)
    }

    IMb0=chol2inv(chol(diag(rep(1,n))+lamb0*Mb))
    Wb0=IMb0%*%(Mb-diag(rep(1,n)))%*%IMb0
    denb0=sum(diag(Wb0%*%(Mb-diag(rep(1,n)))))
    numb0=t(y)%*%Wb0%*%y-sum(diag(Wb0))
    r2b0=as.numeric(numb0/denb0)
    r2b0=min(1,max(0,r2b0))

    IM=chol2inv(chol(diag(rep(1,n))+lama*Ma+lamb*Mb))
    Wa=IM%*%(Ma-diag(rep(1,n)))%*%IM
    Wb=IM%*%(Mb-diag(rep(1,n)))%*%IM

    dena=sum(diag(Wa%*%(Ma-diag(rep(1,n)))))
    denb=sum(diag(Wb%*%(Mb-diag(rep(1,n)))))
    denab=sum(diag(Wa%*%(Mb-diag(rep(1,n)))))
    denba=sum(diag(Wb%*%(Ma-diag(rep(1,n)))))

    numa=t(y)%*%Wa%*%y-sum(diag(Wa))
    numb=t(y)%*%Wb%*%y-sum(diag(Wb))

    r2a=as.numeric((denb*numa-denab*numb)/(dena*denb-denab*denba))
    r2a=min(1,max(0,r2a))
    r2b=as.numeric((dena*numb-denba*numa)/(dena*denb-denab*denba))
    r2b=min(1,max(0,r2b))
  }
  print(c(r2a,r2b,r2b0))

  #3. compute the variance estimate



  #4. output result

  list(c(r2a,r2b,r2b0))
}

