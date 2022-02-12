#' Permutation test for the two-group estimating equation approach (Need further work)
#'
#' One paragraph describing this function...
#'
#' @param y outcome
#' @param xa covariates to be adjusted
#' @param xb covariate effects to be computed
#' @param niter number of iterations for updating lambda. Default is 5
#' @param npm permutation sample size. Default is 1000
#' @param proja matrix projecting to orthogonal complement of xa
#'
#' @details Details of this function...
#'
#' @return Output of this function...
#'
#' @references To be added.
#'
#' @examples \dontrun{R2eeadjPMT(y,xa,xb,niter=5,npm=1000,proja)}
#'
#' @export
R2eeadjPMT=function(y,xa,xb,niter=5,npm=1000,proja){

  # y==outcome
  # xa==covariates to be adjusted
  # xb==covariate effects to be computed
  # niter==number of iterations for updating lambda
  # npm==permutation sample size
  # proja== matrix projecting to orthogonal complement of xa

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

  yr=proja%*%y    # projection
  xrb=proja%*%xb  # projection
  for(j in 1:pb){
    mub=mean(c(xrb[,j]))
    sdb=sd(c(xrb[,j]))
    xrb[,j]=(xrb[,j]-mub)/sdb
  }
  sdy=sd(yr)
  yr=(yr-mean(yr))/sdy

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
  r2=c(r2a,r2b)

  #4. Compute permutation p-value for H0: r2=0.
  result=array(0,c(npm,2))
  for(ii in 1:npm){
    samp=sample(n)
    yy=y[samp]

    numa=t(yy)%*%Wa%*%yy-sum(diag(Wa))
    numb=t(yy)%*%Wb%*%yy-sum(diag(Wb))

    result[ii,1]=as.numeric((denb*numa-denab*numb)/(dena*denb-denab*denba))
    result[ii,1]=min(1,max(result[ii,1],0))

    result[ii,2]=as.numeric((dena*numb-denba*numa)/(dena*denb-denab*denba))
    result[ii,2]=min(1,max(result[ii,2],0))
  }

  crt=2*sqrt(0.05*0.95/npm) # accuracy if truth p-value=0.05
  pvalueESTa=mean(1.0*(result[,1]>r2a))
  acc=2*sqrt(pvalueESTa*(1-pvalueESTa)/npm)
  pvalueBOUNDa=max(pvalueESTa+acc,pvalueESTa+crt)

  pvalueESTb=mean(1.0*(result[,1]>r2b))
  acc=2*sqrt(pvalueESTb*(1-pvalueESTb)/npm)
  pvalueBOUNDb=max(pvalueESTb+acc,pvalueESTb+crt)

  list(matrix(c(pvalueESTa,pvalueESTb,pvalueBOUNDa,pvalueBOUNDb),ncol=2),r2,result)

}
