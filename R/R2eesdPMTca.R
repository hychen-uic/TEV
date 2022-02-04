#' Conditional permutation test fo the estimating equation approach adjusting for covariates with supplementary data
#'
#' One paragraph describing this function
#'
#' @param y outcome
#' @param x covariates
#' @param X supplementary covariate data
#' @param pa dimension of covariates to be adjusted
#' @param lam parameter adjusting the format of the weighting matrix. Default is 0.2
#' @param niter number of iteration for updating lambda
#' @param npm Monte Carlo sample size for permutation
#'
#' @details Details of this function...
#'
#' @return Output of this function...
#'
#' @references Reference 1 here...
#' @references Reference 2 here...
#'
#' @examples \dontrun{R2eesdPMTca(y,x,X,pa,lam=0.2,niter=3,npm=1000)}
#'
#' @export
R2eesdPMTca=function(y,x,X,pa,lam=0.2,niter=3,npm=1000){

  n=dim(x)[1]
  p=dim(x)[2]
  N=dim(X)[1]
  if(dim(X)[2]!=p){
    # print("Stop: supplement data dimension does not match")
    # break
    stop("Supplement data dimension does not match!")
  }
  #1. Standardization

  XX=rbind(x,X) #combine existing and supplement data on covariates
  for(j in 1:p){
    muxx=mean(c(XX[,j]))
    sdxx=sd(c(XX[,j]))
    XX[,j]=(XX[,j]-muxx)/sdxx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Sigular value decomposition
  r2a=lam/(1+lam)
  for(ii in 1:niter){
    if(ii>1 & r2a<1){lam=r2a/(1-r2a)}

    delta=(1+lam*(n+N)/pa)*t(x[,1:pa])%*%x[,1:pa]+t(X[,1:pa])%*%X[,1:pa]
    invd=chol2inv(chol(delta))
    xy=as.vector(t(x[,1:pa])%*%y)
    invdxy=as.vector(invd%*%xy)
    xx=t(x[,1:pa])%*%x[,1:pa]
    invdxx=invd%*%xx
    sd1=sum(diag(invdxx))
    sd2=sum(diag(invdxx%*%invdxx))
    den=n-2*sd1*(1+lam)*(n+N)/pa+sd2*(n+N)^2*(1+lam)^2/pa^2
    num=-sum(y^2)+(n+N)*(1+2*lam)*sum(xy*invdxy)/pa-lam*(1+lam)*(n+N)^2*as.numeric(t(invdxy)%*%xx%*%invdxy)/pa^2
    num=num+n-sd1*(1+2*lam)*(n+N)/pa+sd2*(n+N)^2*(1+lam)*lam/pa^2
    r2a=as.numeric(num/den)
    r2a=min(1,max(0,r2a))
  }


  r2=lam/(1+lam) #initial value
  for(ii in 1:niter){
    if(ii>1 & r2<1){lam=r2/(1-r2)}

    delta=(1+lam*(n+N)/p)*t(x)%*%x+t(X)%*%X
    invd=chol2inv(chol(delta))
    xy=as.vector(t(x)%*%y)
    invdxy=as.vector(invd%*%xy)
    xx=t(x)%*%x
    invdxx=invd%*%xx
    sd1=sum(diag(invdxx))
    sd2=sum(diag(invdxx%*%invdxx))
    den=n-2*sd1*(1+lam)*(n+N)/p+sd2*(n+N)^2*(1+lam)^2/p^2
    num=-sum(y^2)+(n+N)*(1+2*lam)*sum(xy*invdxy)/p-lam*(1+lam)*(n+N)^2*as.numeric(t(invdxy)%*%xx%*%invdxy)/p^2
    num=num+n-sd1*(1+2*lam)*(n+N)/p+sd2*(n+N)^2*(1+lam)*lam/p^2
    r2=as.numeric(num/den)
    r2=min(1,max(0,r2))
  }

  r2ba=max(0,r2-r2a)
  print(c(r2,r2a,r2ba))

  #4. Compute permutation p-value for H0: r2=0.

  result=rep(0,npm)
  for(jj in 1:npm){
    if(jj-as.integer(jj/50)*50==1){print(c(jj,jj))}

    samp=sample(n+N)
    XXp=cbind(XX[,1:pa],XX[samp,(pa+1):p])
    xp=XXp[1:n,]
    Xp=XXp[(n+1):(n+N),]

    #r2p=r2a
    #for(ii in 1:1){
    #  if(ii>1 & r2p<1){lam=r2p/(1-r2p)}

    delta=(1+lam*(n+N)/p)*t(xp)%*%xp+t(Xp)%*%Xp
    invd=chol2inv(chol(delta))
    xy=as.vector(t(xp)%*%y)
    invdxy=as.vector(invd%*%xy)
    xx=t(xp)%*%xp
    invdxx=invd%*%xx
    sd1=sum(diag(invdxx))
    sd2=sum(diag(invdxx%*%invdxx))
    den=n-2*sd1*(1+lam)*(n+N)/p+sd2*(n+N)^2*(1+lam)^2/p^2
    num=-sum(y^2)+(n+N)*(1+2*lam)*sum(xy*invdxy)/p-lam*(1+lam)*(n+N)^2*as.numeric(t(invdxy)%*%xx%*%invdxy)/p^2
    num=num+n-sd1*(1+2*lam)*(n+N)/p+sd2*(n+N)^2*(1+lam)*lam/p^2
    r2p=as.numeric(num/den)
    r2p=min(1,max(0,r2p))
    # }
    result[jj]=max(0,r2p-r2a)
  }

  pvalueEST=mean(1.0*(result>r2ba))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(pvalueEST,pvalueBOUND),c(r2,r2a,r2ba),result)

}
