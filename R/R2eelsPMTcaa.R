#' Another conditional permutation test...
#'
#' One paragraph descirbing this function...
#'
#' @param y outcome
#' @param x covariates
#' @param pa dimension of covariates to be adjusted
#' @param npm permutation sample size. Default is 1000
#'
#' @details Details of this function...
#'
#' @return Output of this function...
#'
#' @references Reference 1 here...
#' @references Reference 2 here...
#'
#' @examples \dontrun{R2eelsPMTcaa(y,x,pa,npm=1000)}
#'
#' @export
R2eelsPMTcaa=function(y,x,pa,npm=1000){

  # y==outcome
  # x==covariates
  # npm==permutation sample size
  # pa==dimension of covariates to be adjusted

  n=dim(x)[1]
  p=dim(x)[2]

  #1. Standardization

  for(j in 1:p){
    mux=mean(x[,j])
    sdx=sd(x[,j])
    x[,j]=(x[,j]-mux)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy

  #2. Compute the estimator

  delta=chol2inv(chol(t(x[,1:pa])%*%x[,1:pa]))
  xy=t(x[,1:pa])%*%y
  r2a=1-(n-1-t(xy)%*%delta%*%xy)/(n-pa)
  r2a=min(1,max(0,r2a))

  #W=diag(rep(1,n))-xa%*%chol2inv(chol(t(xa)%*%xa))%*%t(xa)
  #r2a=as.numeric(1-t(y)%*%Wa%*%y/(n-pa))
  delta=chol2inv(chol(t(x)%*%x))
  xy=t(x)%*%y
  r2=1-(n-1-t(xy)%*%delta%*%xy)/(n-p)
  r2=min(1,max(0,r2))

  r2ba=max(0,r2-r2a)

  print(c(r2,r2a,r2ba))
  #4. Compute permutation p-value for H0: r2=0.
  result=rep(0,npm)
  for(jj in 1:npm){
    if(jj-as.integer(jj/50)*50==1){print(c(jj,jj))}

    samp=sample(n)
    xx=cbind(x[,1:pa],x[samp,(pa+1):p])

    delta=chol2inv(chol(t(xx)%*%xx))
    xy=t(xx)%*%y
    r2p=1-(n-1-t(xy)%*%delta%*%xy)/(n-p)
    r2p=min(1,max(r2p,0))

    result[jj]=max(0,r2p-r2a)
  }

  pvalueEST=mean(1.0*(result>r2ba))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(pvalueEST,pvalueBOUND),c(r2,r2a,r2ba),result)

}
