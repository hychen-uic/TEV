#' @import stats
NULL

#' Computes the confidence interval with the estimating equation approach
#'
#' One paragraph describing this function...
#'
#' @param y outcome
#' @param x covariates
#' @param lam parameter adjusting the format of the weighting matrix. Default is 0.1
#' @param niter number of iterations for updating lambda
#'
#' @details details about this function...
#'
#' @return The output is the confidence interval...
#'
#' @references reference 1 here...
#' @references reference 2 here...
#'
#' @examples \dontrun{R2ee(y,x,lam=0.1,niter=3)}
#'
#'@export
R2ee=function(y,x,lam=0.1,niter=3){

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

  #2. Sigular value decomposition
  Xsvd=svd(x,nv=0) #nv=0 means not computing v matrix
  # singular value decomposition
  # $u%*%diag($d)%*%t($v)=X, t($u)%*%$u=I, t($v)%*%$v=I
  Mev=Xsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.

  r2=lam/(1+lam) # initial value
  for(ii in 1:niter){ #iteration to update lambda
    if(ii>1 & r2<1){lam=r2/(1-r2)}
    Dev=(Mev-1)/(1+lam*Mev)^2  #vector of eigenvalues of weight matrix
    #3. Compute the estimators
    uy=t(Xsvd$u)%*%y
    u1=t(Xsvd$u)%*%rep(1,n)
    if(n>=p){
      com=sum(u1^2*(Dev+1))/n-1
      num=sum(uy^2*(Dev+1))-sum(y^2)-sum(Dev)+n-p
      den=sum(Dev*(Mev-1))+n-p
    }else{
      com=sum(u1^2*Dev)/n
      num=sum(uy^2*Dev)-sum(Dev)
      den=sum(Dev*(Mev-1))
    }

    r2=(num+com)/(den+com)
    r2=min(1,max(0,r2))
  }
  #4. estimate variance
  #4(a). With normal random error assumption

  Wev=Mev*Dev
  if(n>=p){
    W=Xsvd$u%*%diag(Dev+1)%*%t(Xsvd$u)-diag(rep(1,n))
    vest=2*p*var(Wev)*r2^2
    S1=sum(Dev^2*Mev)
    S=sum(Dev^2)+n-p
    vest=vest+4*r2*(1-r2)*S1
    vest=vest+2*(1-r2)^2*S
  }else{
    W=Xsvd$u%*%diag(Dev)%*%t(Xsvd$u)
    vest=2*n*var(Wev)*r2^2 #because the length of Wev is n now.
    S1=sum(Dev^2*Mev)
    S=sum(Dev^2)
    vest=vest+4*r2*(1-r2)*S1+2*(1-r2)^2*S
  }

  #4(b). Without normal random error assumption

  M=Xsvd$u%*%diag(Mev)%*%t(Xsvd$u)
  T=sum(diag(W^2))
  veps2=mean((y^2-1-(diag(M)-1)*r2)^2)-4*r2*(1-r2)-2.0*r2^2
  vest1=vest+T*(max(veps2,0)-2*(1-r2)^2)

  vest1=vest1/den^2
  vest=vest/den^2

  ci=r2+sqrt(vest)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci1=r2+sqrt(vest1)*qnorm(c(0.005,0.995,0.025,0.975,0.05,0.95))
  ci[2*c(1:3)-1]=ci[2*c(1:3)-1]*(ci[2*c(1:3)-1]>0)
  ci[2*c(1:3)]=ci[2*c(1:3)]*(ci[2*c(1:3)]<1)+1.0*(ci[2*c(1:3)]>=1)
  ci1[2*c(1:3)-1]=ci1[2*c(1:3)-1]*(ci1[2*c(1:3)-1]>0)
  ci1[2*c(1:3)]=ci1[2*c(1:3)]*(ci1[2*c(1:3)]<1)+1.0*(ci1[2*c(1:3)]>=1)

  #5. output result

  list(c(r2,vest,vest1),ci,ci1)
  #[[1]]==> estimate, estimated variance (under normal), estimated variance.
  #[[2]]==>confidence interval under normal
  #[[3]]==>confidence interval without normal
}
