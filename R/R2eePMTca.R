#' Conditional permutation test adjusting for covariates for estimating equation approach
#'
#' This function performs test of no extra variation explained by the additional covariates
#' by permuting the additional covariates using the estimating equation approach.
#'
#' @param y outcome
#' @param x covariates
#' @param pa dimension of covariates to be adjusted
#' @param lam parameter adjusting the format of the weighting matrix. Default is 0.2
#' @param niter number of iteration for updating lambda. Default is 3
#' @param npm Monte Carlo sample size for permutation. Default is 1000
#'
#' @details The computation permutes the second part of the covariates and is computationally slow.
#'
#' @return Output includes the p-values (estimate and bound) of the test,
#' estimate of proportion of the extra explained variation, and simulation results.
#'
#' @references Chen, H.Y.; Li, H.; Argos, M.; Persky, V.; Turyk, M.
#' Statistical methods for assessing explained variations of a health outcome by mixtures of exposures.
#' Under review for Prep. Spec. Issue Int. J. Environ. Res. Public Health 2022.
#'
#' @examples \dontrun{R2eePMTca(y, x, pa, lam = 0.2, niter = 3, npm = 1000)}
#'
#' @export
R2eePMTca=function(y,x,pa,lam=0.2,niter=3,npm=1000){

  #1. Standardization

  n=dim(x)[1]
  p=dim(x)[2]
  for(j in 1:p){
    mux=mean(c(x[,j]))
    sdx=sd(c(x[,j]))
    x[,j]=(x[,j]-mux)/sdx
  }
  sdy=sd(y)
  y=(y-mean(y))/sdy


  # for xa
  xa=x[,1:pa]
  xasvd=svd(xa,nv=0) #nv=0 means not computing v matrix
  Mev=xasvd$d^2/pa #Vector of eigenvalues of matrix XX'/p.
  r2a=lam/(1+lam) # initial value
  for(ii in 1:niter){ #iteration to update lambda
    if(ii>1 & r2a<1){lam=r2a/(1-r2a)}
    Dev=(Mev-1)/(1+lam*Mev)^2  #vector of eigenvalues of weight matrix

    uy=t(xasvd$u)%*%y
    u1=t(xasvd$u)%*%rep(1,n)
    if(n>=pa){
      com=sum(u1^2*(Dev+1))/n-1
      num=sum(uy^2*(Dev+1))-sum(y^2)-sum(Dev)+n-pa
      den=sum(Dev*(Mev-1))+n-pa
    }else{
      com=sum(u1^2*Dev)/n
      num=sum(uy^2*Dev)-sum(Dev)
      den=sum(Dev*(Mev-1))
    }

    r2a=(num+com)/(den+com)
    r2a=min(1,max(0,r2a))
    #print(r2a)
  }

  xsvd=svd(x,nv=0) #nv=0 means not computing v matrix
  # singular value decomposition
  # $u%*%diag($d)%*%t($v)=x, t($u)%*%$u=I, t($v)%*%$v=I
  Mev=xsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.
  r2=lam/(1+lam) # initial value
  for(ii in 1:niter){ #iteration to update lambda
    if(ii>1 & r2<1){lam=r2/(1-r2)}
    Dev=(Mev-1)/(1+lam*Mev)^2  #vector of eigenvalues of weight matrix

    uy=t(xsvd$u)%*%y
    u1=t(xsvd$u)%*%rep(1,n)
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
    # print(r2)
  }

  r2ba=max(0,r2-r2a)

  print(c(r2,r2a,r2ba))


  #4. Compute permutation p-value for H0: r2=0.

  result=rep(0,npm)
  for(jj in 1:npm){
    if(jj-as.integer(jj/50)*50==1){print(c(jj,jj))}

    samp=sample(n)
    xp=cbind(x[,1:pa],x[samp,(pa+1):p])

    xpsvd=svd(xp,nv=0) #nv=0 means not computing v matrix
    Mev=xpsvd$d^2/p #Vector of eigenvalues of matrix XX'/p.
    #r2p=r2 # initial value
    #for(ii in 1:niter){ #iteration to update lambda
    #  if(ii>1 & r2p<1){lam=r2p/(1-r2p)}
    Dev=(Mev-1)/(1+lam*Mev)^2  #vector of eigenvalues of weight matrix

    uy=t(xpsvd$u)%*%y
    u1=t(xpsvd$u)%*%rep(1,n)
    if(n>=p){
      com=sum(u1^2*(Dev+1))/n-1
      num=sum(uy^2*(Dev+1))-sum(y^2)-sum(Dev)+n-p
      den=sum(Dev*(Mev-1))+n-p
    }else{
      com=sum(u1^2*Dev)/n
      num=sum(uy^2*Dev)-sum(Dev)
      den=sum(Dev*(Mev-1))
    }

    r2p=(num+com)/(den+com)
    r2p=min(1,max(0,r2p))
    #}

    result[jj]=max(0,r2p-r2a)
  }

  pvalueEST=mean(1.0*(result>r2ba))
  acc=2*sqrt(pvalueEST*(1-pvalueEST)/npm)
  crt=2*sqrt(0.05*0.95/npm) # if truth p-value=0.05, how accurate the estimation is
  pvalueBOUND=max(pvalueEST+acc,pvalueEST+crt)

  list(c(pvalueEST,pvalueBOUND),c(r2,r2a,r2ba),result)

}
