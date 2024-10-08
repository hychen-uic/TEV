#' The GCTA approach with bootstrap variance estimate
#'
#' This function implements the genetic complex trait
#' analysis of Yang et al (2010) approach assuming independent
#' covariates and normally distributed errors, with bootstrap variance
#' estimate of Schweiger et al (2016).
#'
#' @param y outcome: a vector of length n
#' @param x covariates: a matrix of nxp dimension
#' @param alpha a vector of type I error for create the confidence intervals.
#' @param niter number of iterations
#' @param bt variable specifying whether to compute bootstrap variance. Default is FALSE
#' @param nbt bootstrap sample size
#'
#' @details The function uses the singular value decomposition for estimation and bootstrap
#' sampling approach for variance estimation under normal random errors.
#'
#' @return Estimate of proportion of the explained variation, and bootstrap estimate and
#' variance and confidence intervals if bt=T
#'
#' @references Schweiger, R., Kaufman, S., Laaksonen, R., Kleber, M. E., Marz, W., Eskin, E., Rosset, S., Halperin, E. (2016).
#' Fats and ac-curate construction of confidence intervals for heritability.
#' *The American Journal of Human Genetics*, **98**, 1181-1192.
#' @references Yang, J., Lee, S. H., Wray, N. R., Goddard, M. E., Visscher, P. (2016).
#' GCTA-GREML accounts for linkage disequilibrium when estimating genetic variance from genome-wide SNPs.
#' *Proceedings of the National Academy of Sciences*, **113**, E4579-E4580.
#'
#' @examples \dontrun{R2GCTA(y, x)}
#'
#' @export
R2GCTA=function(y,x,alpha=c(0.05),niter=10,bt=NULL,nbt=1000){
  # Fast estimate r^2 without additional weighting
  # y==outcome
  # x==covariates
  # lam==parameter adjusting the format of weighting matrix

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
  uy=t(Xsvd$u)%*%y
  u1=t(Xsvd$u)%*%rep(1,n)

  #3. compute estimates
  lam=1;r2=0.5 # initial value
  for (iter in 1:niter){ # Perform iteration for updating lambda

    if(iter>1 & r2<1){lam=r2/(1-r2)} # update penalty parameter
    Dev=(Mev-1)/(1+lam*Mev)^2  #vector of eigenvalues of weight matrix

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

  } #end iteration

  est=c(r2*sdy^2,(1-r2)*sdy^2,r2)


  #4. bootstrap variance estimation using normal approximation
  if(bt==TRUE){
    btresult=array(0,c(nbt,3))

    for(ii in 1:nbt){

      if((ii-as.integer(ii/50)*50)==0){print(rep(ii,3))}
      #4.1. resample X
      btx=array(0,c(n,p))
      for(j in 1:p){
        samp=sample(n)
        btx[,j]=x[samp,j]
      }
      #4.2. compute singular value decomposition
      btXsvd=svd(btx,nv=0) # do not compute v matrix
      btMev=btXsvd$d^2/p

      btu1=t(btXsvd$u)%*%rep(1,n)
      btlam=lam

      #4.3. Simulate error parametric bootstrap sample
      bty0=rnorm(n)  #preliminary outcome
      btuy0=t(btXsvd$u)%*%bty0

      bty=btXsvd$u%*%((sqrt(btMev*r2+(1-r2))-sqrt(1-r2))*btuy0)+bty0*sqrt(1-r2)
      btsdy=sd(bty)
      bty=(bty-mean(bty))/btsdy
      btuy=t(btXsvd$u)%*%bty

      for(iter in 1:3){
        btDev=(btMev-1)/(1+btlam*btMev)^2
        if(n>=p){
          btcom=sum(btu1^2*(btDev+1))/n-1
          btnum=sum(btuy^2*(btDev+1))-sum(bty^2)-sum(btDev)+n-p
          btden=sum(btDev*(btMev-1))+n-p
        }else{
          btcom=sum(btu1^2*btDev)/n
          btnum=sum(btuy^2*btDev)-sum(btDev)
          btden=sum(btDev*(btMev-1))
        }
        btr2=(btnum+btcom)/(btden+btcom)
        if(btr2>1){btr2=1}
        if(btr2<0){btr2=0}
        btlam=btr2/(1-btr2)
      }

      btresult[ii,]=c(btr2,btr2*btsdy^2,(1-btr2)*btsdy^2)

      #print(c(ii,btresult[ii,]))

    }

    ord=order(btresult[,1])
    cir2=btresult[c(ord[as.integer(nbt*alpha/2)],ord[as.integer(nbt*(1-alpha/2))]),3]

    ord=order(btresult[,2])
    cis2=btresult[c(ord[as.integer(nbt*alpha/2)],ord[as.integer(nbt*(1-alpha/2))]),3]

    #ciL95=btresult[ord[as.integer(btn*0.025)],3]
    #ciU95=btresult[ord[as.integer(btn*0.975)],3]
    #ci=c(ciL95,ciU95)
    btest=apply(btresult,2,mean)
    btvest=apply(btresult,2,var)
  }

  #5. output result

  if(bt=='T'){
    list(est,btest,btvest,cir2,cis2)
  }else{
    list(est)
  }
  #[[1]]==> estimate
  #[[2]]==> bootstrap estimate
  #[[3]]==> bootstrap variance estimate
  #[[4]]==> confidence interval for r2
  #[[5]]==> confidence interval for s2
}
