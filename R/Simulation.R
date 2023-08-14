#' @import graphics
#' @import cccp
NULL
#' Perform simulation comparing methods including
#' EigenPrism
#' RVee
#' RVeels if (n>p)
#' RVeesd if supplementary data are available
#' projection+RVee
#'
#' @param p exposure dimension
#' @param n sample size of original data
#' @param Nsd sample size of the supplementary covariates.
#' @param nrep Monte Carlo sample size for computinging KV.
#' @param rept repetition of simulations
#' @param cindep independet exposures or not
#' @param powx power of normal for x to simulation non-normal covariates
#' @param powy power of normal for y to simulation non-normal outcome
#' @param beta0 non-zero regression coefficient
#' @param p1 number of non-zero regression coefficients
#'
#' @details  perform simulation to compare different methods under varying situations
#'
#' @return simulation results
#'
#' @references Chen, H. Y., Li, H., Argos, M., Persky, V. W., and Turyk, M. (2022). Statistical Methods
#' for Assessing Explained Variation of a Health Outcome by Mixture of Exposures. International Journal
#' of Environmental Research and Public Health.
#' @references Chen H. Y, Zhang, B. and Pan, G.(2023).Estimation and inference on explained variation
#'with possible supplementary data, Manuscript.
#'
#' @examples \dontrun{simulation(p=200,n=400,Nsd=200, rept=1000,nrep=1000,
#'               cindep="T",powx=1,powy=1,beta0=1.5,p1=100)}
#'
#' @export
#
simulation=function(p=200, n=400, Nsd=200, rept=1000,nrep=1000,
                    cindep="T",powx=1,powy=1,beta0=1.5,p1=100){

#cindep="T" # covariate independent ("T") or not ("F")
#fact=1     # rescale factor
#p=200*fact # covariate dimension
#n=400*fact # sample size of the main data
#Nsd=p    # Sample size of the supplementary covariate data

#rept=1000  # simulation replicates
#nrep=1000  # replicates for variance approximation
###bt="F"     # compute bootstrap variance ("T") or not ("F")
###btn=1000   # bootstrap sample size

if(cindep=="T"){
   sqrtsig=diag(rep(1,p))      # independent covariates
  }else{
   sqrtsig=sdgen(2,0,0,p)[[1]] # dependent covariates
   #max(sqrtsig%*%sqrtsig-diag(rep(1,p)))
   #min(sqrtsig%*%sqrtsig)
  }

 #p1=p/2    # regression parameters
 beta=c(rep(beta0/sqrt(p1),p1),rep(0,p-p1))

 xsig=1.0  # random error parameter
 errsig=sqrt(1)

# powx=1    # covariate power transformation
# powy=1    # random error power transformation

# 2. Find true r2 by simulation approximation

fit=trueR2(n*100,p,beta,xsig,errsig,powx,powy,sqrtsig)
r20=fit[[1]][1]
s20=fit[[1]][2]
print(c(r20,s20))
lambda=as.numeric(r20/(1-r20))
ilam=lambda # initial value setup
iiter=1

# 3. Simulations

resultEP=array(0,c(rept,7))  #Eigenprism
resultPS=array(0,c(rept,30)) #Estimating equation
resultLS=array(0,c(rept,30)) #Least-square (n>p)
resultES=array(0,c(rept,30)) #Estimating equation with supplementary data
resultTS=array(0,c(rept,30)) #Transformed data estimating equation approach

for(i in 1:rept){
  print(c(i,i,i))

  xy=datgen(n,p,beta,xsig,errsig,powx,powy,sqrtsig)
  x=xy[[1]]
  y=xy[[2]]

  if(1==1){ # check covariate correlation
    D=cor(x)
    max(D-diag(diag(D)))
    hist(D-diag(diag(D)))
    }

#3.1. EigenPrism approach
    aa=EigenPrismFull(y,x,alpha=c(0.01,0.05,0.1))
    resultEP[i,1]=aa[[1]]   # estimator
    resultEP[i,2:7]=aa[[2]] # 99%, 95%, 90% confidence intervals

#3.2. Estimating equation approach with direct variance estimation
    aa=RVee(y,x,lam=ilam,niter=iiter)

    resultPS[i,1:3]=aa[[1]]   # estimator, variance estimate under normal, variance estimate
    resultPS[i,4:9]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    resultPS[i,10:15]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    resultPS[i,16:18]=aa[[4]]   # estimator, variance estimate under normal, variance estimate
    resultPS[i,19:24]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    resultPS[i,25:30]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

#3.3. least-square approach
   if(n>p){ # compute only when n>p
     aa=RVeels(y,x)
     resultLS[i,1:3]=aa[[1]]   # estimator, variance estimate under normal, variance estimate
     resultLS[i,4:9]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
     resultLS[i,10:15]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
     resultLS[i,16:18]=aa[[4]]   # estimator, variance estimate under normal, variance estimate
     resultLS[i,19:24]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
     resultLS[i,25:30]=aa[[6]] # 99%, 95%, 90% confidence intervals in general
     }

#3.4. Estimating equation approach with supplementary covariates
  if(Nsd>0){ # compute only when supplementary data are available
  if(i==1){ # simulation estimate of the variance parts
    u=rep(1,p)/sqrt(p)
    SUZWZU=rep(0,nrep)
    SUZWWZU=rep(0,nrep)
    SUZZU=rep(0,nrep)
    STRWM=rep(0,nrep)
    #STRW=rep(0,nrep)
    for(j in 1:nrep){
      z=matrix(rnorm(n*p),ncol=p)
      z=zscale(z)[[1]]
      zu=z%*%u
      Z=matrix(rnorm(Nsd*p),ncol=p)
      Z=zscale(Z)[[1]]
      SM=z%*%chol2inv(chol(t(z)%*%z+t(Z)%*%Z))%*%t(z)*(n+Nsd)/p
      ISM=chol2inv(chol(diag(rep(1,n))+ilam*SM))
      SW=ISM%*%(SM-diag(rep(1,n)))%*%ISM
      SWzu=SW%*%zu

      SUZZU[j]=sum(zu^2)
      SUZWZU[j]=t(zu)%*%SW%*%zu
      SUZWWZU[j]=t(SWzu)%*%SWzu
      #STRW[j]=sum(diag(SW))/n
      STRWM[j]=sum(diag(SW%*%SM))/n
    }
    KV=rep(0,3)
    KV[1]=var(SUZWZU-STRWM)
    KV[2]=cov(SUZWZU-STRWM,SUZZU-n)
    KV[3]=var(SUZZU-n)
    }

    xy=datgen(Nsd,p,beta,xsig,errsig,powx,powy,sqrtsig)
    X=xy[[1]]

    aa=RVeesd(y,x,X,lam=ilam,niter=iiter,KV=KV,know="yes",nrep=10000)
    resultES[i,1:3]=aa[[1]]    #  estimator, variance estimate under normal, variance estimate
    resultES[i,4:9]=aa[[2]]    # 99%, 95%, 90% confidence intervals under normal
    resultES[i,10:15]=aa[[3]]  # 99%, 95%, 90% confidence intervals in general
    resultES[i,16:18]=aa[[4]]  # estimator, variance estimate under normal, variance estimate
    resultES[i,19:24]=aa[[5]]  # 99%, 95%, 90% confidence intervals under normal
    resultES[i,25:30]=aa[[6]]  # 99%, 95%, 90% confidence intervals in general
   }

#3.5. Transformation approach with supplementary covariates
    #transform correlated covariates before applying estimating equation approach
    if(Nsd>0){
        z=transf(rbind(x,X))[[1]][1:n,]
     }else{
        z=transf(x)[[1]]
     }
    aa=RVee(y,z,lam=ilam,niter=5)
    resultTS[i,1:3]=aa[[1]]
    resultTS[i,4:9]=aa[[2]]
    resultTS[i,10:15]=aa[[3]]
    resultTS[i,16:18]=aa[[4]]   # estimator, variance estimate under normal, variance estimate
    resultTS[i,19:24]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    resultTS[i,25:30]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

}

# 4. Result output
 c(r20,s20)

#4.1. print("Results on the EigenPrism approach")
 print("Results on the EigenPrism approach")
 print(apply(resultEP,2,mean))
 print(apply(resultEP,2,var))
 print(mean(resultEP[,3]-resultEP[,2]))
 print(mean(resultEP[,5]-resultEP[,4]))
 print(mean(resultEP[,7]-resultEP[,6]))

 print(mean(100*(resultEP[,2]<=rep(r20,rept))*(resultEP[,3]>=rep(r20,rept))))
 print(mean(100*(resultEP[,4]<=rep(r20,rept))*(resultEP[,5]>=rep(r20,rept))))
 print(mean(100*(resultEP[,6]<=rep(r20,rept))*(resultEP[,7]>=rep(r20,rept))))

#4.2. print("Results on the estimating equation assuming independent exposures")
 print("Under indepedent exposure assumption: proportion of explained variation.")
 Soutput(resultPS[,1:15],r20,rept)
 print("Under indepedent exposure assumption: explained variation.")
 Soutput(resultPS[,16:30],s20,rept)

#4.3. print("Results on the least-square approach")
 if(n>p){
   print("Least-square: proportion of explained variation.")
   Soutput(resultLS[,1:15],r20,rept)
   print("Least-square: explained variation.")
   Soutput(resultLS[,16:30],s20,rept)
   }
#4.4. print("Results on the estimating equation with supplementary covariate data")
 if(Nsd>0){
   print("For suplementary data: proportion of explained variation.")
   Soutput(resultES[,1:15],r20,rept)
   print("For suplementary data: explained variation.")
   Soutput(resultES[,16:30],s20,rept)
   }
#4.5. print("Results on the transformation approach with supplementary covariate data")
   Soutput(resultTS[,1:15],r20,rept)
   Soutput(resultTS[,16:30],s20,rept)

#importFrom("grDevices", "dev.new")
#dev.new()
 hist(resultLS[,3])
 hist(resultES[,3],add=T,col='red')
 hist(resultPS[,3],add=T,col='yellow')

}
