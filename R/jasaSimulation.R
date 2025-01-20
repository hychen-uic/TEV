# upload package TEV
#if(1==1){
 #install.packages("devtools")
 #library(devtools)
 #devtools::install_github("hychen-uic/TEV",force=T)
 #devtools::install_github("xliusufe/RidgeVar")

# library(TEV)
# library(RidgeVar)
#}
#' @import RidgeVar
NULL
#' Simulation studies for JASA submission
#'
#' @param cindep determine if covariates are independent: cindep="T" if independent,
#'                                                        cindep="F" if dependent.
#' @param sparse determine whether effects are sparse: spare="T" if effects are sparse,
#'                                                     sparse="F" if effects are dense
#' @param n  sample size of main data, e.g., n=400.
#' @param N  sample size for supplementary data, e.g., N=200, it is set to N=p in the simulation.
#' @param p dimension of the covariates, for example, p=200,
#' @param p1 subdimension of covariates having nonzero effects p1<=p, e.g. p1=100 when p=200.
#'           in the simulation of sparse effects, p1=4.
#' @param cs the constant for determining regression coefficients. it needs to be adjusted to achieve R2 wanted.
#' @param powx powertansformation of normal covariates, powx=1 means normal, powx=2 means chi-square
#' @param powy powertansformation of normal random error, powx=1 means normal, powx=2 means chi-square
#' @param xsig standard deviation of covariates
#' @param errsig standard deviation of random error.
#' @param method method for dependence generate: method="our" means our approach, otherwise Cai and Guo's approach.
#' @param rho correlation coefficient in the AR(1) covariate dependence model. This has effects only when method!="our".
#' @param nrep number of replicates in simulation sample for variance estimation nrep=1000 by default.
#' @param rept number of replicates in simulation rept=1000 by default.
#'
#' @return Estimate of the proportion of the explained variation and confidence intervals for the proportion.
#'
#' @examples \dontrun{jasaSIMULATION(cindep="T",sparse="F", p=200, p1=100, n=400,N=200,
#'                                   powx=1,powy=1,xsig=1,errsig=1,rho=0.9,rept=100,nrep=1000)}
#'
#' @export
#'
jasaSimulation=function(cindep="T",sparse="F",n=400,N=200,p=200,p1=100, cs=1.0,
                powx=1,powy=1,xsig=1,errsig=1,method="our",rho=0.9,nrep=1000,rept=1000){
# 1. Parameter setup in the simulation study
#cindep="F" # covariate independent ("T") or not ("F")

#sparse="F" # if sparse effects are true, otherwise dense effects
#p=200 # covariate dimension
#n=400 # sample size of the main data
#N=p        # Sample size of the supplementary covariate data

#rept=1000   # simulation replicates
#nrep=1000  # replicates for variance approximation

if(cindep=="T"){
   sqrtsig=diag(rep(1,p))      # independent covariates
  }else{
    if(method=='our'){
      sqrtsig=sqrtpdm(makecora(1.0,0,0,p)[[1]])[[1]] # half of the dependence matrix of covariates
    }else{
      sqrtsig=sqrtpdm(makecorb(rho=0.9,fix=F,p)[[1]])[[1]] #alternative dependence matrix of x
    }
  }

if(sparse=="T"){
 #p1=4    # sparse regression parameters
 sgn=1 #sign(runif(p1)-0.5)
 beta=c(sgn*rep(cs/sqrt(p1),p1),rep(0,p-p1))
}else{
 #p1=p/2    # dense regression parameters
 sgn=1.0 #sign(runif(p1)-0.5)
 beta=c(sgn*rep(cs/sqrt(p1),p1),rep(0,p-p1)) #0.18, 0.35 (corr. normal)                                           #0.26, 0.52 (corr. nonnormal)
 #beta=c(sgn*1.0/c(1:p1), rep(0,p-p1)) #alternative dense effects
}

 #xsig=1.0  # random error parameter
 #errsig=sqrt(1)

# powx=1    # covariate power transformation
# powy=1    # random error power transformation

# 2. Find true r2 by simulation approximation

fit=TEV::trueRV(n*1000,p,beta,xsig,errsig,powx,powy,sqrtsig)
r20=fit[[1]]
v20=fit[[2]]
print(c(r20,v20))
#lambda=r20/(1-r20)
#ilam=lambda # initial value setup
ilam=0.1
iiter=0
palpha=c(0.05) #c(0.1,0.05,0.01)
pa=length(palpha)


# 3. Simulations

R2EP=array(0,c(rept,1+2*pa))  #Eigenprism of Janson et al (2017)

R2EE0=array(0,c(rept,3+4*pa)) #Estimating equation with independent covariates
V2EE0=R2EE0                   #Our proposal

R2EE1=array(0,c(rept,3+4*pa)) #Estimating equation with independent covariates
V2EE1=R2EE1                   #One iteration

R2EE20=array(0,c(rept,3+4*pa)) #Estimating equation with independent covariates
V2EE20=R2EE20                  #20 iterations


R2EET=array(0,c(rept,3+4*pa)) #Estimating equation with independent covariates
V2EET=R2EET                   #Our proposal


R2LS=array(0,c(rept,3+4*pa)) #Least-square (computed only when n>p)
V2LS=R2LS                    #Our proposal with theoretical variance estimate

R2LSa=array(0,c(rept,3+4*pa)) #Least-square (computed only when n>p)
V2LSa=R2LSa                   #Our proposal with simulated variance estimate

R2GRE=array(0,c(rept,4)) # generalized regression estimator of Hou et al (2019)
                         # least-square approach with a different variance estimate
                         # works for n>p only.

R2SD0=array(0,c(rept,3+4*pa)) #Estimating equation with supplementary data
V2SD0=R2SD0                    #Our proposal

R2SD1=array(0,c(rept,3+4*pa)) #Estimating equation with supplementary data
V2SD1=R2SD1                   #one iteration

R2SD20=array(0,c(rept,3+4*pa)) #Estimating equation with supplementary data
V2SD20=R2SD20                  #20 iteration

R2SDT=array(0,c(rept,3+4*pa)) #Estimating equation with supplementary data
V2SDT=R2SD20

R2TS=array(0,c(rept,3+4*pa)) #Transformed data estimating equation approach
V2TS=R2TS                    #An ad hoc approach

R2MLDE=array(0,c(rept,2+2*pa))# MLE from Dicker and Erdogdu (2016)
V2MLDE=R2MLDE                 # use MLE for parameter estimate but spectrum for variance estimate

R2MLRE=R2MLDE        # maximum likelihood for all parameters estimate
V2MLRE=R2MLRE        # full information matrix inversion (Dicker and Erdogdu(2016))

R2EEML=R2MLDE   # estimate r2 by standard outcome,sigma_y=1, GCTA+REML
V2EEML=R2EEML    # variance estimate inverting information for r2 only.

R2CHIVE=array(0,c(rept,2+2*pa)) # Calibration approach of Cai and Guo (2020)
V2CHIVE=R2CHIVE                 # with supplemental covariates

R2CHIVE0=array(0,c(rept,2+2*pa)) # Calibration approach of Cai and Guo (2020)
V2CHIVE0=R2CHIVE0                # without supplemental covariates

R2OT=array(0,c(rept,2)) #Other approaches include
V2OT=array(0,c(rept,2)) #RR(ridge regression,Liu et al, 2020, Biometrika)
                        #MM(Dicker, 2014, Biometrika)

for(i in 1:rept){
  if(i==1){print(c(r20,v20))}
  print(c(i,i,i))

  xy=TEV::makedat(n,p,beta,xsig,errsig,powx,powy,sqrtsig)
  x=xy[[1]]
  y=xy[[2]]

  if(1==1){ # check covariate correlation
    D=cor(x)
    D=D-diag(diag(D))
    #print(c(max(D),min(D)))
    hist(D)
    }


#3.1. EigenPrism approach
    aa=TEV::RVep(y,x,alpha=palpha)
    R2EP[i,1]=aa[[1]]   # estimator
    R2EP[i,2:(1+2*pa)]=aa[[2]] # 99%, 95%, 90% confidence intervals
  #print("ee")
#3.2. Estimating equation approach with direct variance estimation
    aa=TEV::RVee(y,x,lam=ilam,niter=0,alpha=palpha)
    R2EE0[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2EE0[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2EE0[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2EE0[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2EE0[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2EE0[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

    aa=TEV::RVee(y,x,lam=ilam,niter=1,alpha=palpha)
    R2EE1[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2EE1[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2EE1[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2EE1[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2EE1[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2EE1[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

    aa=TEV::RVee(y,x,lam=ilam,niter=20,alpha=palpha)
    R2EE20[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2EE20[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2EE20[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2EE20[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2EE20[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2EE20[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

    aa=TEV::RVee(y,x,lam=r20/(1-r20),niter=50,alpha=palpha)
    R2EET[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2EET[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2EET[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2EET[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2EET[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2EET[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general



#3.3. least-square approach
   if(n>p){ # compute only when n>p
     #print("ls")
     aa=TEV::RVls(y,x,alpha=palpha)
     R2LS[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
     R2LS[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
     R2LS[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
     V2LS[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
     V2LS[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
     V2LS[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

    #Alternative least-square variance estimate

    if(i==1){ # simulation estimate of the variance parts
        nr2=100
        u = rep(1, p)/sqrt(p)
        SUZZU = rep(0, nrep)
        SUZWZU = array(0,c(nrep,nr2)) #rep(0, nrep)
        #SUZWWZU = array(0,c(nrep,nr2)) #rep(0, nrep)
        STRWM = array(0,c(nrep,nr2))#rep(0, nrep)
        for (j in 1:nrep) {
            z = matrix(rnorm(n * p), ncol = p)
            z = zscale(z)[[1]]
            zsvd=svd(z,nu=0)
            vu = t(zsvd$v)%*% u
            SUZZU[j] = sum(vu^2*zsvd$d^2)
            for(k in 1:nr2){
              lam=(k-1)/(nr2-k+1)
              SUZWZU[j,k] = SUZZU[j]*(n-p)*p/(p+lam*n)^2
              #SUZWWZU[j,k] = SUZZU[j]*((n-p)*p)^2/(p+lam*n)^4
              STRWM[j,k] = (n-p)/(p+lam*n)^2
            }
        }

        K1 = var(SUZZU - n)
        kvls=array(0,c(3,nr2))
        for(k in 1:nr2){
          K2 = var(SUZWZU[,k] - STRWM[,k])
          K3 = cov(SUZWZU[,k] - STRWM[,k], SUZZU - n)

          kvls[,k]=c(K1,K2,K3)
       }
    }

    aa=TEV::RVsd(y, x, lam =ilam, alpha=palpha, niter = iiter, know="yes", KV=kvls)
    R2LSa[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2LSa[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2LSa[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2LSa[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2LSa[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2LSa[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

#print("GRE")
#another version of LSE by Hou et al (2019) nature genetics
    fit=TEV::GRE(y,x)
    R2GRE[i,]=c(fit[[1]],fit[[2]],fit[[3]])

     }
#print("sd")
#3.5. Estimating equation approach with supplementary covariates
  if(N>0){ # compute only when supplementary data are available
    xy=TEV::makedat(N,p,beta,xsig,errsig,powx,powy,sqrtsig)
    xsup=xy[[1]]

    if(i==1){ # simulation estimate of the variance parts
        nr2=100 # divided r2 in [0,1] into 100 segments.
        u = rep(1, p)/sqrt(p)
        SUZZU = rep(0, nrep) #array(0,c(nrep,nr2))#rep(0, nrep)
        SUZWZU = array(0,c(nrep,nr2))# rep(0, nrep)
        #SUZWWZU = array(0,c(nrep,nr2))# rep(0, nrep)
        STRWM = array(0,c(nrep,nr2))#rep(0, nrep)
        for (j in 1:nrep) {
            z = matrix(rnorm(n * p), ncol = p)
            z = zscale(z)[[1]]
            zu = z %*% u
            Z = matrix(rnorm(N * p), ncol = p)
            Z = zscale(Z)[[1]]
            SM = z %*% chol2inv(chol(t(z) %*% z + t(Z) %*% Z)) %*%
                t(z) * (n + N)/p

            SUZZU[j] = sum(zu^2)
            SMsvd=svd(SM,nv=0)
            QZU=t(SMsvd$u)%*%zu
            for( k in 1:nr2){
              lam=(k-1)/(nr2-k+1)
              SUZWZU[j,k] =sum(QZU^2*(SMsvd$d-1)/(1+lam*SMsvd$d)^2)
                          #t(zu) %*% SW %*% zu
              #SUZWWZU[j,k] =sum(QZU^2*(SMsvd$d-1)^2/(1+lam*SMsvd$d)^4)
                          #t(SWzu) %*% SWzu
              STRWM[j,k] = sum(SMsvd$d*(SMsvd$d-1)/(1+lam*SMsvd$d)^2)/n
                          # sum(diag(SW %*% SM))/n # is it M or ZZ^t/p?
           }
        }

        K1 = var(SUZZU - n)
        K2=rep(0,nr2)
        K3=rep(0,nr2)

        kv=array(0,c(3,nr2))
        for(k in 1:nr2){
          K2[k] = var(SUZWZU[,k] - STRWM[,k])
          K3[k] = cov(SUZWZU[,k] - STRWM[,k], SUZZU - n)
          kv[,k]=c(K1,K2[k],K3[k])
        }
    }

    aa=TEV::RVsd(y, x, xsup, lam =ilam, alpha=palpha, niter =iiter, know="yes", KV=kv )
    R2SD0[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2SD0[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2SD0[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2SD0[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2SD0[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2SD0[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

    aa=TEV::RVsd(y, x, xsup, lam =ilam, alpha=palpha, niter =1, know="yes", KV=kv )
    R2SD1[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2SD1[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2SD1[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2SD1[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2SD1[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2SD1[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

    aa=TEV::RVsd(y, x, xsup, lam =ilam, alpha=palpha, niter =20, know="yes", KV=kv )
    R2SD20[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2SD20[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2SD20[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2SD20[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2SD20[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2SD20[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

    aa=TEV::RVsd(y, x, xsup, lam =r20/(1-r20), alpha=palpha, niter =0, know="yes", KV=kv )
    R2SDT[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2SDT[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2SDT[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2SDT[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2SDT[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2SDT[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

    aa=TEV::CHIVE(y,x,xsup,alpha=palpha)
    R2CHIVE[i,1:2]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2CHIVE[i,3:(2+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    V2CHIVE[i,1:2]=aa[[3]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2CHIVE[i,3:(2+2*pa)]=aa[[4]]   # 99%, 95%, 90% confidence intervals under normal
}

#3.6. Transformation approach with supplementary covariates
    #transform correlated covariates before applying estimating equation approach
    if(N>0){
        z=transf(rbind(x,xsup))[[1]][1:n,]
     }else{
        z=transf(x)[[1]]
     }
    aa=TEV::RVee(y,z,alpha=palpha,lam=ilam,niter=iiter)
    R2TS[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2TS[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2TS[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2TS[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2TS[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2TS[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

#3.7. DEMLE method of Dicker and Erdogdu (2016)

    aa=TEV::RVmlde(y,x,alpha=palpha)
    R2MLDE[i,1:2]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2MLDE[i,3:(2+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    V2MLDE[i,1:2]=aa[[3]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2MLDE[i,3:(2+2*pa)]=aa[[4]]   # 99%, 95%, 90% confidence intervals under normal

    aa=TEV::FULLREML(y,x,alpha=palpha)
    R2MLRE[i,1:2]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2MLRE[i,3:(2+2*pa)]=aa[[2]]
    V2MLRE[i,1:2]=aa[[3]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2MLRE[i,3:(2+2*pa)]=aa[[4]]   # 99%, 95%, 90% confidence intervals under normal

    aa=TEV::GCTAREML(y,x,alpha=palpha)
    R2EEML[i,1:2]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2EEML[i,3:(2+2*pa)]=aa[[2]]
    V2EEML[i,1:2]=aa[[3]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2EEML[i,3:(2+2*pa)]=aa[[4]]   # 99%, 95%, 90% confidence intervals under normal

#3.8. CHIVE by Cai and Guo (2020) JRSSB with sparsity assumption

    aa=TEV::CHIVE(y,x,alpha=palpha)
    R2CHIVE0[i,1:2]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2CHIVE0[i,3:(2+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    V2CHIVE0[i,1:2]=aa[[3]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2CHIVE0[i,3:(2+2*pa)]=aa[[4]]   # 99%, 95%, 90% confidence intervals under normal

#3.9. Other methods include the ridge regression (RidgeVar::VAR_RR, Liu et al, 2020)
#                            the mle (RidgeVar::VAR_MLE, Dicker and Erdogdu, 2016)
#                            the mm (RidgeVar::VAR_MM, Dicker, 2014)
# The functions all estimate the residual variance in the output $sigma2
# it is converted to explained variation by   var(y)-$sigma2
#   and the explained proportion of variation by    R2=1-$sigma2/var(y)
#
 fit=RidgeVar::VAR_RR(y,x)
  R2OT[i,1]=min(1,max(0,1-fit$sigma2/var(y)))
  V2OT[i,1]=max(0,var(y)-fit$sigma2)
 fit=RidgeVar::VAR_MM(y,x)
  R2OT[i,2]=min(1,max(0,1-fit$sigma2/var(y)))
  V2OT[i,2]=max(0,var(y)-fit$sigma2)

}

# 4. Result output
 print("True R2 and V2")
 print(c(r20,v20))

#4.1. print("Results on the EigenPrism approach from Jansen et al (2017)")
 print("Row 1: estimate and CI; Row 2: empirical variance of Row 1")
 print("Row 3: CI length; Row 4: CI coverage")
 print("EP R2")
 print(apply(R2EP,2,mean))
 print(apply(R2EP,2,var))
 for(k in 1:pa){
   print(mean(R2EP[,1+2*k]-R2EP[,2*k]))
   print(mean(100*(R2EP[,2*k]<=rep(r20,rept))*(R2EP[,1+2*k]>=rep(r20,rept))))
   }

#4.2. print("Results on the estimating equation with permutation variance estimation")
 print("Row 1: estimate, e.var, and CI; Row 2: empirical variance of Row 1;Row 3: CI length; Row 4: CI coverage")
 print("For aditional rows, Row 5: CI length, under normal; Row 6: CI coverage, under normal" )
 print("The rest are similar to this famat.")
 print("EE R2")
 Soutput(R2EE0,pa,r20,rept)
 print("EE V2")
 Soutput(V2EE0,pa,v20,rept)
 print("EE1 R2")
 Soutput(R2EE1,pa,r20,rept)
 print("EE1 V2")
 Soutput(V2EE1,pa,v20,rept)
 print("EE20 R2")
 Soutput(R2EE20,pa,r20,rept)
 print("EE20 V2")
 Soutput(V2EE20,pa,v20,rept)
 print("EET R2")
 Soutput(R2EET,pa,r20,rept)
 print("EET V2")
 Soutput(V2EET,pa,v20,rept)

#4.3. print("Results on the estimating equation with least square weight")
 if(n>p){
   print("LS R2")
   Soutput(R2LS,pa,r20,rept)
   print("LS V2")
   Soutput(V2LS,pa,v20,rept)

   print("LSa R2")
   Soutput(R2LSa,pa,r20,rept)
   print("LSa V2")
   Soutput(V2LSa,pa,v20,rept)

   print("GRE R2")
   print(apply(R2GRE,c(2),mean))
   print(apply(R2GRE,c(2),var))
   print(mean(R2GRE[,4]-R2GRE[,3]))
   print(100*mean((R2GRE[,3]<=as.numeric(r20))*1.0*(R2GRE[,4]>=as.numeric(r20))))
   }

#4.4. print("Results on the estimating equation with supplementary covariate data")
 if(N>0){
   print("SD R2")
   Soutput(R2SD0,pa,r20,rept)
   print("SD V2")
   Soutput(V2SD0,pa,v20,rept)
   print("SD1 R2")
   Soutput(R2SD1,pa,r20,rept)
   print("SD1 V2")
   Soutput(V2SD1,pa,v20,rept)
   print("SD20 R2")
   Soutput(R2SD20,pa,r20,rept)
   print("SD20 V2")
   Soutput(V2SD20,pa,v20,rept)
   print("SDT R2")
   Soutput(R2SDT,pa,r20,rept)
   print("SDT V2")
   Soutput(V2SDT,pa,v20,rept)
   }
#4.5. print("Results on the transformation approach with supplementary covariate data")
 print("Transform R2")
 Soutput(R2TS,pa,r20,rept)
 print("Transform V2")
 Soutput(V2TS,pa,v20,rept)

#4.7  print("Results on the MLE method by Dicker and Erdogdu (2016)")
    print("MLDE R2")
    print(apply(R2MLDE, 2, mean))
    print(apply(R2MLDE, 2, var))
    print(mean(R2MLDE[, 4] - R2MLDE[, 3]))
    print(mean(100 * (R2MLDE[, 3] <= rep(r20, rept)) * (R2MLDE[, 4] >= rep(r20, rept))))

    print("MLDE V2")
    print(apply(V2MLDE, 2, mean))
    print(apply(V2MLDE, 2, var))
    print(mean(V2MLDE[, 4] - V2MLDE[, 3]))
    print(mean(100 * (V2MLDE[, 3] <= rep(v20, rept)) * (V2MLDE[, 4] >= rep(r20, rept))))

    print("MLRE R2")
    print(apply(R2MLRE, 2, mean))
    print(apply(R2MLRE, 2, var))
    print(mean(R2MLRE[, 4] - R2MLRE[, 3]))
    print(mean(100 * (R2MLRE[, 3] <= rep(r20, rept)) * (R2MLRE[, 4] >= rep(r20, rept))))

    print("MLRE V2")
    print(apply(V2MLRE, 2, mean))
    print(apply(V2MLRE, 2, var))
    print(mean(V2MLRE[, 4] - V2MLRE[, 3]))
    print(mean(100 * (V2MLRE[, 3] <= rep(r20, rept)) * (V2MLRE[, 4] >= rep(r20, rept))))

    print("EEML R2")
    print(apply(R2EEML, 2, mean))
    print(apply(R2EEML, 2, var))
    print(mean(R2EEML[, 4] - R2EEML[, 3]))
    print(mean(100 * (R2EEML[, 3] <= rep(r20, rept)) * (R2EEML[, 4] >= rep(r20, rept))))

    print("EEML V2")
    print(apply(V2EEML, 2, mean))
    print(apply(V2EEML, 2, var))
    print(mean(V2EEML[, 4] - V2EEML[, 3]))
    print(mean(100 * (V2EEML[, 3] <= rep(r20, rept)) * (V2EEML[, 4] >= rep(r20, rept))))

#4.8  print("Results on the CHIVE method by Cai and Guo (2020)")
    print("CHIVE R2")
    print(apply(R2CHIVE, 2, mean))
    print(apply(R2CHIVE, 2, var))
    print(mean(R2CHIVE[, 4] - R2CHIVE[, 3]))
    print(mean(100 * (R2CHIVE[, 3] <= rep(r20, rept)) * (R2CHIVE[, 4] >= rep(r20, rept))))

    print("CHIVE V2")
    print(apply(V2CHIVE, 2, mean))
    print(apply(V2CHIVE, 2, var))
    print(mean(V2CHIVE[, 4] - V2CHIVE[, 3]))
    print(mean(100 * (V2CHIVE[, 3] <= rep(r20, rept)) * (V2CHIVE[, 4] >= rep(r20, rept))))

    print("CHIVE0 R2")
    print(apply(R2CHIVE0, 2, mean))
    print(apply(R2CHIVE0, 2, var))
    print(mean(R2CHIVE0[, 4] - R2CHIVE0[, 3]))
    print(mean(100 * (R2CHIVE0[, 3] <= rep(r20, rept)) * (R2CHIVE0[, 4] >= rep(r20, rept))))

    print("CHIVE0 V2")
    print(apply(V2CHIVE0, 2, mean))
    print(apply(V2CHIVE0, 2, var))
    print(mean(V2CHIVE0[, 4] - V2CHIVE0[, 3]))
    print(mean(100 * (V2CHIVE0[, 3] <= rep(r20, rept)) * (V2CHIVE0[, 4] >= rep(r20, rept))))

#4.9 Results on ridge regression approach of Liu et al (2019) from their RidgeVar package
 print("RR and MM:rows 1 and 2: R2 Mean and variance; rows 3 and 4: v2 mean and variance")
 print(apply(R2OT,c(2),mean))
 print(apply(R2OT,c(2),var))
 print(apply(V2OT,c(2),mean))
 print(apply(V2OT,c(2),var))
}
