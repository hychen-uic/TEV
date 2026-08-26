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
#' @param p dimension of the covariates, for example, p=200,
#' @param p1 subdimension of covariates having nonzero effects p1<=p, e.g. p1=100 when p=200.
#'           in the simulation of sparse effects, p1=4.
#' @param cs the constant for determining regression coefficients. it needs to be adjusted to achieve R2 wanted.
#' @param powx powertansformation of normal covariates, powx=1 means normal, powx=2 means chi-square
#' @param powy powertansformation of normal random error, powx=1 means normal, powx=2 means chi-square
#' @param xsig standard deviation of covariates
#' @param errsig standard deviation of random error.
#' @param covrate rate of combination for two different covariance generators (to create vary correlation structure)
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
NGsimulation=function(cindep="T",sparse="F",n=400,p=200,p1=100, cs=1.0,
                powx=1,powy=1,xsig=1,errsig=1,covrate=0.5,rho=0.9,nrep=1000,rept=1000){
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
  sqrtsig1=sqrtpdm(makecora(1.0,0.01,0.3,p)[[1]])[[1]] # our dependence matrix of x
  sqrtsig2=sqrtpdm(makecorb(rho=rho,fix=F,p)[[1]])[[1]] #alternative dependence matrix of x
  sqrtsig=sqrtsig1*covrate+sqrtsig2*(1-covrate) # combine two covariance matrices
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

R2EE0=array(0,c(rept,3+4*pa)) #Estimating equation with independent covariates
V2EE0=R2EE0                   #Our proposal

R2EE1=array(0,c(rept,3+4*pa)) #Estimating equation with independent covariates
V2EE1=R2EE1                   #One iteration

R2EE20=array(0,c(rept,3+4*pa)) #Estimating equation with independent covariates
V2EE20=R2EE20                  #20 iterations

R2GRE=array(0,c(rept,4)) # generalized random effects estimator of Hou et al (2019)
R2LS=array(0,c(rept,3+4*pa))  # least-square approach with a different variance estimate
V2LS=array(0,c(rept,3+4*pa))
R2LSa=array(0,c(rept,3+4*pa))  # works for n>p only
V2LSa=array(0,c(rept,3+4*pa))

R2TR=array(0,c(rept,3+4*pa))#Estimating equation with independent covariates
V2TR=R2TR                   #transformation by true covariance matrix

R2TS=array(0,c(rept,3+4*pa)) #Sparse precision matrix estimate
V2TS=R2TS

R2TS10=array(0,c(rept,3+4*pa)) #Transformed data estimating equation approach
V2TS10=R2TS10                 # Size of block=0.1*n

R2TS20=array(0,c(rept,3+4*pa)) #Transformed data estimating equation approach
V2TS20=R2TS20                  # Size of block=0.2*n

R2TS30=array(0,c(rept,3+4*pa)) #Transformed data estimating equation approach
V2TS30=R2TS30                  # Size of block=0.3*n

R2TS40=array(0,c(rept,3+4*pa)) #Transformed data estimating equation approach
V2TS40=R2TS40                  # Size of block=0.4*n

R2TS50=array(0,c(rept,3+4*pa)) #Transformed data estimating equation approach
V2TS50=R2TS50                  # Size of block=0.5*n

R2CHIVE=array(0,c(rept,2+2*pa)) # Calibration approach of Cai and Guo (2020)
V2CHIVE=R2CHIVE                 # with supplemental covariates

for(i in 1:rept){
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


#3.1. Estimating equation approach with direct variance estimation
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


#3.2. least-square approach
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

  }# end of least-square

#3.3. Transformation approach with supplementary covariates
    #transform correlated covariates before applying estimating equation approach
    z=x%*%solve(sqrtsig)  #transformation by true covariance matrix
    aa=TEV::RVee(y,z,alpha=palpha,lam=ilam,niter=1)
    R2TR[i,1:3]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2TR[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    R2TR[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general
    V2TR[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2TR[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2TR[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general

if(1==1){ # transformation by latent variables and estimated sparse precision matrix

    maxpc=10 # maximum number of PCs to be identified
    svdx=svd(x)
    # print(svdx$d[1:10])
    npc=0
    for(ijk in 1:maxpc){
      if(svdx$d[ijk]>2*svdx$d[ijk+1]-svdx$d[ijk+10]){npc=ijk}
    }
    # print(npc)

    if(npc>0){
      pc=x%*%svdx$v[,1:npc]
      fit=lm(y~pc)
      yres=as.numeric(fit$residual)
      exvar=var(y-yres)/var(y) #proportion of variation explained by all pcs.
      coef=as.vector(fit$coefficients)[-1]
      residvar=(summary(fit)$sigma)^2/n
      vexvar=2*t(coef)%*%(t(pc)%*%pc/n)%*%coef*residvar/(var(y))^2/n # need to be updated
      # since diag(a number)=identity matrix with that dimension, it has to be cautious
      x=x-as.matrix(svdx$u[,1:npc])%*%diag(svdx$d)[1:npc,1:npc]%*%t(as.matrix(svdx$v[,1:npc]))
      # remove the PC components from x
    }else{
      yres=y
      exvar=0
      vexvar=0
    }


  # identify sparse precision matrix

    fit=SILGGM(x,alpha=palpha,global=TRUE)

    #Omega=fit$precision # precision matrix

    network=fit$global_decision[[1]]
    # print(sum(network))
    #Omega=nodewisereg(dat=x, fixstruct=network)
    Omega=lezhong(dat=x, network=network)

    sqrtOmega=sqrtpdm(Omega)[[1]]

    z=x%*%sqrtOmega
    aa=TEV::RVee(yres,z,alpha=palpha,lam=ilam,niter=1)
    R2TS[i,1]=aa[[1]][1]*(1-exvar)+exvar   # estimator of R2, variance estimate under normal, variance estimate
    R2TS[i,2]=aa[[1]][2]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2
    R2TS[i,3]=aa[[1]][3]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2

    R2TS[i,4]=R2TS[i,1]-1.96*sqrt(R2TS[i,2])
    R2TS[i,5]=R2TS[i,1]+1.96*sqrt(R2TS[i,2])
    R2TS[i,6]=R2TS[i,1]-1.96*sqrt(R2TS[i,3])
    R2TS[i,7]=R2TS[i,1]+1.96*sqrt(R2TS[i,3])

    #R2TS[i,4:(3+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    #R2TS[i,(4+2*pa):(3+4*pa)]=aa[[3]] # 99%, 95%, 90% confidence intervals in general

    V2TS[i,1:3]=aa[[4]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2TS[i,4:(3+2*pa)]=aa[[5]]   # 99%, 95%, 90% confidence intervals under normal
    V2TS[i,(4+2*pa):(3+4*pa)]=aa[[6]] # 99%, 95%, 90% confidence intervals in general
    # print(exvar)
}


    z=decor(dat=x,decorate=0.1,inter=1) #transformed by 10%n block covariance
    aa=TEV::RVee(yres,z,alpha=palpha,lam=ilam,niter=1)
    R2TS10[i,1]=aa[[1]][1]*(1-exvar)+exvar   # estimator of R2, variance estimate under normal, variance estimate
    R2TS10[i,2]=aa[[1]][2]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2
    R2TS10[i,3]=aa[[1]][3]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2

    R2TS10[i,4]=R2TS10[i,1]-1.96*sqrt(R2TS10[i,2])
    R2TS10[i,5]=R2TS10[i,1]+1.96*sqrt(R2TS10[i,2])
    R2TS10[i,6]=R2TS10[i,1]-1.96*sqrt(R2TS10[i,3])
    R2TS10[i,7]=R2TS10[i,1]+1.96*sqrt(R2TS10[i,3])

    z=decor(dat=x,decorate=0.2,inter=1) #transformed by 20%n block covariance
    aa=TEV::RVee(yres,z,alpha=palpha,lam=ilam,niter=1)
    R2TS20[i,1]=aa[[1]][1]*(1-exvar)+exvar   # estimator of R2, variance estimate under normal, variance estimate
    R2TS20[i,2]=aa[[1]][2]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2
    R2TS20[i,3]=aa[[1]][3]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2

    R2TS20[i,4]=R2TS20[i,1]-1.96*sqrt(R2TS20[i,2])
    R2TS20[i,5]=R2TS20[i,1]+1.96*sqrt(R2TS20[i,2])
    R2TS20[i,6]=R2TS20[i,1]-1.96*sqrt(R2TS20[i,3])
    R2TS20[i,7]=R2TS20[i,1]+1.96*sqrt(R2TS20[i,3])

    z=decor(dat=x,decorate=0.3,inter=1) #transformed by 30%n block covariance
    aa=TEV::RVee(yres,z,alpha=palpha,lam=ilam,niter=1)
    R2TS30[i,1]=aa[[1]][1]*(1-exvar)+exvar   # estimator of R2, variance estimate under normal, variance estimate
    R2TS30[i,2]=aa[[1]][2]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2
    R2TS30[i,3]=aa[[1]][3]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2

    R2TS30[i,4]=R2TS30[i,1]-1.96*sqrt(R2TS30[i,2])
    R2TS30[i,5]=R2TS30[i,1]+1.96*sqrt(R2TS30[i,2])
    R2TS30[i,6]=R2TS30[i,1]-1.96*sqrt(R2TS30[i,3])
    R2TS30[i,7]=R2TS30[i,1]+1.96*sqrt(R2TS30[i,3])

    z=decor(dat=x,decorate=0.4,inter=1) #transformed by 40%n block covariance
    aa=TEV::RVee(yres,z,alpha=palpha,lam=ilam,niter=1)
    R2TS40[i,1]=aa[[1]][1]*(1-exvar)+exvar   # estimator of R2, variance estimate under normal, variance estimate
    R2TS40[i,2]=aa[[1]][2]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2
    R2TS40[i,3]=aa[[1]][3]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2

    R2TS40[i,4]=R2TS40[i,1]-1.96*sqrt(R2TS40[i,2])
    R2TS40[i,5]=R2TS40[i,1]+1.96*sqrt(R2TS40[i,2])
    R2TS40[i,6]=R2TS40[i,1]-1.96*sqrt(R2TS40[i,3])
    R2TS40[i,7]=R2TS40[i,1]+1.96*sqrt(R2TS40[i,3])

    z=decor(dat=x,decorate=0.5,inter=1) #transformed by 50%n block covariance
    aa=TEV::RVee(yres,z,alpha=palpha,lam=ilam,niter=1)
    R2TS50[i,1]=aa[[1]][1]*(1-exvar)+exvar   # estimator of R2, variance estimate under normal, variance estimate
    R2TS50[i,2]=aa[[1]][2]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2
    R2TS50[i,3]=aa[[1]][3]*(1-exvar)^2+vexvar*(1-aa[[1]][1])^2

    R2TS50[i,4]=R2TS50[i,1]-1.96*sqrt(R2TS50[i,2])
    R2TS50[i,5]=R2TS50[i,1]+1.96*sqrt(R2TS50[i,2])
    R2TS50[i,6]=R2TS50[i,1]-1.96*sqrt(R2TS50[i,3])
    R2TS50[i,7]=R2TS50[i,1]+1.96*sqrt(R2TS50[i,3])

#3.4. CHIVE by Cai and Guo (2020) JRSSB with sparsity assumption

    aa=TEV::CHIVE(y,x,alpha=palpha)
    R2CHIVE[i,1:2]=aa[[1]]   # estimator of R2, variance estimate under normal, variance estimate
    R2CHIVE[i,3:(2+2*pa)]=aa[[2]]   # 99%, 95%, 90% confidence intervals under normal
    V2CHIVE[i,1:2]=aa[[3]]   # estimatorof V2, variance estimate under normal, variance estimate
    V2CHIVE[i,3:(2+2*pa)]=aa[[4]]   # 99%, 95%, 90% confidence intervals under normal

}

# 4. Result output
 print("True R2 and V2")
 print(c(r20,v20))

#4.1. print("Results on the estimating equation with permutation variance estimation")
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
 print("True covariance, R2")
 Soutput(R2TR,pa,r20,rept)
 print("V2")
 Soutput(V2TR,pa,v20,rept)

#4.2. print("Results on the estimating equation with least square weight")
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

#4.3  print("Results on the sparse precision matrix decorrelation")

 print("R2 from Transformation (sparse)")
 Soutput(R2TS,pa,r20,rept)
 print("V2")
 Soutput(V2TS,pa,v20,rept)


#4.4  print("Results on the simple block wise decorrelation")

 print("R2 from 10%n block-wise transformation")
 Soutput(R2TS10,pa,r20,rept)
 print("V2")
 Soutput(V2TS10,pa,v20,rept)

 print("R2 from 20%n block-wise transformation")
 Soutput(R2TS20,pa,r20,rept)
 print("EE V2")
 Soutput(V2TS20,pa,v20,rept)

 print("R2 from 30%n block-wise transformation")
 Soutput(R2TS30,pa,r20,rept)
 print("V2")
 Soutput(V2TS30,pa,v20,rept)

 print("R2 from 40%n block-wise transformation")
 Soutput(R2TS40,pa,r20,rept)
 print("EE V2")
 Soutput(V2TS40,pa,v20,rept)

 print("R2 from 50%n block-wise transformation")
 Soutput(R2TS50,pa,r20,rept)
 print("V2")
 Soutput(V2TS50,pa,v20,rept)


#4.4  print("Results on the CHIVE method by Cai and Guo (2020)")
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
}
