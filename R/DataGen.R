datgen=function(n,p,beta,xsig,errsig,powx,powy,sqrtsig){
  #generate data from linear model
  # powx, powy parameters making non-normal x and random error
  # sd, square-root of correlation matrix, making correlated x.

  x=array(0,c(n,p))
  for(j in 1:p){
    x[,j]=(rnorm(n,mean=0,sd=xsig))
    }
   x=sign(x)*abs(x)^powx # simulate correlated covariates
   if(powx==2){
     x=abs(x)
     }
   for(j in 1:p){ #standardization
     mu=mean(x[,j])
     sdx=sd(x[,j])
     x[,j]=(x[,j]-mu)/sdx
     }

   x=x%*%sqrtsig  # create correlated covariates
   x=x*(x<200)+200*(x>=200)


   err=rnorm(n,mean=0,sd=errsig)
   err=sign(err)*abs(err)^powy
   if(powy==2){
     err=abs(err)
     }
   err=err*(err<200)+200*(err>=200)

   y=x%*%beta+err
   y=as.numeric(y)

  list(x,y)
}

sdgen=function(amply,tilt,shrink,p){
   # generate a square root of a correlation matrix
   # amply,tilt,and shrink(>=0) are means of controling the correlation
   # p is the dimension of the matrix

   x=matrix(rnorm(p*p,mean=amply,sd=1),ncol=p)
   y=matrix(runif(p*p),ncol=p)-0.5+tilt
   xy=x%*%y
   covar=t(xy)%*%xy
   covar=covar+diag(shrink*diag(covar))
   covar=abs(covar)
   corr=diag(1/sqrt(diag(covar)))%*%covar%*%diag(1/sqrt(diag(covar)))
            # generate correlation matrix
   if(1==2){
     for(i in 1:(p-1)){ #select correlation to keep
       for(j in (i+1):p){
         corr[i,j]=0.5
         if(j>i+1){
           corr[i,j]=0
          }
         corr[j,i]=corr[i,j]
        }
       }
   }
   svdcorr=svd(corr) # singular value decomposition
   sqrtsig=svdcorr$u%*%diag(sqrt(svdcorr$d))%*%t(svdcorr$v)# square-root a matrix

   list(sqrtsig)
}

trueR2=function(nrep,p,beta,xsig,errsig,powx,powy,sd){
  # generate data from linear model to calculate R2
  # powx, powy parameters making non-normal x and random error
  # sd, square-root of correlation matrix, making correlated x.

  xy=datgen(nrep,p,beta,xsig,errsig,powx,powy,sd)
  s2=var(xy[[1]]%*%beta)
  r2=s2/var(xy[[2]])

  list(c(r2,s2))
}

