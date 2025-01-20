#'  Select data with index in dat1, but not in dat2. Assuming dat2 is part of dat1.
#'
#' @param dat1 The first data set to select data from
#' @param dat2 The second data set to exclude data from data1 whose indices are in dat2.
#' @param vind1 The index variable in dat1 used to determine the selection
#' @param vind2 The index variable in dat2 used to exclude the data from dat1 whose indices match.
#'        As such, vind1 in dat1 and vind2 in dat2 should be the same variable.
#'
#' @details Separate dat1 into two data sets with one having indices of vind1 in dat1, but not in dat2,
#'          and the other having indices match those of vind2 in dat2 ()
#'
#' @return Two data sets
#' Note that we can use merge function to select data with common index data.
#'
#' @examples \dontrun{seldifdat(dat1, dat2, vind1=1,vind2=1)}
#'
#' @export
#'
seldifdat=function(dat1, dat2, vind1,vind2){

  m=dim(dat1)[1]
  n=dim(dat2)[1]
  N=m-n
  select=rep(0,N)
  xf=as.numeric(dat2[1,vind2]) # initialize
  j=1
  k=0
  for(i in 1:m){
    if(as.numeric(dat1[i,vind1])==xf){
      if(j!=n){
        j=j+1
        xf=as.numeric(dat2[j,vind2])
      }
    }else{
      k=k+1
      select[k]=i
    }
  }
  x=dat1[select,]
  x1=dat1[setdiff(c(1:m),select),]

  list(x,x1)
}

#' project data matrix x to the orthogonal complement of data matrix z.
#'
#' This function uses linear projector (I-z(z^tz)^(-1)z^t) to project.
#'
#' @param x matrix (nxp) to be projected.
#' @param z matrix (nxq) the projection is based on.
#'
#' @details compute (I-z(z^tz)^(-1)z^t)x
#'
#' @return the projected x.
#'
#' @examples \dontrun{projed(x,z)}
#'
#' @export
#'

projed=function(x,z){ # with supplementary data
  #projection of x on the orthogonal complement of z
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(z)[2]

  for(j in 1:p){
    x[,j]=(x[,j]-mean(c(x[,j])))/sd(c(x[,j]))
  }
  for(j in 1:q){
    z[,j]=(z[,j]-mean(c(z[,j])))/sd(c(z[,j]))
  }

  projz=diag(rep(1,n))-z%*%chol2inv(chol(t(z)%*%z))%*%t(z)
  pjx=projz%*%x

  list(pjx)
}

#' Perform data analysis to estimate the proportion of explained variation from
#' outcome data: dataseta1; covariate data: dataseta2;
#' and the suplementary covariate data: dataseta3.
#'
#' This function performance data analysis with multiply imputed datasets
#'
#' @param dataseta1 multiply imputed data set for outcomes and confounders.
#' @param dataseta2 multiply imputed data set for covariates.
#' @param dataseta3 multiply imputed data set for supplementary covariates.
#' @param nimpute copies of imputed data sets.
#' @param alpha to determine confidence as 1-alpha.
#'
#' @details perform both unadjusted and adjusted analyses and using Rubin's rule to combine results
#'
#' @return estimators, variance estimates, and confidence intervals.
#'
#' @examples \dontrun{jasaDataAnalysis(outcome,covariates,supplecov)}
#'
#' @export

jasaDataAnalysis=function(dataseta1,dataseta2,dataseta3,nimpute=6,alpha=0.05){
## 1. read data
 #file1="F:/seagatebackupDec3-2024/NIEHSpapers/Combined/JASA/DataAnalysis/BPs"
 #file2="F:/seagatebackupDec3-2024/NIEHSpapers/Combined/JASA/DataAnalysis/CCPCBs"
 #file3="F:/seagatebackupDec3-2024/NIEHSpapers/Combined/JASA/DataAnalysis/IMPallcov"

 #a1=read.table(file1,header =F,sep = " ", quote = "\"", dec = ".", fill = TRUE)
 #a2=read.table(file2,header =F,sep = " ", quote = "\"", dec = ".", fill = TRUE)
 #a3=read.table(file3,header =F,sep = " ", quote = "\"", dec = ".", fill = TRUE)

## 2. Full data analysis
  fulldata=merge(dataseta1,dataseta2, by="V1",)
  fulldata=as.matrix(fulldata)
  y=fulldata[,2]-fulldata[,3] # blood pressure difference
  x=fulldata[,4:65]

  # nimpute=6
  m=dim(dataseta3)[1]/nimpute
  palpha=alpha
  ilam=0.3
  aa=array(0,c(4,12,nimpute))
  bb=aa
  for(k in 1:nimpute){ # Five imputed data sets
    print(c(k,k,k))
    X=seldifdat(as.matrix(dataseta3[((k-1)*m+1):(k*m),]),fulldata,1,1)
    XX=X[[1]];xx=X[[2]]
    CF=XX[,3:9];xsup=XX[,10:71]  #supplementary data part
    cf=xx[,3:9];x=xx[,10:71]  #full data part, x1 should be x in the full data


  #1. Unadjusted analysis

    aa1=TEV::RVep(y,x,alpha=palpha)
    aa[,1,k]=c(aa1[[1]],0,aa1[[2]])

    aa2=TEV::RVee(y,x,lam=ilam, niter=0, alpha=palpha)
    aa[,2,k]=c(aa2[[1]][1:2],aa2[[2]])

    aa3=TEV::RVee(y,x,lam=ilam,niter=20,alpha=palpha)
    aa[,3,k]=c(aa3[[1]][1:2],aa3[[2]])

    aa4=TEV::RVls(y,x,alpha=palpha)
    aa[,4,k]=c(aa4[[1]][1:2],aa4[[2]])

    aa5=TEV::GRE(y,x)
    aa[,5,k]=c(aa5[[1]],aa5[[2]],aa5[[3]])

    aa6=TEV::RVsd(y,x,xsup, lam =ilam, alpha=palpha, niter =0, know="no")
    aa[,6,k]=c(aa6[[1]][1:2],aa6[[2]])

    aa7=TEV::RVsd(y,x,xsup, lam =ilam, alpha=palpha, niter =20, know="no")
    aa[,7,k]=c(aa7[[1]][1:2],aa7[[2]])

    aa8=TEV::CHIVE(y,x,xsup,alpha=palpha)
    aa[,8,k]=c(aa8[[1]][1:2],aa8[[2]])

    aa9=TEV::CHIVE(y,x,alpha=palpha)
    aa[,9,k]=c(aa9[[1]][1:2],aa9[[2]])

    aa10=TEV::RVmlde(y,x,alpha=palpha)    #MLDE
    aa[,10,k]=c(aa10[[1]][1:2],aa10[[2]])

    aa11=TEV::FULLREML(y,x,alpha=palpha) #MLRE
    aa[,11,k]=c(aa11[[1]][1:2],aa11[[2]])

    aa12=TEV::GCTAREML(y,x,alpha=palpha) #EEML
    aa[,12,k]=c(aa12[[1]],aa12[[2]])


  #2. Adjusted analysis

    xsup1=projed(xsup,CF)[[1]] # confounder adjusted for supplementary data
    x1=projed(x,cf)[[1]]       # confounder adjusted for original data
    fit=lm(y~cf)
    yr=as.numeric(fit$resid)   #adjusted outcome

    bb1=TEV::RVep(yr,x1,alpha=palpha)
    bb[,1,k]=c(bb1[[1]],0,bb1[[2]])

    bb2=TEV::RVee(yr,x1,lam=ilam, niter=0, alpha=palpha)
    bb[,2,k]=c(bb2[[1]][1:2],bb2[[2]])

    bb3=TEV::RVee(yr,x1,lam=ilam,niter=20,alpha=palpha)
    bb[,3,k]=c(bb3[[1]][1:2],bb3[[2]])

    bb4=TEV::RVls(yr,x1,alpha=palpha)
    bb[,4,k]=c(bb4[[1]][1:2],bb4[[2]])

    bb5=TEV::GRE(yr,x1)
    bb[,5,k]=c(bb5[[1]],bb5[[2]],bb5[[3]])

    bb6=TEV::RVsd(yr,x1,xsup1, lam =ilam, alpha=palpha, niter =0, know="no")
    bb[,6,k]=c(bb6[[1]][1:2],bb6[[2]])

    bb7=TEV::RVsd(yr,x1,xsup1, lam =ilam, alpha=palpha, niter =20, know="no")
    bb[,7,k]=c(bb7[[1]][1:2],bb7[[2]])

    bb8=TEV::CHIVE(yr,x1,alpha=palpha)
    bb[,8,k]=c(bb8[[1]][1:2],bb8[[2]])

    bb9=TEV::CHIVE(yr,x1,xsup1,alpha=palpha)
    bb[,9,k]=c(bb9[[1]][1:2],bb9[[2]])

    bb10=TEV::RVmlde(yr,x1,alpha=palpha)    #MLDE
    bb[,10,k]=c(bb10[[1]],bb10[[2]])

    bb11=TEV::FULLREML(yr,x1,alpha=palpha) #MLRE
    bb[,11,k]=c(bb11[[1]],bb11[[2]])

    bb12=TEV::GCTAREML(yr,x1,alpha=palpha) #EEML
    bb[,12,k]=c(bb12[[1]],bb12[[2]])
  }

  maa=apply(aa,c(1,2),mean)
  vaa=apply(aa,c(1,2),var)
  print("Unadjusted analysis")
  print("Methods: 1.EP, 2.EE, 3.EEit, 4.LS, 5.GRE, 6.SD, 7.SDit, 8.CHIVE0, 9.CHIVE, 10.MLDE, 11.MLRE, 12.EEML")
  print("Rubin's combination for estimator and variance")
  print( cbind(maa[1,],maa[2,]+vaa[1,]) )
  print("CI(general) and CI(normal)")
  print( cbind(maa[1,]-1.96*sqrt(maa[2,]+vaa[1,]),maa[1,]+1.96*sqrt(maa[2,]+vaa[1,]),
        maa[3,],maa[4,]) )

  print("Adjusted analysis")
  mbb=apply(bb,c(1,2),mean)
  vbb=apply(bb,c(1,2),var)
  print( cbind(mbb[1,],mbb[2,]+vbb[1,]) )
  print( cbind(mbb[1,]-1.96*sqrt(mbb[2,]+vbb[1,]),mbb[1,]+1.96*sqrt(mbb[2,]+vbb[1,]),
        mbb[3,],mbb[4,]) )
}
