#-----------------------------------------------------------------------
# Project exposure to the orthogonal complement of confounders  space
#-----------------------------------------------------------------------
#' @import stats
NULL
#' Project x to the orthogonal complement of c.
#'
#' This function project exposures x to the orthogonal complement of confounders c.
#'
#' @param x a matrix of nxp dimensional exposures
#' @param cfd a matrix of nxq dimensional confounders
#'
#' @details This function uses c to form projection P=I-c(c'c)^(-)c'
#' and then compute projection z=Px
#'
#' @return projtocorth, the projected x.
#'
#' @export
#'
project=function(cfd,x){

  n=dim(cfd)[1]
  q=dim(cfd)[2]
  for(j in 1:q){
    muc=mean(cfd[,j])
    sdc=sd(cfd[,j])

    cfd[,j]=(cfd[,j]-muc)/sdc
  }

  csvd=svd(cfd,nv=0)
  projxtocorth=x-csvd$u%*%diag(1/csvd$d^2)%*%t(csvd$u)%*%x
    #project x to the orthogonal linear space of cfd column vectors.

  list(projxtocorth)
}

#----------------------------------------#
#    decorrelation transformation        #
#----------------------------------------#
#' Decorrelation transformation
#'
#' This function performs a decorrelation transformation of a data set.
#'
#' @param x a matrix of nxp dimension
#'
#' @details This function use the input data to estimate the covariance matrix,
#' and then use the estimated covariate matrix to decorrelate
#' the input data matrix. require n>p.
#'
#' @return Decorrelated data matrix.
#'
#' @export
transf=function(x){

  xsvd=svd(x)
  tran=xsvd$v%*%diag(1/sqrt(xsvd$d))%*%t(xsvd$v)
  z=as.matrix(x)%*%tran

  list(z)
}

#' Decorrelation
#'
#' @param x a matrix of nxp dimension. Input matrix to be transformed.
#' @param X a matrix of Nxp dimension. Supplementary data for covariance matrix estimation.
#'
#' @details This function use both the input data and the supplementary data to estimate
#' the covariance matrix, and then use the estimated covariate matrix to decorrelate
#' the input data matrix. require (n+N)>p.
#'
#' @return Decorrelated data matrix.
#'
#' @export
transfsd=function(x,X){

  Xsvd=svd(rbind(x,X))
  tran=Xsvd$v%*%diag(1/sqrt(Xsvd$d))%*%t(Xsvd$v)
  z=as.matrix(x)%*%tran

  list(z)
}


#' Rescale data matrix
#'
#' Rescale data matrix by subtracting mean and divide the standard deviation for each variable.
#'
#' @param z input matrix of nxp dimension
#'
#' @return rescaled data matrix of nxp dimension
#'
#' @export
zscale=function(z){
  n=dim(z)[1]
  p=dim(z)[2]
  for(j in 1:p){
    mu=mean(c(z[,j]))
    sd=sd(c(z[,j]))
    z[,j]=(z[,j]-mu)/sd
  }
  list(z)
}

#' Output simulation results.
#'
#' This function summarize the simulation results with the input of simulation results.
#'
#' @param result a vector of simulation results
#' @param pa length of alpha (type I error vector)
#' @param truth  true (proportion of) the explained variation
#' @param rept the number of simulations
#'
#' @return  empirical confidence interval length and coverage.
#'
#' @export
#'
Soutput=function(result,pa,truth,rept){
  print(apply(result,2,mean))
  print(apply(result,2,var))

  for(k in 1:(2*pa)){
    print(mean(result[,3+2*k]-result[,2+2*k]))
    print(mean(100*(result[,2+2*k]<=rep(truth,rept))*(result[,3+2*k]>=rep(truth,rept))))
  }
}

