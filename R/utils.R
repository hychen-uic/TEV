#----------------------------------------#
#    decorrelation transformation        #
#----------------------------------------#

#' Decorrelation transformation
#'
#' This function performs a decorrelation transformation of a data set.
#'
#' @param x a matrix of nxp dimension.
#'
#' @details This function use the input data to estimate the covariance matrix,
#' and then use the estimated covariate matrix to decorrelate
#' the input data matrix. require n>p.
#'
#' @return decorrelation data matrix.
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
#' Decorrelation with the aid of supplementary data.
#'
#' @param x a matrix of nxp dimension. Input matrix to be transformed.
#' @param X a matrix of Nxp dimension. Supplementary data for covariance matrix estimation.
#'
#' @details This function use both the input data and the supplementary data to estimate
#' the covariance matrix, and then use the estimated covariate matrix to decorrelate
#' the input data matrix. require (n+N)>p.
#'
#' @return decorrelation data matrix.
#'
#' @export
transfsd=function(x,X){

  Xsvd=svd(rbind(x,X))
  tran=Xsvd$v%*%diag(1/sqrt(Xsvd$d))%*%t(Xsvd$v)
  z=as.matrix(x)%*%tran

  list(z)
}

### projection when n>p ###
#' Compute projection matrix
#'
#' @param x input matrix of nxp dimension.
#'
#' @return projection matrix of nxn.
#'
#' @export
proj1=function(x){ # without supplementary data

  n=dim(x)[1]
  p=dim(x)[2]
  for(j in 1:p){
    mux=mean(c(x[,j]))
    sdx=sd(c(x[,j]))

    x[,j]=(x[,j]-mux)/sdx
  }

  proj=diag(rep(1,n))-x%*%chol2inv(chol(t(x)%*%x))%*%t(x)

  list(proj)
}

#' Compute projection matrices
#'
#' @param x input matrix of nxp dimension.
#' @param X input matrix of Nxp dimension.
#'
#' @return projection matrix of nxn and projection matrix of NxN based on (x, X).
#'
#' @export
proj2=function(x,X){ # with supplementary data

  n=dim(x)[1]
  p=dim(x)[2]
  N=dim(X)[1]
  xx=rbind(x,X)
  for(j in 1:p){
    muxx=mean(c(xx[,j]))
    sdxx=sd(c(xx[,j]))

    xx[,j]=(xx[,j]-muxx)/sdxx
  }

  projx=diag(rep(1,n))-xx[1:n,]%*%chol2inv(chol(t(xx)%*%xx))%*%t(xx[1:n,])*(n+N)/n
  projX=diag(rep(1:N))-xx[(n+1):(n+N),]%*%chol2inv(chol(t(xx)%*%xx))%*%t(xx[(n+1):(n+N),])*(n+N)/N

  list(projx,projX)
}

#' Compute projection of x on the orthogonal complement of z
#'
#' @param x input matrix of nxp dimension.
#' @param z input matrix of nxq dimension.
#'
#' @return a nxp matrix, the projection of x on the orthogonal complement of z.
#'
#' @export
projed=function(x,z){ # with supplementary data

  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(z)[2]

  for(j in 1:p){
    xx[,j]=(x[,j]-mean(c(x[,j])))/sd(c(x[,j]))
  }
  for(j in 1:q){
    zz[,j]=(z[,j]-mean(c(z[,j])))/sd(c(z[,j]))
  }

  projz=diag(rep(1,n))-zz%*%chol2inv(chol((t(zz)%*%zz)))%*%t(zz)

  projedx=projz%*%xx

  list(projedx)
}

#' Rescale data matrix
#'
#' Rescale data matrix by subtracting mean and divide the standard deviation for each variable.
#'
#' @param z input matrix of nxp dimension.
#'
#' @return rescaled data matrix of nxp dimension.
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
