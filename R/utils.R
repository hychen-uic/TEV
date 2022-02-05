#----------------------------------------#
#    decorrelation transformation        #
#----------------------------------------#

#' Decorrelation transformation
#'
#' @param x input
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
#' @param x input 1
#' @param X input 2
#'
#' @export
transfsd=function(x,X){

  Xsvd=svd(rbind(x,X))
  tran=Xsvd$v%*%diag(1/sqrt(Xsvd$d))%*%t(Xsvd$v)
  z=as.matrix(x)%*%tran

  list(z)
}

### projection when n>p ###
proj1=function(x){ # without supplementary data

  n=dim(x)[1]
  p=dim(x)[2]
  for(j in 1:p){
    mux=mean(c(x[,j]))
    sdx=sd(c(x[,j]))

    x[,j]=(x[,j]-mux)/sdx
  }

  proj=diag(rep(1,n))-x%*%solve(t(x)%*%x)%*%t(x)

  list(proj)
}

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

  projx=diag(rep(1,n))-xx[1:n,]%*%solve(t(xx)%*%xx)%*%t(xx[1:n,])*(n+N)/n
  projX=diag(rep(1:N))-xx[(n+1):(n+N),]%*%solve(t(xx)%*%xx)%*%t(xx[(n+1):(n+N),])*(n+N)/N

  list(projx,projX)
}

#' Transformation function
#'
#' @param z input argument
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
