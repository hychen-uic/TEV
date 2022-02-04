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
