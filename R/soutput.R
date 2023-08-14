#' Output simulation results.
#'
#' This function summarize the simulation results with the input of simulation results.
#'
#' @param result a vector of simulation results
#' @param r20  true proportion of the explained variation
#' @param rept the number of simulations
#'
#' @return  empirical confidence interval length and coverage.
#'
#' @export
#'
Soutput=function(result,r20,rept){
  print(apply(result,2,mean))
  print(apply(result,2,var))
  print(mean(result[,5]-result[,4]))
  print(mean(result[,7]-result[,6]))
  print(mean(result[,9]-result[,8]))

  print(mean(100*(result[,4]<=rep(r20,rept))*(result[,5]>=rep(r20,rept))))
  print(mean(100*(result[,6]<=rep(r20,rept))*(result[,7]>=rep(r20,rept))))
  print(mean(100*(result[,8]<=rep(r20,rept))*(result[,9]>=rep(r20,rept))))

  print(mean(result[,11]-result[,10]))
  print(mean(result[,13]-result[,12]))
  print(mean(result[,15]-result[,14]))

  print(mean(100*(result[,10]<=rep(r20,rept))*(result[,11]>=rep(r20,rept))))
  print(mean(100*(result[,12]<=rep(r20,rept))*(result[,13]>=rep(r20,rept))))
  print(mean(100*(result[,14]<=rep(r20,rept))*(result[,15]>=rep(r20,rept))))
}
