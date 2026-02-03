#' Summary of simulation results
#'
#' @description
#' Basic summary of the confusion matrix.
#'
#'
#' @param results   matrix given from the c_matrix function. Its columns are the counting of true and
#' false positive, false and true negatives.
#'
#' @return a vector of mean and standard deviation ,
#' false positive, false negative and true negative counts for each iteration of the simulation.


#' @export

sum_results<-function(results){
  avg<-colMeans(results)
  sd<-sqrt(diag(var(results)))
  return(data.frame("mean"= avg, "sd"= sd))
}
