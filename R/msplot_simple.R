#' Simplification of msplot
#'
#' @description
#' Function to simplify the msplot
#'
#' @param data   The data input. If you are inputting a data given from sim_out, then remember to drop the last column
#'
#' @return data and outliers
#' false positive, false negative and true negative counts for each iteration of the simulation.

#' @importFrom fdaoutlier msplot

#' @export

msplot_simple<-function(data){
  outliers<-fdaoutlier::msplot(data, plot = FALSE)$outliers
  return(list("data"=data, "outliers"=outliers))
}


