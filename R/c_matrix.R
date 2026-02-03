#' Confusion matrix for the outlier detection algorithm in simulation
#'
#' @description
#' This is an intermediate function to present the results of the simulation exercise
#' for the heterogeneous outlier detection in functional data.
#'
#'
#' @param data   The data input. If you are inputting a data given from sim_out, then remember to drop the last column
#' @param class  A string vector containing the types of outliers.
#' @param moderate   (LOGICAL) Is it the moderate or the extreme scenario?
#' @param perc_out a vector of two entries  c(amplitude, shape) for the proportion of
#' observations outlying in each case.
#' @return A named matrix which columns are, respectively, the true positive,
#' false positive, false negative and true negative counts for each iteration of the simulation.


#' @export

c_matrix<- function(data, method){
  output<-do.call(rbind, lapply(ldata, function(data, class=NULL, method ){
    if (is.null(class)){
      class<-apply(data[,ncol(data)],2, as.numeric)
      data<-data[,(1:ncol(data)-1)]
    }
    tp<-NA
    fp<-NA
    tn<-NA
    fn<-NA
    N = nrow(data)
    outliers<-method(data)$outliers
    tp<-sum(as.numeric(outliers %in% which(class!="1.Base")))
    fp<-sum(as.numeric(outliers %in% which(class=="1.Base")))
    tn<-sum(as.numeric(which(c(1:N) %in% outliers == FALSE) %in%  which(class=="1.Base")))
    fn<-sum(as.numeric(which(c(1:N) %in% outliers == FALSE) %in% which(class!="1.Base")))

    return(c(tp, fp, fn, tn))
  }))
  colnames(output) = c("tp", "fp", "fn", "tn")
  return(output)}
