#' Basic functions for univariate depth calculation
#'
#' @description
#' This function calculates column-wise univariate depths that would serve as inputs for integrate and geometric
#' multivariate and functional depths.
#'
#' @param data    An NxP Matrix, from which a column-wise depth is calculated
#' @param depth   What depth is going to be used. Supported are Mahalanobis, Projection, Univariate Tukey and Simplicial
#' @param out     Should the associated outlyingness be also included in the result?
#' @export


# This generates the man/*.Rd help files and updates NAMESPACE

uni_depth <- function(data, depth = "maha"){
  depth = match.arg(depth, c("mahalanobis", "simplicial", "halfspace", "tukey", "projection"))
  ### Mahalanobis depth
  mahadepth<-function(vec){
      m<-mean(vec)
      s2<-var(vec)
      return (1/(1+(vec-m)^2/s2))
  }
  ## Simplicial
  usdepth<-function(vec){
    possible<-0
    actual<-0
    possiblevec<-NA
    actualvec<-NA
    ratevec<-NA
    for (x in 1:length(vec)){
      for(i in 1:length(vec)){
        for (j in 1:length(vec)){
          possible<-possible+1
          if((vec[x]>min(vec[i],vec[j])&(vec[x]<max(vec[i],vec[j])))){
            actual <- actual + 1
          }
        }
      }
      possiblevec[x]<-possible
      actualvec[x]<-actual
      ratevec[x]<-actual/possible
      possible <- 0
      actual <- 0
    }
    d<-ratevec/max(ratevec)
    return(d)
  }
  ### Projection depth
  pd<-function(vec){
    vec<-as.numeric(as.matrix(vec))
    median<-median(vec)
    mad<-mad(vec)
    sdout<-abs((vec-median))/mad
    pd<-1/(1+sdout)
    return(pd)
  }
  ### Tukey depth
  udepth<-function(vec){
    vec1<-rank(vec)
    N <-length(vec)
    vec2<-NA
    for (i in 1:N){
      vec2[i]<-min(vec1[i]/N, 1-(vec1[i]/N))}
    return(vec2)
  }
  if (depth == "mahalanobis"){
    d = apply(data, 2, mahadepth)
  } else if (depth == "simplicial"){
    d = apply(data, 2, usdepth)
  } else if (depth == "projection"){
    d = apply(data, 2, pd)
  } else if (depth == "tukey"){
    d = apply(data, 2, udepth)
  }
  out = ((d + (10^-10))^-1)-1
  return(list("depth" = d, "outlyingness" = out, "mean depth" = rowMeans(d), "mean outlyingness" = rowMeans(out)))
}
