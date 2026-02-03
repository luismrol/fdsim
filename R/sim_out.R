#' Simulation algorithm for multivariate functional data for analyzing heterogeneous outliers
#'
#' @description
#' This function uses funData package to generate univariate functional observations
#' with different types of outliers and different degrees of contamination
#'
#' @param nrep    Number of iterations to be performed
#' @param N       Number of observations at each iteration
#' @param p       Number of discretizations for t
#' @param moderate   (LOGICAL) Is it the moderate or the extreme scenario?
#' @param perc_out a vector of two entries  c(amplitude, shape) for the proportion of
#' observations outlying in each case.
#' @return A list of dimension nrep of simulated functional
#' @importFrom funData simFunData
#' @export


# This generates the man/*.Rd help files and updates NAMESPACE

sim_out<-function(nrep, N, p, moderate = TRUE, perc_out = c(0,0)){
  if(length(perc_out)!= 2){
    stop("Percentage of outliers should be a numeric vector of dimension 2")
  }
  if(!is.numeric(perc_out[1]) | !is.numeric(perc_out[2])){
    stop("both elements in perc_out should be numeric")
  }
  argvals<-seq(0,1, by = 1/(p-1))
  mu_t<-30*argvals*(1-argvals)^(2)
  if (moderate == TRUE){
    amplitude<-2*mu_t
    shape<-5-10*((argvals)^(3/2))*(1-(argvals))^(3/2)
  }else{
    amplitude<-3*mu_t
    shape<--20*(argvals^(3/2))*((1-argvals)^(3/2))
  }
  perc_amp<-perc_out[1]
  perc_shape<-perc_out[2]
  class<-array()
  ldata<-list()
  for (i in 1:nrep){
    fdata1<-funData::simFunData(argvals, M=10, "Fourier", eValType = "exponential", N=N)
    sdata<-fdata1$simData@X
    for (j in 1:nrow(sdata)){
      if (j<=(1-perc_amp-perc_shape)*nrow(sdata)){
        sdata[j,]<-sdata[j,]+mu_t
        class[j]<-"1.Base"
      }  else if (j<=(1-perc_shape)*nrow(sdata)){
        sdata[j,]<-sdata[j,]+amplitude
        class[j]<-"2.Amplitude"
      }
      else{
        sdata[j,]<-sdata[j,]+shape
        class[j]<-"3.Shape"
      }    }
    ldata[[i]] <- cbind(sdata, class)
  }
  return(ldata)
}
