#' Simulation algorithm for univariate functional data with heterogeneous outliers
#'
#' @description
#' This function uses `funData` package to generate univariate functional observations
#' with different types of outliers and different degrees of contamination.
#'
#' @param N     Number of functional observations.
#' @param Nt    Number of discrete time points.
#' @param mean  A list of mean functions for base data and outlier types.
#' @param prop  A vector of proportions for base data and outlier types.
#' @param M     (Numeric) the number of eignevalues / eigenfunctions to be considered in simFunData
#' @param eFunType (Character) the type of eigenfunction to be simulated for the random term. Types include "Poly", "PolyHigh", "Fourier", "FourierLin", and "Wiener"
#' @param eValType (Character) The type of decay for the eigenvalues in the random component. Could be "linear", "exponential" or "wiener"
#' @return      An Nx(P+1) matrix containing functional data and a class membership indicator for base data and outlier types.
#' @importFrom funData simFunData
#' @export
#'
simFunDataOut <- function(N, Nt, mean, prop, M=10, eFunType = "Fourier", eValType = "exponential"){
  ## Generate the grid T
  argvals <- seq(0,1, by = 1/(Nt-1))
  ## translate the mu to the coordinates of T
  f<-mean
  mus_t <- lapply(f, function(f) f(argvals))
  ## Put it in matrix
  aux <- prop*N
  mus_m <- do.call(
    rbind,
    Map(function(x, k) {
      matrix(rep(x, each = k), ncol = length(x), byrow = TRUE)
    }, mus_t, aux)
  )
  ## Generate epsilon
  eps_m <- funData::simFunData(argvals, M=M, eFunType = eFunType, eValType = eValType, N=N)$simData@X
  X <- as.matrix(mus_m + eps_m)
  cl <- rep(seq_along(aux), times = aux)
  X <- cbind(X, cl)
  return(X)
}
