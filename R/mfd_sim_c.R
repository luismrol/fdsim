#' Simulation algorithm for multivariate functional data from Claeskens et al 2014
#'
#' @description
#' This function simulates multivariate functional data using either a mixture of function without
#' covariance specification, or a heteroskedastic function from a mean and random noise
#'
#' @param N       Number of functional observations
#' @param grid,grid2    Grids of discrete points in which the function is mapped (1) and the heteroskedasticity component (2)
#' @param der     Boolean indicating if numerical derivatives should be calculating
#' @param covar   Covariance method. Could be absolute exponential decay ("abs") or cuadratic exponential decay ("sq")
#' @param delta   Delta parameter of the covariance kernel
#' @param method  Simulation method
#' @return        A list of simulated functional data
#' @export


# This generates the man/*.Rd help files and updates NAMESPACE

mfd_sim_c <- function(N, grid1, grid2, der = FALSE, covar = "sq",
                    method = c("svd", "chol", "eigen")){
  # Columns for each dimension
  if (covar == "sq"){
    kernel <- function(a, b, delta) {
      outer(a, b, function(i, j) exp(-((i - j)^2) / ((2 * delta)^2)))
    }
  }
  else if (covar == "abs"){
    kernel <- function(a, b, delta) {
      outer(a, b, function(i, j) exp(-(abs(i - j)) / ((2 * delta)^2)))
    }
  }
  P <- ncol(mu)
  N <- floor(N)
  # Number of dimensions
  L <- nrow(mu)
  Kxx <- kernel(grid2, grid2, delta)
  Ktx <- kernel(grid1, grid2, delta)
  Ktt <- kernel(grid1, grid1, delta)
  D <- diag(pmin(((mean(grid1)) - x)^2, 1))
  K_m <- Ktx %*% solve(Kxx + D)
  mu_x <- as.matrix(runif(1, -2, 2)) %*% t(sin(x)) + as.matrix(runif(1, -1, 1)) %*% t(cos(x))
  mu_t <- as.vector(K_m %*% t(mu_x))
  Sigma <- K_tt - K_tx %*% solve(K_xx + D) %*% t(K_tx)
  Y <- list()
  Y[[1]] <- as.matrix(exp(mvtnorm::rmvnorm(n = N, mean = mu_t, sigma = Sigma, method = method)))
  if (der == TRUE){
    dY<-matrix(NA, ncol= ncol(Y[[1]])-1, nrow = nrow(Y[[1]]))
    for (i in 2:ncol(Y[[1]])){
      dY[,i-1]<-(Y[[1]][,i]-Y[[1]][,(i-1)])/0.06346652
    }
    Y[[2]]<- dY
    }
  return(Y)
}
