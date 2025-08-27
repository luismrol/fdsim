#' Simulation algorithm for multivariate functional data from Claeskens et al 2014
#'
#' @description
#' This function simulates multivariate functional data based on a mean function and a covariance matrix
#' @param N       Number of functional observations
#' @param P       Grid for the final data
#' @param grid1,grid2    Numeric (Seq) grids of discrete points in which the function is mapped (1) and the heteroskedasticity component (2)
#' @param der     Boolean indicating if numerical derivatives should be calculated
#' @param covar   Covariance method. Could be absolute exponential decay ("abs") or cuadratic exponential decay ("sq")
#' @param delta   Delta parameter of the covariance kernel
#' @param het     Heteroskedasticity. Default is TRUE. If False, D is an identity matrix
#' @param method  Simulation method
#' @param lognormal if TRUE, the distribution is assumed to be lognormal. If FALSE, assumed to be normal.
#' @return        A list of simulated functional data
#' @export


# This generates the man/*.Rd help files and updates NAMESPACE

mfd_sim_h <- function(N, P, grid, mu, der = FALSE, covar = "sq",
                    method = c("svd", "chol", "eigen"), het = TRUE, lognormal = TRUE, delta = 0.25){


  # Depending on which type of decay we are inputting for the simulation,
  # we define the variance kernel.
  # the parameter delta is associated to the decay rate. A higher delta means a lower
  # correlation.
  if (!is.numeric(delta) || delta <= 0) {
    stop("`delta` must be a positive number.")
  }
  if (!(method %in% c("svd", "chol", "eigen"))) {
    stop("`Method should be either svd, chol or eigen")
  }

  if (covar == "sq"){
    kernel <- function(a, b, delta) {
      as.matrix(outer(a, b, function(i, j) exp(-((i - j)^2) / ((2 * delta)^2))))
    }
  }
  else if (covar == "abs"){
    kernel <- function(a, b, delta) {
      as.matrix(outer(a, b, function(i, j) exp(-(abs(i - j)) / ((2 * delta)^2))))
    }
  }
  der <- der
  ## Generate grid1 for the whole domain (thicker grid), grid2 for a less dense
  ## grid for the first estimation.
  grid2<-grid
  grid1<-seq(min(grid), max(grid), by = max(grid)/(P-1))

  # Partial kernels
  Kxx <- kernel(grid2, grid2, delta)
  Ktx <- kernel(grid1, grid2, delta)
  Ktt <- kernel(grid1, grid1, delta)

  # For each element, give the minimum between 1 and median(grid1)-grid2^2. This is constraining
  # the effect of Kxx to a subset given by min (pi-grid2^2). The original was with Pi instead of median(grid1), but this
  # function makes it more general to any domain of grid1
  if (isTRUE(het)){
    D <- diag(pmin(((0.5*median(grid1)) - grid2)^2, 1))
  }
  else {
    D <- diag(1, nrow = length(grid2))
  }

  K_m <- Ktx %*% solve(Kxx + D)
  mu_x <- mu
  mu_t <- as.vector(K_m %*% t(mu_x))
  Sigma <- ( Ktt - Ktx %*% solve(Kxx + D) %*% t(Ktx))
  Y <- list()
  # Generate log-normal process
  if (isTRUE(lognormal)){
    Y[[1]] <- as.matrix(exp(mvtnorm::rmvnorm(n = N, mean = mu_t, sigma = Sigma, method = method)))
  }
  else{
    Y[[1]] <- as.matrix(mvtnorm::rmvnorm(n = N, mean = mu_t, sigma = Sigma, method = method))
  }
  den<-as.numeric(grid1[2]-grid1[1])
  if (isTRUE(der)){
    dY<-matrix(NA, ncol= ncol(Y[[1]])-1, nrow = nrow(Y[[1]]))
    for (i in 2:ncol(Y[[1]])){
      dY[,i-1]<-(Y[[1]][,i]-Y[[1]][,(i-1)])/den
    }
    Y[[2]]<- as.matrix(dY)
    }
  return(Y)
}
