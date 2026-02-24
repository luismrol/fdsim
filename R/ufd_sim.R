#' Simulation algorithm for multivariate functional data based on Ieva and Paganoni (2017)
#'
#' @description
#' This function simulates multivariate functional data using Eigen, Cholesky and Singular
#' Value Decomposition, based on the method proposed by Ieva and Paganoni (2017) and implemented
#' initially in the library 'roahd', but with some changes enabling data exploration.
#'
#' @param N       Number of functional observations
#' @param mu      Matrix with 1 row and P columns containing the mean curve
#' @param grid    If you want to specify your own grid, the vector of the grid. Length must coincide with the number of elements. Default is FALSE, in which case a [0,1] grid with even spacing
#'                over P points is implemented.
#' @param covar   Either a covariance matrix or a type of exponential decay, either "sq" for square or "abs" for absolute exponential decay
#' @param method  Factorization method for the covariance matrix. Either "eigen" for eigendecomposition or "chol" for Choleski factorization.
#' @return a list with "x", the simulated data, "covar", the covariance matrix that was used as input, "factorized covar" the covariance matrix after factorization, "Z_matrix", each of the elements of factorization.
#' @export


# This generates the man/*.Rd help files and updates NAMESPACE

ufd_sim <- function(N, mu, grid = FALSE, covar = "sq",
                    method = c("chol", "eigen"), delta = 0.5){
  if (!is.numeric(mu) && !is.matrix(mu)){
    stop(paste("Mu should be Px1 matrix or a numeric array. An object of class", class(mu), "is not valid"))
  }
  if (!is.numeric(N)){
    stop(paste("Sample size must be a positive integer. Class", class(N)[1], "is not valid"))
  }
  P <- length(mu)
  rand <- matrix(rnorm(N*P), N, P)
  if (isFALSE(grid)){
    grid = seq(0,1, by = 1/(P-1))
  } else if (length(grid)!=ncol(mu)){
    stop("size of grid and size of mean vector are not equal")
  }
  # Columns for each dimension
  if (is.matrix(covar)){
    covar = covar
  } else {
    covar1<-covar
    if (covar1 == "sq"){
      kernel <- function(a, b, delta) {
        as.matrix(outer(a, b, function(i, j) exp(-((i - j)^2) / ((delta)^2))))
      }
    }
    else if (covar1 == "abs"){
      kernel <- function(a, b, delta) {
        as.matrix(outer(a, b, function(i, j) exp(-(abs(i - j)) / ((delta)))))
      }
    }
    covar = kernel(grid, grid, delta)
  }
  c<-covar
  if (method == "chol"){
      cv <- chol(as.matrix(Matrix::nearPD(covar, keepDiag = TRUE)$mat), pivot = TRUE)
      pivot<-attr(cv, "pivot")
      Z = t((cv[,order(pivot)]))
      c = Z%*%t(Z)
      x<- (matrix(1, nrow = N, ncol = 1) %*% mu) +(rand  %*% t(Z))
      }
  if (method == "eigen"){
      ed <- eigen(covar)
      Z <- ed$vectors%*%diag(sqrt(pmax(ed$values, 0)))
      c<-Z%*%t(Z)
      x<- (matrix(1, nrow = N, ncol = 1) %*% mu) + (rand  %*% t(Z))
  }
  data <- list("x" = x, "covar_input" = covar, "factorized_covar" = c, "Z_matrix" = t(Z))
  return(data)
}
