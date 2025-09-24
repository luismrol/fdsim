#' Simulation algorithm for multivariate functional data based on Ieva and Paganoni (2017)
#'
#' @description
#' This function simulates multivariate functional data using Eigen, Cholesky and Singular
#' Value Decomposition, based on the method proposed by Ieva and Paganoni (2017) and implemented
#' initially in the library 'roahd', but with some changes enabling data exploration.
#'
#' @param N       Number of functional observations
#' @param mu      Matrix with L rows and P columns containing the mean curves for the L components at discrete grid points
#' @param grid    If you want to specify your own grid, the vector of the grid. Length must coincide with the number of columns in mu. Default is FALSE, in which case a [0,1] grid with even spacing
#'                over P points is implemented.
#' @param covar   List of L covariance matrices, one for each variable, each with P x P dimensions
#' @param rho     Vector of cross correlations among L components
#' @param method  Simulation method
#' @return A list of simulated functional data
#' @export


# This generates the man/*.Rd help files and updates NAMESPACE

mfd_sim <- function(N, mu, grid = FALSE, covar = "sq", rho = 0,
                    method = c("svd", "chol", "eigen"), delta = 0.5){
  if (!is.matrix(mu)){
    stop(paste("Mu should be a LxP matrix. An object of class", class(mu), "is not valid"))
  }
  if (!is.numeric(N)){
    stop(paste("Sample size must be a positive integer. Class", class(N)[1], "is not valid"))
  }
  # Columns for each dimension
  P <- ncol(mu)
  # Number of dimensions
  L <- nrow(mu)
  method <- method
  if (length(covar)==1 && is.matrix(covar)){
    covar = rep_len(covar, L)
  }
  # Create empty list to store each of the dimensions
  l_rand <- list()
  rand <- matrix(rnorm(N*P*L), N*P, L)
  if (isFALSE(grid)){
    grid = seq(0,1, by = 1/(P-1))
  }
  if(!is.list(covar)){
    covar1<-covar
    covar = list()
    for (i in 1:L){
      if (covar1[i] == "sq"){
        kernel <- function(a, b, delta) {
          as.matrix(outer(a, b, function(i, j) exp(-((i - j)^2) / ((delta)^2))))
        }
      }
      else if (covar1[i] == "abs"){
        kernel <- function(a, b, delta) {
          as.matrix(outer(a, b, function(i, j) exp(-(abs(i - j)) / ((delta)))))
        }
      }
      # if (isTRUE(het)){
      #  covar[[i]] = kernel(grid, grid, delta) - diag(abs((grid)-0.5))
      # }
      # else{
         covar[[i]] = kernel(grid, grid, delta)
      # }
    }
  }
  # Generate the correlation matrix between dimensions
  m_rho <- matrix(1, ncol = L, nrow = L)

  # Careful: The upper triangle vector should be inputed column-wise.

  m_rho[upper.tri(m_rho)] <- rho
  m_rho[lower.tri(m_rho)] <- rho

  # Generate random normal uncorrelated normal values as a base

  c<-covar

  if (method == "chol"){
    cholm_rho<-chol(Matrix::nearPD(m_rho)$mat, pivot = TRUE)
    f_rho<-t(cholm_rho[,attr(cholm_rho, "pivot")])
    # Apply correlation to the random values
    # For each variable, apply autocorrelation to the random values.
    for (i in 1:L){
      if(is.null(covar[[i]])){
        covar[[i]] <- diag(,P,P)
        message(paste("Running with Identity matrix of dimension", P, "as covariance for component", (i)))
      }
      # I used nearPD to force the positive-definiteness of the covariance function
      cv <- chol(as.matrix(Matrix::nearPD(covar[[i]], keepDiag = TRUE)$mat), pivot = TRUE)
      pivot<-attr(cv, "pivot")
      Z = t((cv[,order(pivot)]))
      c[[i]] = Z%*%t(Z)
      ed <- eigen(covar[[i]])
      l_rand[[i]]<- (matrix(1, nrow = N, ncol = 1) %*% mu[i,]) +(matrix(rand%*% t(f_rho), N, P)  %*% t(Z))
    }
  }
  if (method == "eigen"){
    e_rho <- eigen(m_rho)
    f_rho <- e_rho$vectors %*% diag(sqrt(pmax(e_rho$values, 0)))
    #Apply correlation to the random values
    for (i in 1:L){
      if(is.null(covar[[i]])){
        covar[[i]] <- diag(,P,P)
        message(paste("Running with Identity matrix of dimension", P, "as covariance for component", (i)))
      }
      ed <- eigen(covar[[i]])
      Z <- ed$vectors%*%diag(sqrt(pmax(ed$values, 0)))
      c[[i]]<-Z%*%t(Z)
      l_rand[[i]]<- (matrix(1, nrow = N, ncol = 1) %*% mu[i,]) + (matrix(rand%*% t(f_rho), N, P)  %*% t(Z))
    }
  }
  return(l_rand)
}
