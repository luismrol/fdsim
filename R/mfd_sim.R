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
#' @param het     Logical, if a heteroskedasticity element should be taken into account. If TRUE, the main diagonal of the covariance
#'                matrix will minimize at
#' @param covar   List of L covariance matrices, one for each variable, each with P x P dimensions
#' @param rho     Vector of cross correlations among L components
#' @param method  Simulation method
#' @return A list of simulated functional data
#' @export


 # This generates the man/*.Rd help files and updates NAMESPACE

mfd_sim <- function(N, mu, grid = FALSE, covar = "sq", rho = 0,
                    method = c("svd", "chol", "eigen"), delta = 0.5, het = TRUE){
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
  if (length(covar)==1){
    covar = rep_len(covar, L)
  }
  # Create empty list to store each of the dimensions
  l_rand <- list()
  if (isFALSE(grid)){
    grid = seq(0,1, by = 1/(P-1))
  }
  if(!is.list(covar)){
    covar1<-covar
    covar = list()
    for (i in 1:L){
      if (covar1[i] == "sq"){
        kernel <- function(a, b, delta) {
          as.matrix(outer(a, b, function(i, j) 2*exp(-((i - j)^2) / ((delta)^2))))
        }
      }
      else if (covar1[i] == "abs"){
        kernel <- function(a, b, delta) {
          as.matrix(outer(a, b, function(i, j) 2*exp(-(abs(i - j)) / ((delta)))))
        }
      }
      if (isTRUE(het)){
       covar[[i]] = kernel(grid, grid, delta) - diag(abs((grid)-0.5))
      }
      else{
        covar[[i]] = kernel(grid, grid, delta)
      }
    }
  }
  # Generate the correlation matrix between dimensions
  m_rho <- matrix(1, ncol = L, nrow = L)

  # Careful: The upper triangle vector should be inputed column-wise.

  m_rho[upper.tri(m_rho)] <- rho
  m_rho[lower.tri(m_rho)] <- rho

  # Generate random normal uncorrelated normal values as a base
  rand <- matrix(rnorm(N*P*L), N*P, L)

  if (method == "chol"){

    # Apply correlation to the random values
    m_rand <- rand %*% chol(m_rho)

    # For each variable, apply autocorrelation to the random values.
    for (i in 1:L){
      if(is.null(covar[[i]])){
        covar[[i]] <- diag(,P,P)
        message(paste("Running with Identity matrix of dimension", P, "as covariance for component", (i)))
      }
      # I used nearPD to force the positive-definiteness of the covariance function
      cv <- chol(as.matrix(Matrix::nearPD(covar[[i]], keepDiag = TRUE)$mat))
      l_rand[[i]]<- (matrix(1, nrow = N, ncol = 1) %*% mu[i,]) +(matrix(m_rand[,i], N, P) %*% cv)
    }
  }
  if (method == "svd"){

    # Apply correlation to the random values
    m_rand <- rand %*% t(svd(m_rho)$v %*% (t(svd(m_rho)$u) * sqrt(pmax(svd(m_rho)$d, 0))))
    #print(str(m_rand))
    for (i in 1:L){
      if(is.null(covar[[i]])){
        covar[[i]] <- diag(,P,P)
        message(paste("Running with Identity matrix of dimension", P, "as covariance for component", (i)))
      }
      svd <- svd(covar[[i]])
      cv <- t(svd$v %*% (t(svd$u) * sqrt(pmax(svd$d, 0))))
      l_rand[[i]]<- (matrix(1, nrow = N, ncol = 1) %*% mu[i,])+ (matrix(m_rand[,i], N, P) %*% cv)
    }
  }
  if (method == "eigen"){
    e_rho <- eigen(m_rho, symmetric <- TRUE)
    m_rand <- rand %*% (e_rho$vectors %*% diag(sqrt(pmax (e_rho$values, 0))) %*% t(e_rho$vectors))
    for (i in 1:L){
      if(is.null(covar[[i]])){
        covar[[i]] <- diag(,P,P)
        message(paste("Running with Identity matrix of dimension", P, "as covariance for component", (i)))
      }
      ed <- eigen(covar[[i]], symmetric <- TRUE)
      cv <- (ed$vectors %*% diag(sqrt(pmax (ed$values, 0))) %*% t(ed$vectors))
      l_rand[[i]]<- (matrix(1, nrow = N, ncol = 1) %*% mu[i,]) + (matrix(m_rand[,i], N, P) %*% cv)
    }
  }
  return(l_rand)
}
