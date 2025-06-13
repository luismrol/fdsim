#' Simulation algorithm for multivariate functional data from Ieva and Paganoni (2017)
#'
#' @description
#' This function simulates multivariate functional data using Eigen, Cholesky and Singular
#' Value Decomposition, based on the method proposed by Ieva and Paganoni (2017) and implemented
#' In the library 'roahd
#'
#' @param N       Number of functional observations
#' @param mu      Matrix with L rows and P columns containing the mean curves for the L components at discrete grid points
#' @param covar   List of L covariance matrices, one for each variable, each with P x P dimensions
#' @param rho     Vector of cross correlations among L components
#' @param method  Simulation method
#' @return A list of simulated functional data
#' @export


 # This generates the man/*.Rd help files and updates NAMESPACE

mfd_sim <- function(N, mu, covar = NULL, rho = 0,
                    method = c("svd", "chol", "eigen")){
  if (!is.matrix(mu)){
    stop(paste("Mu should be a LxP matrix. An object of class", class(mu), "is not valid"))
  }
  if (!is.null(covar) && !is.list(covar)){
    stop(paste("If provided, Covariance matrix must be a list of symmetric, positive definite matrices. Class", class(covar)[1], "is not valid"))
  }
  if (!is.numeric(N)){
    stop(paste("Sample size must be a positive integer. Class", class(N)[1], "is not valid"))
  }
  # Columns for each dimension
  P <- ncol(mu)
  N <- floor(N)
  # Number of dimensions
  L <- nrow(mu)
  method <- match.arg(method, c("svd", "chol", "eigen"))

  # Create empty list to store each of the dimensions
  l_rand <- list()

  if (is.null(covar)) {
    covar <- list()
    for (i in 1:L){
      covar[i]<-list(NULL)
    }
  }

  # Generate the correlation matrix between dimensions
  m_rho <- matrix(1, ncol <- L, nrow <- L)
  m_rho[upper.tri(m_rho)] <- rho
  m_rho[lower.tri(m_rho)] <- rho



  # Generate random uncorrelated normal values as a base
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
      cv <- chol(covar[[i]])
      l_rand[[i]]<- (matrix(1, nrow = N, ncol = 1) %*% mu) +(matrix(m_rand[,i], N, P) %*% cv)
    }
  }
  if (method == "svd"){

    # Apply correlation to the random values
    m_rand <- rand %*% t(svd(m_rho)$v %*% (t(svd(m_rho)$u) * sqrt(pmax(svd(m_rho)$d, 0))))

    for (i in 1:L){
      if(is.null(covar[[i]])){
        covar[[i]] <- diag(,P,P)
        message(paste("Running with Identity matrix of dimension", P, "as covariance for component", (i)))
      }
      svd <- svd(covar[[i]])
      cv <- t(svd$v %*% (t(svd$u) * sqrt(pmax(svd$d, 0))))
      l_rand[[i]]<- (matrix(1, nrow = N, ncol = 1) %*% mu)+ (matrix(m_rand[,i], N, P) %*% cv)
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
      l_rand[[i]]<- (matrix(1, nrow = N, ncol = 1) %*% mu) + (matrix(m_rand[,i], N, P) %*% cv)
    }
  }
  return(l_rand)
}
