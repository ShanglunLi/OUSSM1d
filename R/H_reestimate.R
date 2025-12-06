#' Solve for H Matrix using NNLS
#'
#' Estimates the diagonal elements of the H matrix (observation variance) using
#' Non-Negative Least Squares (NNLS) to ensure positivity. It assumes a structure
#' where off-diagonal elements approximate the sum of diagonal components.
#'
#' @param result.H.mat Numeric matrix. The input H matrix to be solved/regularized.
#' @param k Integer. The dimension of the matrix (number of rows/columns).
#' @return A numeric matrix `result.H.mat.hat` with estimated diagonals and reconstructed off-diagonals.
#' @importFrom nnls nnls
#' @export
solve_H <- function(result.H.mat, k) {
  d <- dim(result.H.mat)
  
  Hij <- rep(NA, d[1] * (d[2] - 1) / 2)
  Hii <- diag(result.H.mat)
  A.mat <- matrix(0, nrow = d[1] * (d[2] - 1) / 2, ncol = length(Hii))
  
  count <- 1
  for (i in 1:(k-1)) {
    for (j in (i + 1):k) {
      Hij[count] <- result.H.mat[i,j]
      A.mat[count, i] <- 1
      A.mat[count, j] <- 1
      count <- count + 1
    }
  }
  
  fit <- nnls(A.mat, Hij)
  Hii.hat <- fit$x
  Hii.hat[Hii.hat == 0] <- 1e-06
  
  result.H.mat.hat <- result.H.mat
  diag(result.H.mat.hat) <- Hii.hat
  
  for (i in 1:k) {
    for (j in 1:k) {
      if (i != j){
        result.H.mat.hat[i,j] <- Hii.hat[i] + Hii.hat[j]
      }
    }
  }
  
  return(result.H.mat.hat)
}