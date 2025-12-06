#' Calculate Negative Log-Likelihood for OUSSM
#'
#' This function calculates the negative log-likelihood of the Ornstein-Uhlenbeck
#' State-Space Model using the Kalman Filter.
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(theta)
#'     \item \code{para[2]}: log(H)
#'     \item \code{para[3]}: log(sigma)
#'     \item \code{para[4]}: mu
#'   }
#' @param N Integer. The total number of observations.
#' @param ytn Numeric vector. The observed values.
#' @param tn.diff Numeric vector. Time differences between observations.
#' @return Numeric. The negative log-likelihood value.
#' @importFrom hash hash
#' @export
neglogL <- function(para, N, ytn, tn.diff) {
  theta <- exp(para[1])
  H <- exp(para[2])
  sigma <- exp(para[3])
  mu <- para[4]
  
  atn <- rep(NA, N)
  Ptn <- rep(NA, N)
  atn[1] <- mu
  Ptn[1] <- sigma / (2*theta)
  
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  Ktn <- rep(NA, N)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, N)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, N)
  vtn[1] <- ytn[1] - mu - atn[1]
  
  for (i in 1:(N-1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn[i+1] - mu - atn[i+1]
  }
  
  # !!!!!!
  logL <- -(N + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  return(-logL)
}

#' Calculate Negative Log-Likelihood (Fixed H)
#'
#' Calculates -logL for OUSSM when the observation noise H is fixed/known.
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(theta)
#'     \item \code{para[2]}: log(sigma)
#'     \item \code{para[3]}: mu
#'   }
#' @param H Numeric. The fixed observation variance.
#' @param N Integer. Number of observations.
#' @param ytn Numeric vector. Observed values.
#' @param tn.diff Numeric vector. Time differences.
#' @return Numeric. The negative log-likelihood.
#' @importFrom hash hash
#' @export
neglogL.H.fixed <- function(para, H, N, ytn, tn.diff) {
  theta <- exp(para[1])
  sigma <- exp(para[2])
  mu <- para[3]
  
  atn <- rep(NA, N)
  Ptn <- rep(NA, N)
  atn[1] <- mu
  Ptn[1] <- sigma / (2*theta)
  
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  Ktn <- rep(NA, N)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, N)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, N)
  vtn[1] <- ytn[1] - mu - atn[1]
  
  for (i in 1:(N-1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn[i+1] - mu - atn[i+1]
  }
  
  # !!!!!!
  logL <- -(N + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  return(-logL)
}

#' Calculate Negative Log-Likelihood (Fixed Theta)
#'
#' Calculates -logL for OUSSM when the mean reversion rate theta is fixed/known.
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(H)
#'     \item \code{para[2]}: log(sigma)
#'     \item \code{para[3]}: mu
#'   }
#' @param theta Numeric. The fixed mean reversion rate.
#' @param N Integer. Number of observations.
#' @param ytn Numeric vector. Observed values.
#' @param tn.diff Numeric vector. Time differences.
#' @return Numeric. The negative log-likelihood.
#' @importFrom hash hash
#' @export
neglogL.theta.fixed <- function(para, theta, N, ytn, tn.diff) {
  H <- exp(para[1])
  sigma <- exp(para[2])
  mu <- para[3]
  
  atn <- rep(NA, N)
  Ptn <- rep(NA, N)
  atn[1] <- mu
  Ptn[1] <- sigma / (2*theta)
  
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  Ktn <- rep(NA, N)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, N)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, N)
  vtn[1] <- ytn[1] - mu - atn[1]
  
  for (i in 1:(N-1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn[i+1] - mu - atn[i+1]
  }
  
  # !!!!!!
  logL <- -(N + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  return(-logL)
}

#' Calculate Negative Log-Likelihood (Fixed Theta and H)
#'
#' Calculates -logL for OUSSM when both theta and H are fixed/known.
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(sigma)
#'     \item \code{para[2]}: mu
#'   }
#' @param H Numeric. Fixed observation variance.
#' @param theta Numeric. Fixed mean reversion rate.
#' @param N Integer. Number of observations.
#' @param ytn Numeric vector. Observed values.
#' @param tn.diff Numeric vector. Time differences.
#' @return Numeric. The negative log-likelihood.
#' @importFrom hash hash
#' @export
neglogL.theta.H.fixed <- function(para, H, theta, N, ytn, tn.diff) {
  sigma <- exp(para[1])
  mu <- para[2]
  
  atn <- rep(NA, N)
  Ptn <- rep(NA, N)
  atn[1] <- mu
  Ptn[1] <- sigma / (2*theta)
  
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  Ktn <- rep(NA, N)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, N)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, N)
  vtn[1] <- ytn[1] - mu - atn[1]
  
  for (i in 1:(N-1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn[i+1] - mu - atn[i+1]
  }
  
  # !!!!!!
  logL <- -(N + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  return(-logL)
}

#' Calculate Negative Log-Likelihood with Regime-Switching Mean
#'
#' Calculates -logL for OUSSM where the mean shifts from mu1 to mu2 at `split.pt`.
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(theta)
#'     \item \code{para[2]}: log(H)
#'     \item \code{para[3]}: log(sigma)
#'     \item \code{para[4]}: mu1
#'     \item \code{para[5]}: mu2
#'   }
#' @param N Integer. Number of observations.
#' @param ytn Numeric vector. Observed values.
#' @param tn.diff Numeric vector. Time differences.
#' @param split.pt Integer. Index where the regime shift occurs.
#' @return Numeric. The negative log-likelihood.
#' @importFrom hash hash
#' @export
neglogL.2mu <- function(para, N, ytn, tn.diff, split.pt) {
  theta <- exp(para[1])
  H <- exp(para[2])
  sigma <- exp(para[3])
  mu1 <- para[4]
  mu2 <- para[5]
  
  ytn1 <- ytn[1:split.pt]
  ytn2 <- ytn[(split.pt + 1):N]
  N1 <- split.pt
  N2 <- N - split.pt
  tn.diff1 <- tn.diff[1:(split.pt-1)]
  tn.diff2 <- tn.diff[(split.pt+1):(N-1)]
  
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  atn <- rep(NA, split.pt)
  Ptn <- rep(NA, split.pt)
  atn[1] <- mu1
  Ptn[1] <- sigma / (2*theta)
  
  Ktn <- rep(NA, split.pt)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, split.pt)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, split.pt)
  vtn[1] <- ytn1[1] - mu1 - atn[1]
  
  for (i in 1:(split.pt-1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff1[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff1[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff1[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff1[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn1[i+1] - mu1 - atn[i+1]
  }
  
  # !!!!!!
  logL.p1 <- -(split.pt + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  atn <- rep(NA, (N - split.pt))
  Ptn <- rep(NA, (N - split.pt))
  atn[1] <- mu2
  Ptn[1] <- sigma / (2*theta)
  
  Ktn <- rep(NA, (N - split.pt))
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, (N - split.pt))
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, (N - split.pt))
  vtn[1] <- ytn2[1] - mu2 - atn[1]
  
  for (i in 1:(N - split.pt - 1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff2[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff2[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff2[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff2[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn2[i+1] - mu2 - atn[i+1]
  }
  
  logL.p2 <- -(N - split.pt + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  logL <- logL.p1 + logL.p2
  
  return(-logL)
}


#' Calculate Negative Log-Likelihood with Regime-Switching Mean (Fixed H)
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(theta)
#'     \item \code{para[2]}: log(sigma)
#'     \item \code{para[3]}: mu1
#'     \item \code{para[4]}: mu2
#'   }
#' @param H Numeric. Fixed observation variance.
#' @param N Integer. Number of observations.
#' @param ytn Numeric vector. Observed values.
#' @param tn.diff Numeric vector. Time differences.
#' @param split.pt Integer. Split point index.
#' @return Numeric. The negative log-likelihood.
#' @importFrom hash hash
#' @export
neglogL.H.fixed.2mu <- function(para, H, N, ytn, tn.diff, split.pt) {
  theta <- exp(para[1])
  sigma <- exp(para[2])
  mu1 <- para[3]
  mu2 <- para[4]
  
  ytn1 <- ytn[1:split.pt]
  ytn2 <- ytn[(split.pt + 1):N]
  N1 <- split.pt
  N2 <- N - split.pt
  tn.diff1 <- tn.diff[1:(split.pt-1)]
  tn.diff2 <- tn.diff[(split.pt+1):(N-1)]
  
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  atn <- rep(NA, split.pt)
  Ptn <- rep(NA, split.pt)
  atn[1] <- mu1
  Ptn[1] <- sigma / (2*theta)
  
  Ktn <- rep(NA, split.pt)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, split.pt)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, split.pt)
  vtn[1] <- ytn1[1] - mu1 - atn[1]
  
  for (i in 1:(split.pt-1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff1[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff1[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff1[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff1[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn1[i+1] - mu1 - atn[i+1]
  }
  
  # !!!!!!
  logL.p1 <- -(split.pt + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  atn <- rep(NA, (N - split.pt))
  Ptn <- rep(NA, (N - split.pt))
  atn[1] <- mu2
  Ptn[1] <- sigma / (2*theta)
  
  Ktn <- rep(NA, (N - split.pt))
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, (N - split.pt))
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, (N - split.pt))
  vtn[1] <- ytn2[1] - mu2 - atn[1]
  
  for (i in 1:(N - split.pt - 1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff2[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff2[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff2[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff2[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn2[i+1] - mu2 - atn[i+1]
  }
  
  logL.p2 <- -(N - split.pt + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  logL <- logL.p1 + logL.p2
  
  return(-logL)
}

#' Calculate Negative Log-Likelihood with Regime-Switching Mean (Fixed Theta and H)
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(sigma)
#'     \item \code{para[2]}: mu1
#'     \item \code{para[3]}: mu2
#'   }
#' @param H Numeric. Fixed observation variance.
#' @param theta Numeric. Fixed mean reversion rate.
#' @param N Integer. Number of observations.
#' @param ytn Numeric vector. Observed values.
#' @param tn.diff Numeric vector. Time differences.
#' @param split.pt Integer. Split point index.
#' @return Numeric. The negative log-likelihood.
#' @importFrom hash hash
#' @export
neglogL.theta.H.fixed.2mu <- function(para, H, theta, N, ytn, tn.diff, split.pt) {
  sigma <- exp(para[1])
  mu1 <- para[2]
  mu2 <- para[3]
  
  ytn1 <- ytn[1:split.pt]
  ytn2 <- ytn[(split.pt + 1):N]
  N1 <- split.pt
  N2 <- N - split.pt
  tn.diff1 <- tn.diff[1:(split.pt-1)]
  tn.diff2 <- tn.diff[(split.pt+1):(N-1)]
  
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  atn <- rep(NA, split.pt)
  Ptn <- rep(NA, split.pt)
  atn[1] <- mu1
  Ptn[1] <- sigma / (2*theta)
  
  Ktn <- rep(NA, split.pt)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, split.pt)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, split.pt)
  vtn[1] <- ytn1[1] - mu1 - atn[1]
  
  for (i in 1:(split.pt-1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff1[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff1[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff1[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff1[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn1[i+1] - mu1 - atn[i+1]
  }
  
  # !!!!!!
  logL.p1 <- -(split.pt + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  atn <- rep(NA, (N - split.pt))
  Ptn <- rep(NA, (N - split.pt))
  atn[1] <- mu2
  Ptn[1] <- sigma / (2*theta)
  
  Ktn <- rep(NA, (N - split.pt))
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, (N - split.pt))
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, (N - split.pt))
  vtn[1] <- ytn2[1] - mu2 - atn[1]
  
  for (i in 1:(N - split.pt - 1)) {
    atn[i+1] <- Ctn.dic[[as.character(tn.diff2[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff2[i])]] * Ktn[i] * vtn[i]
    Ptn[i+1] <- Ctn.dic[[as.character(tn.diff2[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff2[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn2[i+1] - mu2 - atn[i+1]
  }
  
  logL.p2 <- -(N - split.pt + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  logL <- logL.p1 + logL.p2
  
  return(-logL)
}

#' Calculate Negative Log-Likelihood for Local Level Model
#'
#' This function calculates the negative log-likelihood of a Local Level Model
#' (Random Walk plus Noise) using the Kalman Filter.
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(H) (Observation variance)
#'     \item \code{para[2]}: log(sigma) (State variance)
#'     \item \code{para[3]}: mu (Initial state/intercept)
#'   }
#' @param N Integer. The total number of observations.
#' @param ytn Numeric vector. The observed values.
#' @param tn.diff Numeric vector. Time differences between observations.
#' @return Numeric. The negative log-likelihood value.
#' @importFrom hash hash
#' @export
neglogL.local <- function(para, N, ytn, tn.diff) {
  H <- exp(para[1])
  sigma <- exp(para[2])
  mu <- para[3]
  
  atn <- rep(NA, N)
  Ptn <- rep(NA, N)
  atn[1] <- mu
  Ptn[1] <- sigma
  
  unique.t <- sort(unique(tn.diff))
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*t
  }
  
  Ktn <- rep(NA, N)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, N)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, N)
  vtn[1] <- ytn[1] - mu - atn[1]
  
  for (i in 1:(N-1)) {
    atn[i+1] <- atn[i] + Ktn[i] * vtn[i]
    Ptn[i+1] <- Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn[i+1] - mu - atn[i+1]
  }
  
  # !!!!!!
  logL <- -(N + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  return(-logL)
}

#' Calculate Negative Log-Likelihood for Local Level Model with Regime-Switching Mean
#'
#' This function calculates the negative log-likelihood for a Local Level Model
#' where the mean shifts from `mu1` to `mu2` at a specified split point.
#'
#' @param para A numeric vector of parameters to be estimated:
#'   \itemize{
#'     \item \code{para[1]}: log(H)
#'     \item \code{para[2]}: log(sigma)
#'     \item \code{para[3]}: mu1 (Mean before split)
#'     \item \code{para[4]}: mu2 (Mean after split)
#'   }
#' @param N Integer. The total number of observations.
#' @param ytn Numeric vector. The observed values.
#' @param tn.diff Numeric vector. Time differences between observations.
#' @param split.pt Integer. The time index at which the mean shifts.
#' @return Numeric. The negative log-likelihood value.
#' @importFrom hash hash
#' @export
neglogL.local.2mu <- function(para, N, ytn, tn.diff, split.pt) {
  H <- exp(para[1])
  sigma <- exp(para[2])
  mu1 <- para[3]
  mu2 <- para[4]
  
  ytn1 <- ytn[1:split.pt]
  ytn2 <- ytn[(split.pt + 1):N]
  N1 <- split.pt
  N2 <- N - split.pt
  tn.diff1 <- tn.diff[1:(split.pt-1)]
  tn.diff2 <- tn.diff[(split.pt+1):(N-1)]
  
  unique.t <- sort(unique(tn.diff))
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*t
  }
  
  atn <- rep(NA, split.pt)
  Ptn <- rep(NA, split.pt)
  atn[1] <- mu1
  Ptn[1] <- sigma
  
  Ktn <- rep(NA, split.pt)
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, split.pt)
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, split.pt)
  vtn[1] <- ytn1[1] - mu1 - atn[1]
  
  for (i in 1:(split.pt-1)) {
    atn[i+1] <- atn[i] + Ktn[i] * vtn[i]
    Ptn[i+1] <- Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn1[i+1] - mu1 - atn[i+1]
  }
  
  # !!!!!!
  # Corrected from original snippet to use split.pt instead of N
  logL.p1 <- -(split.pt + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  
  atn <- rep(NA, (N-split.pt))
  Ptn <- rep(NA, (N-split.pt))
  atn[1] <- mu2
  Ptn[1] <- sigma
  
  Ktn <- rep(NA, (N-split.pt))
  Ktn[1] <- Ptn[1] / (Ptn[1] + H)
  Ftn <- rep(NA, (N-split.pt))
  Ftn[1] <- Ptn[1] + H
  vtn <- rep(NA, (N-split.pt))
  vtn[1] <- ytn2[1] - mu2 - atn[1]
  
  for (i in 1:(N-split.pt-1)) {
    atn[i+1] <- atn[i] + Ktn[i] * vtn[i]
    Ptn[i+1] <- Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
    Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
    Ftn[i+1] <- Ptn[i+1] + H
    vtn[i+1] <- ytn2[i+1] - mu2 - atn[i+1]
  }
  
  # !!!!!!
  # Corrected from original snippet to use (N - split.pt) instead of N
  logL.p2 <- -(N - split.pt + 1)/2 * log(2 * pi) - 1/2 * sum(log(Ftn) + vtn^2 / Ftn)
  logL <- logL.p1 + logL.p2
  
  return(-logL)
}