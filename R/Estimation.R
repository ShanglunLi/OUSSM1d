#' Estimate OUSSM Parameters (MLE)
#'
#' Performs Maximum Likelihood Estimation for the Ornstein-Uhlenbeck State-Space Model
#' using the Nelder-Mead optimization method.
#'
#' @param theta.init Numeric. Initial guess for the mean reversion rate.
#' @param ytn Numeric vector. The observed time series.
#' @param N Integer. The number of observations.
#' @param tn.diff Numeric vector. Time differences between observations.
#' @param H.par.init Numeric (optional). Initial guess for observation variance H.
#' @param Sigma.par.init Numeric (optional). Initial guess for state volatility sigma.
#' @param mu.init Numeric (optional). Initial guess for the long-run mean.
#' @return A list containing estimated parameters (H, mu, theta, sigma),
#'   derived quantities (Qtn, Qinf), the Hessian matrix, and the negative log-likelihood value.
#' @importFrom stats optim var
#' @export
estimation.OU <- function(theta.init, ytn, N, tn.diff, H.par.init = NULL, Sigma.par.init = NULL, mu.init = NULL) {
  if (is.null(H.par.init)) {
    H.par.init <- var(ytn[2:N] - ytn[1:(N - 1)]) / 3
  }
  if (is.null(Sigma.par.init)) {
    Sigma.par.init <- var(ytn[2:N] - ytn[1:(N - 1)]) / 3
  }
  if (is.null(mu.init)) {
    mu.init <- mean(ytn)
  }
  
  init.par <- c(log(theta.init), log(H.par.init), log(Sigma.par.init), mu.init)
  
  result.sim <- optim(init.par, N = N, ytn = ytn, tn.diff = tn.diff, neglogL, method = "Nelder-Mead", hessian = TRUE)
  theta <- exp(result.sim$par[1])
  H <- exp(result.sim$par[2])
  sigma <- exp(result.sim$par[3])
  mu <- result.sim$par[4]
  
  Qtn <- sigma * (1-exp(-2*theta)) / (2*theta)
  Qinf <- sigma / (2*theta)
  
  Hessian.record <- result.sim$hessian
  L.dist <- result.sim$value
  result.para.dist <- result.sim$par
  
  result <- list("H" = H, "mu" = mu, "theta" = theta, "sigma" = sigma, "Qtn" = Qtn, "Qinf" = Qinf, "Hessian.record" = Hessian.record,
                 "neglogL" = L.dist, "result.para.dist" = result.para.dist)
  
  return(result)
}

#' Estimate OUSSM Parameters with Regime-Switching Mean
#'
#' Performs MLE for an OUSSM where the mean shifts from mu1 to mu2 at `split.pt`.
#'
#' @param theta.init Numeric. Initial guess for theta.
#' @param ytn Numeric vector. Observed time series.
#' @param N Integer. Number of observations.
#' @param tn.diff Numeric vector. Time differences.
#' @param split.pt Integer. Index at which the regime shift occurs.
#' @return A list containing estimated parameters (H, mu1, mu2, theta, sigma) and optimization details.
#' @importFrom stats optim var
#' @export
estimation.OU.2mu <- function(theta.init, ytn, N, tn.diff, split.pt) {
  SigmaHinit <- var(ytn[2:N] - ytn[1:(N - 1)]) / 3
  H.par.init <- SigmaHinit
  Sigma.par.init <- SigmaHinit
  init.par <- c(log(theta.init), log(H.par.init), log(Sigma.par.init), mean(ytn[1:split.pt]), mean(ytn[(split.pt+1):length(ytn)]))
  
  result.sim <- optim(init.par, N = N, ytn = ytn, tn.diff = tn.diff, split.pt = split.pt, neglogL.2mu, method = "Nelder-Mead", hessian = TRUE)
  
  theta <- exp(result.sim$par[1])
  H <- exp(result.sim$par[2])
  sigma <- exp(result.sim$par[3])
  mu1 <- result.sim$par[4]
  mu2 <- result.sim$par[5]
  
  Qtn <- sigma * (1-exp(-2*theta)) / (2*theta)
  Qinf <- sigma / (2*theta)
  
  Hessian.record <- result.sim$hessian
  L.dist <- result.sim$value
  result.para.dist <- result.sim$par
  
  result <- list("H" = H, "mu1" = mu1, "mu2" = mu2, "theta" = theta, "sigma" = sigma, "Qtn" = Qtn, "Qinf" = Qinf, "Hessian.record" = Hessian.record,
                 "neglogL" = L.dist, "result.para.dist" = result.para.dist)
  
  return(result)
}

#' Estimate OUSSM Parameters (Fixed H)
#'
#' Performs MLE for OUSSM assuming the observation variance H is known/fixed.
#'
#' @param theta.init Numeric. Initial guess for theta.
#' @param ytn Numeric vector. Observed time series.
#' @param N Integer. Number of observations.
#' @param tn.diff Numeric vector. Time differences.
#' @param H Numeric. The fixed observation variance.
#' @param Sigma.par.init Numeric (optional). Initial guess for sigma.
#' @param mu.init Numeric (optional). Initial guess for mu.
#' @return A list containing estimated parameters (mu, theta, sigma) and optimization details.
#' @importFrom stats optim var
#' @export
estimation.OU.H.fixed <- function(theta.init, ytn, N, tn.diff, H, Sigma.par.init = NULL, mu.init = NULL) {
  if (is.null(Sigma.par.init)) {
    Sigma.par.init <- var(ytn[2:N] - ytn[1:(N - 1)]) / 3
  }
  if (is.null(mu.init)) {
    mu.init <- mean(ytn)
  }
  
  init.par <- c(log(theta.init), log(Sigma.par.init), mu.init)
  
  result.sim <- optim(init.par, H = H, N = N, ytn = ytn, tn.diff = tn.diff, neglogL.H.fixed, method = "Nelder-Mead", hessian = TRUE)
  theta <- exp(result.sim$par[1])
  sigma <- exp(result.sim$par[2])
  mu <- result.sim$par[3]
  
  Qtn <- sigma * (1-exp(-2*theta)) / (2*theta)
  Qinf <- sigma / (2*theta)
  
  Hessian.record <- result.sim$hessian
  L.dist <- result.sim$value
  result.para.dist <- result.sim$par
  
  result <- list("mu" = mu, "theta" = theta, "sigma" = sigma, "Qtn" = Qtn, "Qinf" = Qinf, "Hessian.record" = Hessian.record,
                 "neglogL" = L.dist, "result.para.dist" = result.para.dist)
  
  return(result)
}

#' Estimate OUSSM Parameters with Regime-Switching Mean (Fixed H)
#'
#' Performs MLE for Regime-Switching OUSSM assuming the observation variance H is known/fixed.
#'
#' @param theta.init Numeric. Initial guess for theta.
#' @param ytn Numeric vector. Observed time series.
#' @param N Integer. Number of observations.
#' @param tn.diff Numeric vector. Time differences.
#' @param split.pt Integer. Regime shift index.
#' @param H Numeric. The fixed observation variance.
#' @return A list containing estimated parameters (mu1, mu2, theta, sigma) and optimization details.
#' @importFrom stats optim var
#' @export
estimation.OU.H.fixed.2mu <- function(theta.init, ytn, N, tn.diff, split.pt, H) {
  Sigma.par.init <- var(ytn[2:N] - ytn[1:(N - 1)]) / 3
  init.par <- c(log(theta.init), log(Sigma.par.init), mean(ytn[1:split.pt]), mean(ytn[(split.pt + 1):length(ytn)]))
  
  result.sim <- optim(init.par, H = H, N = N, ytn = ytn, tn.diff = tn.diff, split.pt = split.pt, neglogL.H.fixed.2mu, method = "Nelder-Mead", hessian = TRUE)
  
  theta <- exp(result.sim$par[1])
  sigma <- exp(result.sim$par[2])
  mu1 <- result.sim$par[3]
  mu2 <- result.sim$par[4]
  
  Qtn <- sigma * (1-exp(-2*theta)) / (2*theta)
  Qinf <- sigma / (2*theta)
  
  Hessian.record <- result.sim$hessian
  L.dist <- result.sim$value
  result.para.dist <- result.sim$par
  
  result <- list("mu1" = mu1, "mu2" = mu2, "theta" = theta, "sigma" = sigma, "Qtn" = Qtn, "Qinf" = Qinf, "Hessian.record" = Hessian.record,
                 "neglogL" = L.dist, "result.para.dist" = result.para.dist)
  
  return(result)
}