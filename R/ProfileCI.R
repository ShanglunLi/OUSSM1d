#' Calculate Profile Likelihood CI for Theta (OUSSM)
#'
#' Computes the confidence interval for the mean reversion parameter (theta)
#' using the Profile Likelihood method.
#'
#' @param l.MLE Numeric. The maximum log-likelihood value (from the MLE result).
#' @param ytn Numeric vector. Observed time series.
#' @param N Integer. Number of observations.
#' @param tn.diff Numeric vector. Time differences between observations.
#' @param theta_MLE Numeric. The MLE point estimate for theta.
#' @param alpha Numeric. Significance level (default 0.05 for 95 percent CI).
#' @param grid.length Integer. Number of grid points for the profile search.
#' @param theta.factor Numeric. Multiplier to define the search range boundaries.
#' @return A list containing `lb` (lower bound) and `ub` (upper bound).
#' @importFrom stats qchisq optim var
#' @export
ProfileL.CI <- function(l.MLE, ytn, N, tn.diff, theta_MLE, alpha = 0.05, grid.length = 50, theta.factor = 5) {
  L0 <- l.MLE - 1/2 * qchisq((1 - alpha), df = 1)
  
  # Construct log-theta grids for increasing and decreasing theta
  log.theta.right <- seq(log(theta_MLE), log(theta_MLE * theta.factor), length.out = grid.length / 2)
  log.theta.left  <- seq(log(theta_MLE), log(theta_MLE / theta.factor), length.out = grid.length / 2)
  
  L.value.right <- rep(NA, length(log.theta.right))
  L.value.left <- rep(NA, length(log.theta.left))
  
  # Initialize optimization parameters
  Sigma.par.init <- var(diff(ytn)) / 3
  H.par.init <- var(diff(ytn)) / 3
  mean.ytn <- mean(ytn)
  init.par <- c(log(H.par.init), log(Sigma.par.init), mean.ytn)
  
  result.sim <- NULL
  
  # Search Right Side (Increasing `theta`)
  for (j in seq_along(log.theta.right)) {
    if (!is.null(result.sim) && !is.null(result.sim$convergence) && result.sim$convergence == 0) {
      init.par <- result.sim$par
    }
    
    result.sim <- tryCatch({
      optim(init.par, 
            theta = exp(log.theta.right[j]), N = N, ytn = ytn, tn.diff = tn.diff, fn = neglogL.theta.fixed, method = "Nelder-Mead")
    }, error = function(e) NULL)
    
    if (!is.null(result.sim) && result.sim$convergence == 0) {
      L.value.right[j] <- -result.sim$value
    } else {
      L.value.right[j] <- NA
    }
    
    if (!is.na(L.value.right[j]) && L.value.right[j] < L0 - 5) {  
      break 
    }
  }
  
  result.sim <- NULL
  init.par   <- c(log(H.par.init), log(Sigma.par.init), mean.ytn)
  
  # Search Left Side (Decreasing `theta`)
  for (j in seq_along(log.theta.left)) {
    if (!is.null(result.sim) && !is.null(result.sim$convergence) && result.sim$convergence == 0) {
      init.par <- result.sim$par
    }
    
    result.sim <- tryCatch({
      optim(init.par, 
            theta = exp(log.theta.left[j]), N = N, ytn = ytn, tn.diff = tn.diff, fn = neglogL.theta.fixed, method = "Nelder-Mead")
    }, error = function(e) NULL)
    
    if (!is.null(result.sim) && result.sim$convergence == 0) {
      L.value.left[j] <- -result.sim$value
    } else {
      L.value.left[j] <- NA
    }
    
    if (!is.na(L.value.left[j]) && L.value.left[j] < L0 - 5) {  
      break 
    }
  }
  
  valid.right.idx <- which(L.value.right >= L0 & !is.na(L.value.right))
  valid.left.idx <- which(L.value.left >= L0 & !is.na(L.value.left))
  
  if (length(valid.left.idx) > 0) {
    theta.lower <- exp(log.theta.left[max(valid.left.idx)]) 
  } else {
    theta.lower <- NA
  }
  
  if (length(valid.right.idx) > 0) {
    theta.upper <- exp(log.theta.right[max(valid.right.idx)]) 
  } else {
    theta.upper <- NA
  }
  
  lb <- theta.lower
  ub <- theta.upper
  
  result <- list("ub" = ub, "lb" = lb)
  
  return(result)
}

#' Calculate Profile Likelihood CI for Theta (Fixed H)
#'
#' Computes the confidence interval for theta for the OUSSM when H is fixed.
#'
#' @param l.MLE Numeric. Maximum log-likelihood.
#' @param ytn Numeric vector. Observed time series.
#' @param H Numeric. Fixed observation variance.
#' @param N Integer. Number of observations.
#' @param tn.diff Numeric vector. Time differences.
#' @param theta_MLE Numeric. MLE point estimate for theta.
#' @param alpha Numeric. Significance level.
#' @param grid.length Integer. Grid size.
#' @param theta.factor Numeric. Search range multiplier.
#' @return A list containing `lb` and `ub`.
#' @importFrom stats qchisq optim var
#' @export
ProfileL.CI.H.fixed <- function(l.MLE, ytn, H, N, tn.diff, theta_MLE, alpha = 0.05, grid.length = 50, theta.factor = 5) {
  L0 <- l.MLE - 1/2 * qchisq((1 - alpha), df = 1)
  
  # Construct log-theta grids for increasing and decreasing theta
  log.theta.right <- seq(log(theta_MLE), log(theta_MLE * theta.factor), length.out = grid.length / 2)
  log.theta.left  <- seq(log(theta_MLE), log(theta_MLE / theta.factor), length.out = grid.length / 2)
  
  L.value.right <- rep(NA, length(log.theta.right))
  L.value.left <- rep(NA, length(log.theta.left))
  
  # Initialize optimization parameters
  Sigma.par.init <- var(diff(ytn)) / 3
  mean.ytn <- mean(ytn)
  init.par <- c(log(Sigma.par.init), mean.ytn)
  
  result.sim <- NULL
  
  # Search Right Side (Increasing `theta`)
  for (j in seq_along(log.theta.right)) {
    if (!is.null(result.sim) && !is.null(result.sim$convergence) && result.sim$convergence == 0) {
      init.par <- result.sim$par
    }
    
    result.sim <- tryCatch({
      optim(init.par, 
            H = H, theta = exp(log.theta.right[j]), N = N, ytn = ytn, tn.diff = tn.diff, fn = neglogL.theta.H.fixed, method = "Nelder-Mead")
    }, error = function(e) NULL)
    
    if (!is.null(result.sim) && result.sim$convergence == 0) {
      L.value.right[j] <- -result.sim$value
    } else {
      L.value.right[j] <- NA
    }
    
    if (!is.na(L.value.right[j]) && L.value.right[j] < L0 - 5) {  
      break 
    }
  }
  
  result.sim <- NULL
  init.par   <- c(log(Sigma.par.init), mean.ytn)
  
  # Search Left Side (Decreasing `theta`)
  for (j in seq_along(log.theta.left)) {
    if (!is.null(result.sim) && !is.null(result.sim$convergence) && result.sim$convergence == 0) {
      init.par <- result.sim$par
    }
    
    result.sim <- tryCatch({
      optim(init.par, 
            H = H, theta = exp(log.theta.left[j]), N = N, ytn = ytn, tn.diff = tn.diff, fn = neglogL.theta.H.fixed, method = "Nelder-Mead")
    }, error = function(e) NULL)
    
    if (!is.null(result.sim) && result.sim$convergence == 0) {
      L.value.left[j] <- -result.sim$value
    } else {
      L.value.left[j] <- NA
    }
    
    if (!is.na(L.value.left[j]) && L.value.left[j] < L0 - 5) {  
      break 
    }
  }
  
  valid.right.idx <- which(L.value.right >= L0 & !is.na(L.value.right))
  valid.left.idx <- which(L.value.left >= L0 & !is.na(L.value.left))
  
  if (length(valid.left.idx) > 0) {
    theta.lower <- exp(log.theta.left[max(valid.left.idx)]) 
  } else {
    theta.lower <- NA
  }
  
  if (length(valid.right.idx) > 0) {
    theta.upper <- exp(log.theta.right[max(valid.right.idx)]) 
  } else {
    theta.upper <- NA
  }
  
  lb <- theta.lower
  ub <- theta.upper
  
  result <- list("ub" = ub, "lb" = lb)
  
  return(result)
}

#' Calculate Profile Likelihood CI for Theta (Fixed H, Regime-Switching Mean)
#'
#' Computes the confidence interval for theta for the Regime-Switching OUSSM when H is fixed.
#'
#' @param l.MLE Numeric. Maximum log-likelihood.
#' @param ytn Numeric vector. Observed time series.
#' @param split.pt Integer. Regime shift index.
#' @param H Numeric. Fixed observation variance.
#' @param N Integer. Number of observations.
#' @param tn.diff Numeric vector. Time differences.
#' @param theta_MLE Numeric. MLE point estimate for theta.
#' @param alpha Numeric. Significance level.
#' @param grid.length Integer. Grid size.
#' @param theta.factor Numeric. Search range multiplier.
#' @return A list containing `lb` and `ub`.
#' @importFrom stats qchisq optim var
#' @export
ProfileL.CI.H.fixed.2mu <- function(l.MLE, ytn, split.pt, H, N, tn.diff, theta_MLE, alpha = 0.05, grid.length = 50, theta.factor = 5) {
  L0 <- l.MLE - 1/2 * qchisq((1 - alpha), df = 1)
  
  # Construct log-theta grids for increasing and decreasing theta
  log.theta.right <- seq(log(theta_MLE), log(theta_MLE * theta.factor), length.out = grid.length / 2)
  log.theta.left  <- seq(log(theta_MLE), log(theta_MLE / theta.factor), length.out = grid.length / 2)
  
  L.value.right <- rep(NA, length(log.theta.right))
  L.value.left <- rep(NA, length(log.theta.left))
  
  # Initialize optimization parameters
  Sigma.par.init <- var(diff(ytn)) / 3
  init.par <- c(log(Sigma.par.init), mean(ytn[1:split.pt]), mean(ytn[split.pt:length(ytn)]))
  
  result.sim <- NULL
  
  # Search Right Side (Increasing `theta`)
  for (j in seq_along(log.theta.right)) {
    if (!is.null(result.sim) && !is.null(result.sim$convergence) && result.sim$convergence == 0) {
      init.par <- result.sim$par
    }
    
    result.sim <- tryCatch({
      optim(init.par, 
            H = H, theta = exp(log.theta.right[j]), N = N, ytn = ytn, tn.diff = tn.diff, split.pt = split.pt, fn = neglogL.theta.H.fixed.2mu, method = "Nelder-Mead")
    }, error = function(e) NULL)
    
    if (!is.null(result.sim) && result.sim$convergence == 0) {
      L.value.right[j] <- -result.sim$value
    } else {
      L.value.right[j] <- NA
    }
    
    if (!is.na(L.value.right[j]) && L.value.right[j] < L0 - 5) {  
      break 
    }
  }
  
  result.sim <- NULL
  # FIXED: Re-initialize correctly for the 2-mu model
  init.par <- c(log(Sigma.par.init), mean(ytn[1:split.pt]), mean(ytn[split.pt:length(ytn)]))
  
  # Search Left Side (Decreasing `theta`)
  for (j in seq_along(log.theta.left)) {
    if (!is.null(result.sim) && !is.null(result.sim$convergence) && result.sim$convergence == 0) {
      init.par <- result.sim$par
    }
    
    result.sim <- tryCatch({
      optim(init.par, 
            H = H, theta = exp(log.theta.left[j]), N = N, ytn = ytn, tn.diff = tn.diff, split.pt = split.pt, fn = neglogL.theta.H.fixed.2mu, method = "Nelder-Mead")
    }, error = function(e) NULL)
    
    if (!is.null(result.sim) && result.sim$convergence == 0) {
      L.value.left[j] <- -result.sim$value
    } else {
      L.value.left[j] <- NA
    }
    
    if (!is.na(L.value.left[j]) && L.value.left[j] < L0 - 5) {  
      break 
    }
  }
  
  valid.right.idx <- which(L.value.right >= L0 & !is.na(L.value.right))
  valid.left.idx <- which(L.value.left >= L0 & !is.na(L.value.left))
  
  if (length(valid.left.idx) > 0) {
    theta.lower <- exp(log.theta.left[max(valid.left.idx)]) 
  } else {
    theta.lower <- NA
  }
  
  if (length(valid.right.idx) > 0) {
    theta.upper <- exp(log.theta.right[max(valid.right.idx)]) 
  } else {
    theta.upper <- NA
  }
  
  lb <- theta.lower
  ub <- theta.upper
  
  result <- list("ub" = ub, "lb" = lb)
  
  return(result)
}