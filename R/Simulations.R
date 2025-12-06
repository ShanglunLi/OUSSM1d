#' Simulate Ornstein-Uhlenbeck State-Space Model
#'
#' This function generates simulated paths for an Ornstein-Uhlenbeck state-space model (OUSSM).
#' It simulates both the latent state process (atn) and the observed process (ytn)
#' based on the specified parameters.
#'
#' @param H Numeric. The observation variance (noise).
#' @param theta Numeric. The mean reversion rate of the OU process.
#' @param sigma Numeric. The volatility (diffusion) parameter of the state process.
#' @param tn.diff Numeric vector. The time differences between observations.
#' @param mu Numeric. The long-run mean of the process.
#' @param N Integer. The number of time points per simulation.
#' @param sim.num Integer. The number of independent simulation paths to generate.
#' @return A list containing:
#' \item{ytn}{A matrix of dimension (sim.num x N) containing the observed values.}
#' \item{atn}{A matrix of dimension (sim.num x N) containing the latent state values.}
#' @importFrom hash hash
#' @importFrom stats rnorm
#' @export
simulation.OU <- function(H, theta, sigma, tn.diff, mu, N, sim.num) {
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  P1 <- sigma / (2*theta)
  
  ytn <- matrix(rep(NA,N*sim.num), nrow = sim.num)
  atn <- matrix(rep(NA,N*sim.num), nrow = sim.num)
  
  # simulate ytn
  for (i in 1:sim.num) {
    # simulate epsilon_tn
    epsilon.tn <- rnorm(n = N, 0, sqrt(H))
    
    atn[i,1] <- rnorm(n = 1, 0, sqrt(P1)) # at0
    for (tn in 1:N) {
      ytn[i,tn] <- mu + atn[i,tn] + epsilon.tn[tn]
      
      if (tn < N) {
        atn[i,(tn + 1)] <- Ctn.dic[[as.character(tn.diff[tn])]] * atn[i,tn] + rnorm(n = 1, mean = 0, sd = sqrt(Qtn.dic[[as.character(tn.diff[tn])]]))
      }
    }
  }
  result <- list("ytn" = ytn, "atn" = atn)
  
  return(result)
}

#' Simulate OUSSM with Regime-Switching Mean
#'
#' This function simulates an Ornstein-Uhlenbeck state-space model (OUSSM) where the
#' long-run mean shifts from `mu1` to `mu2` at a specified split point.
#'
#' @param H Numeric. The observation variance.
#' @param theta Numeric. The mean reversion rate.
#' @param sigma Numeric. The volatility parameter.
#' @param tn.diff Numeric vector. Time differences between observations.
#' @param mu1 Numeric. The long-run mean before the split point.
#' @param mu2 Numeric. The long-run mean after the split point.
#' @param N Integer. The total number of time points.
#' @param sim.num Integer. The number of simulation paths.
#' @param split.pt Integer. The time index at which the mean shifts from `mu1` to `mu2`.
#' @return A list containing `ytn` (observed matrix) and `atn` (latent state matrix).
#' @importFrom hash hash
#' @importFrom stats rnorm
#' @export
simulation.OU.2mu <- function(H, theta, sigma, tn.diff, mu1, mu2, N, sim.num, split.pt) {
  unique.t <- sort(unique(tn.diff))
  
  Ctn.dic <- hash()
  for (t in unique.t) {
    Ctn.dic[[as.character(t)]] <- exp(-theta*t)
  }
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*(1-exp(-2*theta*t)) / (2*theta)
  }
  
  P1 <- sigma / (2*theta)
  
  ytn <- matrix(rep(NA,N*sim.num), nrow = sim.num)
  atn <- matrix(rep(NA,N*sim.num), nrow = sim.num)
  
  # simulate ytn
  for (i in 1:sim.num) {
    # simulate epsilon_tn
    epsilon.tn <- rnorm(n = N, 0, sqrt(H))
    
    atn[i,1] <- rnorm(n = 1, 0, sqrt(P1)) # at0
    for (tn in 1:N) {
      if (tn <= split.pt) {
        ytn[i,tn] <- mu1 + atn[i,tn] + epsilon.tn[tn]
        
        if (tn < N) {
          atn[i,(tn + 1)] <- Ctn.dic[[as.character(tn.diff[tn])]] * atn[i,tn] + rnorm(n = 1, mean = 0, sd = sqrt(Qtn.dic[[as.character(tn.diff[tn])]]))
        }
      }
      else {
        ytn[i,tn] <- mu2 + atn[i,tn] + epsilon.tn[tn]
        
        if (tn < N) {
          atn[i,(tn + 1)] <- Ctn.dic[[as.character(tn.diff[tn])]] * atn[i,tn] + rnorm(n = 1, mean = 0, sd = sqrt(Qtn.dic[[as.character(tn.diff[tn])]]))
        }
      }
    }
  }
  result <- list("ytn" = ytn, "atn" = atn)
  
  return(result)
}


#' Simulate Local Level State-Space Model
#'
#' This function generates simulated paths for a Local Level state-space model
#' (Random Walk plus Noise).
#'
#' @param H Numeric. The observation variance.
#' @param sigma Numeric. The volatility (diffusion) of the random walk.
#' @param mu Numeric. The intercept/mean parameter.
#' @param N Integer. The number of time points.
#' @param sim.num Integer. The number of simulation paths.
#' @param tn.diff Numeric vector. Time differences between observations.
#' @return A list containing `ytn` (observed matrix) and `atn` (latent state matrix).
#' @importFrom hash hash
#' @importFrom stats rnorm
#' @export
simulation.local <- function(H, sigma, mu, N, sim.num, tn.diff) {
  unique.t <- sort(unique(tn.diff))
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*t
  }
  
  ytn <- matrix(rep(NA,N*sim.num), nrow = sim.num)
  atn <- matrix(rep(NA,N*sim.num), nrow = sim.num)
  
  # simulate ytn
  for (i in 1:sim.num) {
    # simulate epsilon_tn
    epsilon.tn <- rnorm(n = N, 0, sqrt(H))
    
    atn[i,1] <- rnorm(n = 1, 0, sqrt(sigma)) # at0
    for (tn in 1:N) {
      ytn[i,tn] <- mu + atn[i,tn] + epsilon.tn[tn]
      
      if (tn < N) {
        atn[i,(tn + 1)] <- atn[i,tn] + rnorm(n = 1, 0, sqrt(Qtn.dic[[as.character(tn.diff[tn])]]))
      }
    }
  }
  result <- list("ytn" = ytn, "atn" = atn)
  
  return(result)
}

#' Simulate Local Level Model with Regime-Switching Mean
#'
#' This function simulates a Local Level state-space model where the
#' mean shifts from `mu1` to `mu2` at a specified split point.
#'
#' @param H Numeric. The observation variance.
#' @param sigma Numeric. The volatility of the random walk.
#' @param mu1 Numeric. The mean before the split point.
#' @param mu2 Numeric. The mean after the split point.
#' @param N Integer. The total number of time points.
#' @param sim.num Integer. The number of simulation paths.
#' @param tn.diff Numeric vector. Time differences between observations.
#' @param split.pt Integer. The time index at which the mean shifts.
#' @return A list containing `ytn` (observed matrix) and `atn` (latent state matrix).
#' @importFrom hash hash
#' @importFrom stats rnorm
#' @export
simulation.local.2mu <- function(H, sigma, mu1, mu2, N, sim.num, tn.diff, split.pt) {
  unique.t <- sort(unique(tn.diff))
  
  Qtn.dic <- hash()
  for (t in unique.t) {
    Qtn.dic[[as.character(t)]] <- sigma*t
  }
  
  ytn <- matrix(rep(NA,N*sim.num), nrow = sim.num)
  atn <- matrix(rep(NA,N*sim.num), nrow = sim.num)
  
  # simulate ytn
  for (i in 1:sim.num) {
    # simulate epsilon_tn
    epsilon.tn <- rnorm(n = N, 0, sqrt(H))
    
    atn[i,1] <- rnorm(n = 1, 0, sqrt(sigma)) # at0
    for (tn in 1:N) {
      if (tn <= split.pt) {
        ytn[i,tn] <- mu1 + atn[i,tn] + epsilon.tn[tn]
        if (tn < N) {
          atn[i,(tn + 1)] <- atn[i,tn] + rnorm(n = 1, 0, sqrt(Qtn.dic[[as.character(tn.diff[tn])]]))
        }
      }
      else {
        ytn[i,tn] <- mu2 + atn[i,tn] + epsilon.tn[tn]
        if (tn < N) {
          atn[i,(tn + 1)] <- atn[i,tn] + rnorm(n = 1, 0, sqrt(Qtn.dic[[as.character(tn.diff[tn])]]))
        }
      }
    }
  }
  result <- list("ytn" = ytn, "atn" = atn)
  
  return(result)
}