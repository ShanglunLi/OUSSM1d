#' Calculate Residuals and Filtered Variances (OUSSM)
#'
#' Computes the prediction errors (vtn) and their variances (Ftn) for the OUSSM
#' using the Kalman Filter. Useful for residual analysis.
#'
#' @param para A numeric vector of parameters (log(theta), log(H), log(sigma), mu).
#' @param N Integer. Number of observations.
#' @param ytn Numeric vector. Observed time series.
#' @param tn.diff Numeric vector. Time differences.
#' @return A list containing `vtn` (residuals), `Ftn` (variances), and `atn` (states).
#' @importFrom hash hash
#' @export
cal.vtn.Ftn <- function(para, N, ytn, tn.diff) {
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
  
  result <- list("vtn" = vtn, "Ftn" = Ftn, "atn" = atn)
  
  return(result)
}

#' Calculate Residuals and Variances (Regime-Switching OUSSM)
#'
#' Computes prediction errors (vtn) and variances (Ftn) for the OUSSM with
#' a structural break in the mean.
#'
#' @param para A numeric vector of parameters (log(theta), log(H), log(sigma), mu1, mu2).
#' @param N Integer. Number of observations.
#' @param ytn Numeric vector. Observed time series.
#' @param tn.diff Numeric vector. Time differences.
#' @param split.pt Integer. Index of the regime shift.
#' @return A list containing `vtn` (residuals), `Ftn` (variances), and `atn` (states).
#' @importFrom hash hash
#' @export
cal.vtn.Ftn.2mu <- function(para, N, ytn, tn.diff, split.pt) {
  theta <- exp(para[1])
  H <- exp(para[2])
  sigma <- exp(para[3])
  mu1 <- para[4]
  mu2 <- para[5]
  
  atn <- rep(NA, N)
  Ptn <- rep(NA, N)
  atn[1] <- mu1
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
  vtn[1] <- ytn[1] - mu1 - atn[1]
  
  for (i in 1:(N-1)) {
    if (i <= split.pt) {
      atn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff[i])]] * Ktn[i] * vtn[i]
      Ptn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
      Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
      Ftn[i+1] <- Ptn[i+1] + H
      vtn[i+1] <- ytn[i+1] - mu1 - atn[i+1]
    }
    else {
      atn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]] * atn[i] + Ctn.dic[[as.character(tn.diff[i])]] * Ktn[i] * vtn[i]
      Ptn[i+1] <- Ctn.dic[[as.character(tn.diff[i])]]^2 * Ptn[i] * (1 - Ktn[i]) + Qtn.dic[[as.character(tn.diff[i])]]
      Ktn[i+1] <- Ptn[i+1] / (Ptn[i+1] + H)
      Ftn[i+1] <- Ptn[i+1] + H
      vtn[i+1] <- ytn[i+1] - mu2 - atn[i+1]
    }
    
  }
  
  result <- list("vtn" = vtn, "Ftn" = Ftn, "atn" = atn)
  
  return(result)
}

#' Diagnostic Test for OUSSM vs Local Level Model
#'
#' Performs a bootstrap Likelihood Ratio Test (LRT) to compare the OUSSM against
#' a Local Level (Random Walk + Noise) model. Also performs a Shapiro-Wilk test
#' on the standardized residuals.
#'
#' @param ytn Numeric vector. Observed time series.
#' @param tn.diff Numeric vector. Time differences.
#' @param sim.num Integer. Number of bootstrap simulations for the LRT.
#' @return A list containing:
#'   \item{test.result}{Boolean. TRUE if OUSSM is significantly better.}
#'   \item{observed.logL.diff}{Numeric. The observed difference in log-likelihoods.}
#'   \item{oussm.p}{Numeric. The p-value of the bootstrap test.}
#'   \item{shapiro.p_value}{Numeric. P-value of the residual normality test.}
#'   \item{vtn}{Vector. Standardized residuals.}
#' @importFrom stats optim var quantile shapiro.test
#' @export
diag.test <- function(ytn, tn.diff, sim.num = 100) {
  N = length(ytn)
  # Step 1: Compute observed log-likelihoods from SSM and Local model
  SigmaHinit <- var(diff(ytn)) / 3
  init.par <- c(log(0.1), log(SigmaHinit), log(SigmaHinit), mean(ytn))
  init.par2 <- c(log(SigmaHinit), log(SigmaHinit), mean(ytn))
  
  # Fit full SSM
  result.ssm <- optim(init.par, fn = neglogL, N = N, ytn = ytn, tn.diff = tn.diff, method = "Nelder-Mead")
  
  # Fit local model
  result.local <- optim(init.par2, fn = neglogL.local, N = N, ytn = ytn, tn.diff = tn.diff, method = "Nelder-Mead")
  
  L.obs.ssm <- result.ssm$value
  L.obs.local <- result.local$value
  observed.diff <- L.obs.local - L.obs.ssm
  
  # Step 2: Simulate under null (local model) and compute test statistic distribution
  H <- exp(result.local$par[1])
  sigma <- exp(result.local$par[2])
  mu <- result.local$par[3]
  log.ratios.null <- rep(NA, sim.num)
  count <- 0
  
  while (count < sim.num) {
    sim.result <- simulation.local(H, sigma, mu, N, sim.num = 1, tn.diff)
    ytn.sim <- sim.result$ytn[1,]
    
    SigmaHinit.sim <- var(diff(ytn.sim)) / 3
    init.sim.ssm <- c(log(0.1), log(SigmaHinit.sim), log(SigmaHinit.sim), mean(ytn.sim))
    init.sim.local <- c(log(SigmaHinit.sim), log(SigmaHinit.sim), mean(ytn.sim))
    
    error.ssm <- FALSE
    error.local <- FALSE
    
    tryCatch({
      fit.ssm <- optim(init.sim.ssm, fn = neglogL, N = N, ytn = ytn.sim, tn.diff = tn.diff, method = "Nelder-Mead")
    }, error = function(e) error.ssm <- TRUE)
    
    tryCatch({
      fit.local <- optim(init.sim.local, fn = neglogL.local, N = N, ytn = ytn.sim, tn.diff = tn.diff, method = "Nelder-Mead")
    }, error = function(e) error.local <- TRUE)
    
    if (!error.ssm && !error.local) {
      count <- count + 1
      log.ratios.null[count] <- fit.local$value - fit.ssm$value
    }
  }
  
  # Step 3: Likelihood ratio test
  test.result <- observed.diff > quantile(log.ratios.null, 0.95)
  oussm.p <- mean(log.ratios.null >= observed.diff)
  
  # Step 4: Residual diagnostics
  para.input <- result.ssm$par
  result.vtn <- cal.vtn.Ftn(para.input, N, ytn, tn.diff)
  
  # qqPlot(result.vtn$vtn, main = "Standardized residuals (vtn)")
  shapiro.result <- shapiro.test(result.vtn$vtn / sqrt(result.vtn$Ftn))
  
  return(list(
    test.result = test.result,
    observed.logL.diff = observed.diff,
    simulated.distribution = log.ratios.null,
    shapiro.p_value = shapiro.result$p.value,
    vtn = result.vtn$vtn,
    atn = result.vtn$atn,
    result.ssm = result.ssm,
    result.local = result.local,
    oussm.p = oussm.p
  ))
}

#' Diagnostic Test for Regime-Switching OUSSM
#'
#' Performs a bootstrap LRT to compare the Regime-Switching OUSSM against
#' a Regime-Switching Local Level model.
#'
#' @param ytn Numeric vector. Observed time series.
#' @param tn.diff Numeric vector. Time differences.
#' @param split.pt Integer. Index of the regime shift.
#' @param sim.num Integer. Number of bootstrap simulations.
#' @return A list containing test results and residual diagnostics.
#' @importFrom stats optim var quantile shapiro.test
#' @export
diag.test.2mu <- function(ytn, tn.diff, split.pt, sim.num = 100) {
  N = length(ytn)
  # Step 1: Compute observed log-likelihoods from SSM and Local model
  SigmaHinit <- var(diff(ytn)) / 3
  init.par <- c(log(0.1), log(SigmaHinit), log(SigmaHinit), mean(ytn[1:split.pt]), mean(ytn[split.pt:length(ytn)]))
  init.par2 <- c(log(SigmaHinit), log(SigmaHinit), mean(ytn[1:split.pt]), mean(ytn[split.pt:length(ytn)]))
  
  # Fit full SSM
  result.ssm <- optim(init.par, fn = neglogL.2mu, N = N, ytn = ytn, tn.diff = tn.diff, split.pt = split.pt, method = "Nelder-Mead")
  
  # Fit local model
  result.local <- optim(init.par2, fn = neglogL.local.2mu, N = N, ytn = ytn, tn.diff = tn.diff, split.pt = split.pt, method = "Nelder-Mead")
  
  L.obs.ssm <- result.ssm$value
  L.obs.local <- result.local$value
  observed.diff <- L.obs.local - L.obs.ssm
  
  # Step 2: Simulate under null (local model) and compute test statistic distribution
  H <- exp(result.local$par[1])
  sigma <- exp(result.local$par[2])
  mu1 <- result.local$par[3]
  mu2 <- result.local$par[4]
  log.ratios.null <- rep(NA, sim.num)
  count <- 0
  
  while (count < sim.num) {
    sim.result <- simulation.local.2mu(H, sigma, mu1, mu2, N, sim.num = 1, tn.diff, split.pt)
    ytn.sim <- sim.result$ytn[1,]
    
    SigmaHinit.sim <- var(diff(ytn.sim)) / 3
    init.sim.ssm <- c(log(0.1), log(SigmaHinit.sim), log(SigmaHinit.sim), mean(ytn.sim[1:split.pt]), mean(ytn.sim[split.pt:length(ytn)]))
    init.sim.local <- c(log(SigmaHinit.sim), log(SigmaHinit.sim), mean(ytn.sim[1:split.pt]), mean(ytn.sim[split.pt:length(ytn)]))
    
    error.ssm <- FALSE
    error.local <- FALSE
    
    tryCatch({
      fit.ssm <- optim(init.sim.ssm, fn = neglogL.2mu, N = N, ytn = ytn.sim, tn.diff = tn.diff, split.pt = split.pt, method = "Nelder-Mead")
    }, error = function(e) error.ssm <- TRUE)
    
    tryCatch({
      fit.local <- optim(init.sim.local, fn = neglogL.local.2mu, N = N, ytn = ytn.sim, tn.diff = tn.diff, split.pt = split.pt, method = "Nelder-Mead")
    }, error = function(e) error.local <- TRUE)
    
    if (!error.ssm && !error.local) {
      count <- count + 1
      log.ratios.null[count] <- fit.local$value - fit.ssm$value
    }
  }
  
  # Step 3: Likelihood ratio test
  test.result <- observed.diff > quantile(log.ratios.null, 0.95)
  oussm.p <- mean(log.ratios.null >= observed.diff)
  
  # Step 4: Residual diagnostics
  para.input <- result.ssm$par
  result.vtn <- cal.vtn.Ftn.2mu(para.input, N, ytn, tn.diff, split.pt)
  
  # qqPlot(result.vtn$vtn, main = "Standardized residuals (vtn)")
  shapiro.result <- shapiro.test(result.vtn$vtn / sqrt(result.vtn$Ftn))
  
  return(list(
    test.result = test.result,
    observed.logL.diff = observed.diff,
    simulated.distribution = log.ratios.null,
    shapiro.p_value = shapiro.result$p.value,
    vtn = result.vtn$vtn,
    atn = result.vtn$atn,
    result.ssm = result.ssm,
    result.local = result.local,
    oussm.p = oussm.p
  ))
}