# OUSSM1d

The **OUSSM1d** package implements the framework described in **"State Space Modeling of the Ornstein-Uhlenbeck Process with Measurement Error: An Application to Microbiome Data"** and **"On the Optimal Sampling Scheme for Ornstein-Uhlenbeck State-Space Models"** by Shanglun Li, Toby Kenney, and Hong Gu.

It provides a likelihood-based framework grounded in Ornstein-Uhlenbeck (OU) state-space modeling to analyze longitudinal data, with a specific focus on microbiome log-ratio transformed abundances. Additionally, it provides tools to generate OUSSM trajectories under various sampling schemes (e.g., with and without repeated measurements) to validate the relationship between sampling intervals and estimation accuracy.

## Key Features

* **OUSSM Implementation:** Explicitly models measurement error ($\sigma_\epsilon^2$) distinct from process noise ($\sigma^2$), which is critical for noisy sequencing data.
* **Profile Likelihood Inference:** Constructs confidence intervals for the mean reversion rate ($\theta$) that are more reliable than asymptotic approximations for short time series.
* **Structural Break Testing:** Includes likelihood ratio tests to detect regime switches (e.g., changes in mean reversion rates or equilibrium levels) following disruptions like antibiotic use or diet changes.
* **Diagnostic Tools:** Checks for normality and independence of residuals to validate model assumptions.

## Installation

You can install the development version of OUSSM1d from [GitHub](https://github.com/ShanglunLi/OUSSM1d) with:

``` r
# install.packages("devtools")
devtools::install_github("ShanglunLi/OUSSM1d")
```

## Quick Start Example

This example demonstrates how to simulate an OUSSM process and estimate its parameters using Maximum Likelihood Estimation (MLE).

``` r
library(OUSSM1d)

# 1. Simulate Data (mimicking a microbiome log-ratio series)
set.seed(123)
N <- 100
tn.diff <- rep(1, N) # Regular time intervals (e.g., daily sampling)

# Parameters based on paper simulations
true_params <- list(H=0.5, theta=0.5, sigma=1.0, mu=0.0)

sim_data <- simulation.OU(
  H = true_params$H, 
  theta = true_params$theta, 
  sigma = true_params$sigma, 
  tn.diff = tn.diff, 
  mu = true_params$mu, 
  N = N, 
  sim.num = 1
)

y <- sim_data$ytn[1,] # Extract the first simulated path

# 2. Estimate Parameters (MLE)
est_result <- estimation.OU(
  theta.init = 1.0, # Initial guess
  ytn = y, 
  N = N, 
  tn.diff = tn.diff
)

# View estimated parameters
print(est_result$result.para.dist)
```


