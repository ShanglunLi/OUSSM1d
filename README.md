# OUSSM1d

The **OUSSM1d** package implements the framework described in **"State Space Modeling of the Ornstein-Uhlenbeck Process with Measurement Error: An Application to Microbiome Data"** by Shanglun Li, Toby Kenney, and Hong Gu.

It provides a likelihood-based framework grounded in Ornstein-Uhlenbeck (OU) state-space modeling to analyze longitudinal data, with a specific focus on microbiome log-ratio transformed abundances.

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
