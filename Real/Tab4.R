# ==============================================================================
## Frequentist DP-FRT calibration on Real Dataset
## ------------------------------------------------------------------------------
## Compares the two frequentist calibration methods (worst-case LFC and data-
## adaptive) on the full real-data cohort. Both calibrations sweep K in
## {0,...,n=n1+n0} with n_sim_calib Monte Carlo draws per K and call dp_frt()
## inside every inner iteration, so the total cost scales as
## O(n * n_sim_calib * dp_frt()).
# ==============================================================================

library(dplyr)
library(tibble)
library(matrixStats)

source("../DP-FRT.R")
source("../DP-Decision.R")

## ================================================================
## Dataset (same as Figures 6 and 7)
## ================================================================
n1  <- 7536
n0  <- 7540
n11 <- 569
n01 <- 590

## ================================================================
## Configuration
## ================================================================
eps_grid       <- c(0.2, 0.5, 1.0)   # privacy budgets to calibrate
n_sim_calib    <- 1000    # Monte Carlo draws per K; sweep covers all K in 0..n
R_calib        <- 2000   # MCMC sample size inside dp_frt()
alpha_test     <- 0.05
alpha_freq_use <- 0.05
zeta_use       <- 0.025  # confidence-set tolerance for data-adaptive

cat(sprintf("[Tab4] full cohort: n1=%d, n0=%d, n11=%d, n01=%d\n",
            n1, n0, n11, n01))

## ================================================================
## Stage-1 fit + both calibrations, per epsilon
## ================================================================
Tab4_rows <- lapply(eps_grid, function(eps_calib) {

  cat(sprintf("[Tab4] --- epsilon = %.2f ---\n", eps_calib))

  fit_calib <- dp_frt(
    n1 = n1, n0 = n0, n11 = n11, n01 = n01,
    epsilon  = eps_calib,
    prior    = list(type = "uniform"),
    method   = "sampling",
    R        = R_calib,
    alpha_CI = 0.05,
    seed     = 42,
    clip     = TRUE
  )

  dec_wc <- dp_frt_decision(
    frt_fit    = fit_calib,
    alpha      = alpha_test,
    framework  = "freq",
    calib      = "worst_case",
    alpha_freq = alpha_freq_use,
    n_sim      = n_sim_calib,
    seed       = 42
  )

  dec_da <- dp_frt_decision(
    frt_fit    = fit_calib,
    alpha      = alpha_test,
    framework  = "freq",
    calib      = "data_adaptive",
    alpha_freq = alpha_freq_use,
    zeta       = zeta_use,
    n_sim      = n_sim_calib,
    seed       = 42
  )

  tibble::tibble(
    epsilon     = eps_calib,
    alpha       = alpha_test,
    alpha_freq  = alpha_freq_use,
    calibration = c("Worst-case (LFC)", "Data-adaptive"),
    psi_obs     = round(c(dec_wc$psi_obs,   dec_da$psi_obs),   4),
    threshold   = round(c(dec_wc$threshold, dec_da$threshold), 4),
    decision    = c(dec_wc$decision, dec_da$decision)
  )
})

Tab4 <- dplyr::bind_rows(Tab4_rows)

cat("\n[Tab4] Frequentist calibration comparison on real data:\n")
print(Tab4)

# readr::write_csv(Tab4, "Tab4.csv")