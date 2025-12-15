# ==============================================================================
# Table 2: Coverage & Length of the Credible Set under DP-FRT
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
source("../DP-FRT.R")

################################################################################
################################################################################
################################################################################
## Utilitys

# (1) FRT right-tail p under sharp null given (a,b) and margins
frt_p <- function(a, b, n1, n0) {
  K <- a + b
  phyper(q = a - 1, m = K, n = (n1 + n0) - K, k = n1, lower.tail = FALSE)
}

# (2) One Monte Carlo run
sim_once_cover <- function(n11, n10, n01, n00, epsilon,
                           prior  = list(type = "uniform"),
                           alpha_CI  = 0.05,
                           method = "enumerate",
                           clip_tilde = FALSE) {
  n1 <- n11 + n10
  n0 <- n01 + n00
  p_true <- frt_p(a = n11, b = n01, n1 = n1, n0 = n0)
  
  out <- tryCatch({
    fit <- dp_frt(
      n1 = n1, n0 = n0, n11 = n11, n01 = n01, epsilon = epsilon,
      prior = prior, method = method, alpha_CI = alpha_CI, clip = clip_tilde
    )
    
    cred_set_vals <- as.numeric(fit$summaries$credible_set)
    ci_bounds     <- as.numeric(fit$summaries$credible_interval)
    
    post_mean <- as.numeric(fit$summaries$mean)
    bias <- post_mean - p_true
    
    p_key   <- round(p_true, 12)
    covered <- any(round(cred_set_vals, 12) == p_key)
    
    set_cardinality <- length(cred_set_vals)
    interval_width  <- diff(ci_bounds)
    
    c(covered = as.numeric(covered),
      set_cardinality = as.numeric(set_cardinality),
      interval_width  = as.numeric(interval_width),
      bias            = as.numeric(bias))
  }, error = function(e) {
    c(covered = NA_real_, 
      set_cardinality = NA_real_, 
      interval_width = NA_real_,
      bias = NA_real_)
  })
  
  as.numeric(out)
}

# (3) Replicate runs to compute coverage
sim_coverage <- function(n11, n10, n01, n00, epsilon, S = 1000,
                         prior  = list(type = "uniform"),
                         alpha_CI  = 0.05,
                         method = "enumerate") {
  res <- vapply(
    X = seq_len(S),
    FUN = function(i) sim_once_cover(n11, n10, n01, n00, epsilon,
                                     prior = prior, alpha_CI = alpha_CI, method = method),
    FUN.VALUE = c(covered = 0, set_cardinality = 0, interval_width = 0, bias = 0)
  )
  covered <- res["covered", ]
  set_len <- res["set_cardinality", ]
  int_w   <- res["interval_width", ]
  bias    <- res["bias", ]
  
  ok <- is.finite(covered)
  cov_hat <- mean(covered[ok], na.rm = TRUE)
  cov_se  <- sqrt(cov_hat * (1 - cov_hat) / sum(ok))
  
  ok_bias <- is.finite(bias)
  bias_hat <- mean(bias[ok_bias], na.rm = TRUE)
  
  c(
    coverage = cov_hat,
    coverage_se = cov_se,
    avg_set_cardinality = mean(set_len[is.finite(set_len)], na.rm = TRUE),
    avg_interval_width  = mean(int_w[is.finite(int_w)], na.rm = TRUE),
    avg_bias = bias_hat
  )
}

scen_defs <- tibble::tibble(
  scenario = factor(c("no_effect","small","medium","large"),
                    levels = c("no_effect","small","medium","large")),
  p_ctrl   = c(0.50, 0.50, 0.50, 0.50),
  p_treat  = c(0.50, 0.55, 0.60, 0.65),
  ratio    = c(1, 1, 1, 1)
)

sizes <- tibble::tibble(
  sample_size = factor(c("small","medium","large"),
                       levels = c("small","medium","large")),
  n_total     = c(100, 500, 1000)
)

design_grid <- tidyr::crossing(scen_defs, sizes) |>
  rowwise() |>
  mutate(
    n1  = round(n_total * ratio / (1 + ratio)),
    n0  = n_total - n1,
    n11 = round(n1 * p_treat),
    n01 = round(n0 * p_ctrl),
    n10 = n1 - n11,
    n00 = n0 - n01
  ) |>
  ungroup()

grouped_scenarios <- split(design_grid, design_grid$sample_size) |>
  purrr::map(\(df) {
    purrr::set_names(
      lapply(seq_len(nrow(df)), function(i) {
        r <- df[i, ]
        c(n11 = r$n11, n10 = r$n10, n01 = r$n01, n00 = r$n00)
      }),
      df$scenario
    )
  })

################################################################################
################################################################################
################################################################################
## Simulation Settings

set.seed(2025)
eps_grid   <- c(0.2, 0.5, 1)
S_reps     <- 1000
prior_use  <- list(type = "uniform")
alpha_use  <- 0.05
method_use <- "sampling"


table2_df <- imap_dfr(grouped_scenarios, function(group_list, size_label) {
  imap_dfr(group_list, function(tab, sc_name) {
    n11 <- tab["n11"]; n10 <- tab["n10"]; n01 <- tab["n01"]; n00 <- tab["n00"]
    map_dfr(eps_grid, function(ep) {
      stats <- sim_coverage(n11, n10, n01, n00, epsilon = ep, S = S_reps,
                            prior = prior_use, alpha_CI = alpha_use, method = method_use)
      tibble(
        Table = "Table 2",
        sample_size = size_label,
        scenario = sc_name,
        n11 = n11, n10 = n10, n01 = n01, n00 = n00,
        epsilon = ep,
        alpha_CI   = alpha_use,
        coverage = stats["coverage"],
        coverage_se = stats["coverage_se"],
        avg_set_cardinality = stats["avg_set_cardinality"],
        avg_interval_width  = stats["avg_interval_width"],
        avg_bias_posterior_mean  = stats["avg_bias"]
      )
    })
  })
}) |>
  mutate(
    sample_size = factor(sample_size, levels = c("small","medium","large"),
                         labels = c("n=100","n=500", "n=1000")),
    scenario = factor(scenario, levels = c("no_effect","small","medium","large"),
                      labels = c("Case 1: No effect",
                                 "Case 2: Small effect",
                                 "Case 3: Medium effect",
                                 "Case 4: Large effect"))
  ) |>
  relocate(Table)

print(table2_df, n = nrow(table2_df))

# readr::write_csv(table2_df, "Tab2.csv")