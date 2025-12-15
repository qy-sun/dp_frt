# ==============================================================================
# Table S2: DP-FRT abstain behavior with DP-TopUp
# Fixed first-stage epsilon = 0.5
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)

source("../DP-FRT.R")
source("../DP-Decision.R")
source("../DP-TopUp.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------------------------
# Science table (fixed potential outcomes)
# ------------------------------------------------------------------------------

science_N <- tibble::tibble(
  scenario = factor(
    rep(c("no_effect","small","medium","large"), each = 3),
    levels = c("no_effect","small","medium","large")
  ),
  n_total  = rep(c(100, 500, 1000), times = 4),
  N11 = c( 50, 250, 500,
           50, 250, 500,
           50, 250, 500,
           50, 250, 500),
  N10 = c(  0,    0,   0,
            5,   25,  50,
            10,  50, 100,
            15,  75, 150),
  N01 = c(  0,   0,   0,
            0,   0,   0,
            0,   0,   0,
            0,   0,   0),
  N00 = c( 50, 250, 500,
           45, 225, 450,
           40, 200, 400,
           35, 175, 350)
)

make_science_from_N <- function(N11, N10, N01, N00) {
  Y1 <- c(rep(1L, N11), rep(1L, N10), rep(0L, N01), rep(0L, N00))
  Y0 <- c(rep(1L, N11), rep(0L, N10), rep(1L, N01), rep(0L, N00))
  list(Y1 = Y1, Y0 = Y0)
}

# ------------------------------------------------------------------------------
# One replicate: stage-1 DP-FRT + epsilon lower bound + DP-TopUp decisions
# ------------------------------------------------------------------------------

sim_one_rep_topup <- function(Y1, Y0,
                              n1,
                              epsilon1 = 0.5,
                              alpha = 0.05,
                              xi = 0.05,
                              prior = list(type = "uniform"),
                              method_use = "sampling",
                              alpha_CI_use = 0.05,
                              clip1 = TRUE,
                              clip2 = TRUE,
                              lambda_0 = 1,
                              lambda_1 = 1,
                              lambda_u = NULL,
                              mults = c(1, 2, 5, 10),
                              seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n  <- length(Y1)
  n0 <- n - n1
  
  Z <- c(rep(1L, n1), rep(0L, n0))
  Z <- sample(Z)
  
  Y_obs <- ifelse(Z == 1L, Y1, Y0)
  
  n11 <- sum(Z == 1L & Y_obs == 1L)
  n01 <- sum(Z == 0L & Y_obs == 1L)
  
  # --- stage 1 DP-FRT ---
  fit1 <- dp_frt(
    n1 = n1, n0 = n0,
    n11 = n11, n01 = n01,
    epsilon  = epsilon1,
    prior    = prior,
    method   = method_use,
    alpha_CI = alpha_CI_use,
    clip     = clip1
  )
  
  dec1 <- dp_frt_decision(
    frt_fit       = fit1,
    alpha         = alpha,
    framework     = "bayes",
    bayes_abstain = TRUE,
    lambda_0      = lambda_0,
    lambda_1      = lambda_1,
    lambda_u      = lambda_u
  )
  
  d1 <- dec1$decision %||% NA_character_
  
  eps_lb <- NA_real_
  d2_vec <- rep(NA_character_, length(mults))
  names(d2_vec) <- paste0("m", mults)
  
  # --- top-up branch ---
  if (!is.na(d1) && d1 == "abstain") {
    
    lb_obj <- dp_frt_epsilon_plus_lower_bound(
      decision1 = dec1,
      frt_fit   = fit1,
      xi        = xi,
      alpha     = alpha
    )
    
    eps_lb <- lb_obj$epsilon_plus_lower_bound %||% NA_real_
    
    if (is.finite(eps_lb) && eps_lb > 0) {
      
      for (k in seq_along(mults)) {
        
        eps_plus <- mults[k] * eps_lb
        
        top <- dp_frt_topup_release(
          n11 = n11, n01 = n01,
          n1 = n1, n0 = n0,
          epsilon_plus = eps_plus,
          clip = clip2
        )
        
        if ("clip_plus" %in% names(formals(dp_frt_update_with_topup_and_decide))) {
          res2 <- dp_frt_update_with_topup_and_decide(
            decision1 = dec1,
            tilde_plus = top$tilde_plus,
            epsilon_plus = eps_plus,
            frt_fit = fit1,
            clip_plus = top$clip %||% clip2,
            alpha = alpha,
            bayes_abstain = TRUE,
            lambda_0 = lambda_0,
            lambda_1 = lambda_1,
            lambda_u = lambda_u,
            method = method_use,
            alpha_CI = alpha_CI_use
          )
        } else {
          res2 <- dp_frt_update_with_topup_and_decide(
            decision1 = dec1,
            tilde_plus = top$tilde_plus,
            epsilon_plus = eps_plus,
            frt_fit = fit1,
            alpha = alpha,
            bayes_abstain = TRUE,
            lambda_0 = lambda_0,
            lambda_1 = lambda_1,
            lambda_u = lambda_u,
            method = method_use,
            alpha_CI = alpha_CI_use
          )
        }
        
        d2_vec[k] <- res2$decision2$decision %||% NA_character_
      }
    }
  }
  
  tibble(
    decision1 = d1,
    eps_lb = eps_lb,
    !!!as.list(d2_vec)
  )
}

# ------------------------------------------------------------------------------
# Summarize one design setting for Table S2
# ------------------------------------------------------------------------------

summarize_one_setting_TabS2 <- function(scenario,
                                        sample_size,
                                        n_total,
                                        N11, N10, N01, N00,
                                        n1,
                                        epsilon1 = 0.5,
                                        S = 1000,
                                        alpha = 0.05,
                                        xi = 0.05,
                                        prior = list(type = "uniform"),
                                        method_use = "sampling",
                                        alpha_CI_use = 0.05,
                                        clip1 = TRUE,
                                        clip2 = TRUE,
                                        lambda_0 = 1,
                                        lambda_1 = 1,
                                        lambda_u = NULL,
                                        mults = c(1, 2, 5, 10),
                                        seed = 2025) {
  
  sci <- make_science_from_N(N11, N10, N01, N00)
  
  set.seed(seed)
  
  reps <- map_dfr(seq_len(S), function(s) {
    sim_one_rep_topup(
      Y1 = sci$Y1,
      Y0 = sci$Y0,
      n1 = n1,
      epsilon1 = epsilon1,
      alpha = alpha,
      xi = xi,
      prior = prior,
      method_use = method_use,
      alpha_CI_use = alpha_CI_use,
      clip1 = clip1,
      clip2 = clip2,
      lambda_0 = lambda_0,
      lambda_1 = lambda_1,
      lambda_u = lambda_u,
      mults = mults,
      seed = sample.int(1e9, 1)
    )
  })
  
  abstain_stage1 <- mean(reps$decision1 == "abstain", na.rm = TRUE)
  
  eps_lb_vec <- reps$eps_lb[reps$decision1 == "abstain"]
  eps_lb_vec <- eps_lb_vec[is.finite(eps_lb_vec) & eps_lb_vec > 0]
  
  eps_lb_median <- if (length(eps_lb_vec) == 0) NA_real_ else median(eps_lb_vec)
  eps_lb_mean   <- if (length(eps_lb_vec) == 0) NA_real_ else mean(eps_lb_vec)
  eps_lb_sd     <- if (length(eps_lb_vec) == 0) NA_real_ else stats::sd(eps_lb_vec)
  
  eligible <- (reps$decision1 == "abstain") &
    is.finite(reps$eps_lb) & reps$eps_lb > 0
  
  out <- tibble(
    Case = scenario,
    sample_size = sample_size,
    n_total = n_total,
    epsilon = epsilon1,
    abstain_stage1 = abstain_stage1,
    eps_lb_median = eps_lb_median,
    eps_lb_mean = eps_lb_mean,
    eps_lb_sd = eps_lb_sd
  )
  
  for (m in mults) {
    col <- paste0("m", m)
    out[[paste0("abstain_after_topup_", m, "x")]] <-
      mean(reps[[col]][eligible] == "abstain", na.rm = TRUE)
  }
  
  for (m in mults) {
    a_col <- paste0("abstain_after_topup_", m, "x")
    out[[paste0("escape_after_topup_", m, "x")]] <- 1 - out[[a_col]]
  }
  
  out
}

# ------------------------------------------------------------------------------
# Design grid and Table S2
# ------------------------------------------------------------------------------

set.seed(2025)

epsilon1_fixed <- 0.5
S_reps <- 1000

lambda_0_use <- 1
lambda_1_use <- 1
lambda_u_use <- NULL

design_grid <- science_N %>%
  mutate(
    n1 = n_total / 2,
    sample_size = factor(
      n_total,
      levels = c(100, 500, 1000),
      labels = c("n=100","n=500","n=1000")
    ),
    Case = factor(
      scenario,
      levels = c("no_effect","small","medium","large"),
      labels = c("Case 5: No effect",
                 "Case 6: Small effect",
                 "Case 7: Medium effect",
                 "Case 8: Large effect")
    )
  )

TabS2 <- pmap_dfr(
  design_grid %>%
    select(Case, sample_size, n_total,
           n1, N11, N10, N01, N00) %>%
    as.list(),
  function(Case, sample_size, n_total,
           n1, N11, N10, N01, N00) {
    
    summarize_one_setting_TabS2(
      scenario = Case,
      sample_size = sample_size,
      n_total = n_total,
      N11 = N11, N10 = N10, N01 = N01, N00 = N00,
      n1 = n1,
      epsilon1 = epsilon1_fixed,
      S = S_reps,
      lambda_0 = lambda_0_use,
      lambda_1 = lambda_1_use,
      lambda_u = lambda_u_use,
      seed = 2025 + n_total
    )
  }
) %>%
  arrange(sample_size, Case) %>%
  mutate(
    across(
      c(starts_with("abstain"), starts_with("escape")),
      ~ round(100 * .x, 1)
    ),
    eps_lb_median = round(eps_lb_median, 3),
    eps_lb_mean   = round(eps_lb_mean, 3),
    eps_lb_sd     = round(eps_lb_sd, 3)
  )

print(TabS2)

# ------------------------------------------------------------------------------
# Export CSV
# ------------------------------------------------------------------------------

readr::write_csv(TabS2, "TabS2.csv")