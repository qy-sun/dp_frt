# ==============================================================================
# FigS5: Frequentist DP-FRT decision behavior (alpha_freq = 0.1)
#        Overlay of two statistics under the same Section 4.2 calibration:
#        Psi (DP-FRT-Bayes) and T = tilde n11/n1 - tilde n01/n0
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(ggplot2)
library(scales)

source("../DP-FRT.R")
source("../DP-Decision.R")

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
  N10 = c(  0,   0,   0,
            5,  25,  50,
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

## T statistic (difference in privatized proportions) + Freq calibration
## (geometric mechanism matches dp_frt in DP-FRT.R)

## Two-sided geometric mechanism
rgeom_two_sided <- function(size, rho) {
  rgeom(size, prob = 1 - rho) - rgeom(size, prob = 1 - rho)
}

## One privatized release of (n11, n01)
privatize_counts <- function(n11, n01, n1, n0, epsilon, clip = TRUE) {
  rho <- exp(-epsilon)
  t11 <- n11 + rgeom_two_sided(1, rho)
  t01 <- n01 + rgeom_two_sided(1, rho)
  if (clip) { t11 <- max(0, min(n1, t11)); t01 <- max(0, min(n0, t01)) }
  c(t11 = as.integer(t11), t01 = as.integer(t01))
}

## Test statistic
t_stat <- function(t11, t01, n1, n0) t11 / n1 - t01 / n0

## Vectorized draw of (t11, t01) and T under Q_K
t_sim_QK <- function(K, n1, n0, epsilon, n_sim, clip = TRUE) {
  n <- n1 + n0; rho <- exp(-epsilon)
  A   <- rhyper(n_sim, m = K, n = n - K, k = n1)
  t11 <- A       + rgeom_two_sided(n_sim, rho)
  t01 <- (K - A) + rgeom_two_sided(n_sim, rho)
  if (clip) { t11 <- pmax(0, pmin(n1, t11)); t01 <- pmax(0, pmin(n0, t01)) }
  list(T = t11 / n1 - t01 / n0, t11 = as.integer(t11), t01 = as.integer(t01))
}

## Worst-case threshold: max over K of the (1 - alpha_freq) quantile of T under Q_K
t_calibrate_worst_case <- function(n1, n0, epsilon, alpha_freq, n_sim = 100, clip = TRUE) {
  n <- n1 + n0
  K_vals <- 0:n
  tK <- numeric(length(K_vals))
  for (idx in seq_along(K_vals)) {
    K <- K_vals[idx]
    if (K == 0 || K == n) n_sim_K <- max(20, ceiling(n_sim / 20)) else n_sim_K <- n_sim
    sim <- t_sim_QK(K, n1, n0, epsilon, n_sim_K, clip)
    tK[idx] <- as.numeric(stats::quantile(sim$T, probs = 1 - alpha_freq, type = 1))
  }
  max(tK, na.rm = TRUE)
}

## Data-adaptive threshold: max over K in the (1 - zeta) confidence set of the
## (1 - alpha_freq + zeta) quantile of T (confidence set via a 1-D integer code)
t_calibrate_data_adaptive <- function(n1, n0, epsilon, alpha_freq, zeta,
                                      t11_obs, t01_obs, n_sim = 100, clip = TRUE) {
  n <- n1 + n0
  K_vals <- 0:n
  t1K     <- rep(NA_real_, length(K_vals))
  in_conf <- rep(FALSE,    length(K_vals))
  alpha1 <- alpha_freq - zeta
  if (alpha1 <= 0) stop("alpha_freq must be larger than zeta.")
  n0i      <- as.integer(n0)
  obs_code <- if (clip) as.integer(t11_obs) * (n0i + 1L) + as.integer(t01_obs) else NA_integer_
  obs_key  <- paste(t11_obs, t01_obs, sep = "_")
  for (idx in seq_along(K_vals)) {
    K <- K_vals[idx]
    sim <- t_sim_QK(K, n1, n0, epsilon, n_sim, clip)
    if (clip) {
      code <- as.integer(sim$t11) * (n0i + 1L) + as.integer(sim$t01)
      tab  <- sort(table(code), decreasing = TRUE)
      pos  <- match(as.character(obs_code), names(tab))
    } else {
      tab <- sort(table(paste(sim$t11, sim$t01, sep = "_")), decreasing = TRUE)
      pos <- match(obs_key, names(tab))
    }
    cum <- cumsum(as.numeric(tab)) / sum(tab)
    in_conf[idx] <- !is.na(pos) && (pos <= which(cum >= (1 - zeta))[1])
    t1K[idx] <- as.numeric(stats::quantile(sim$T, probs = 1 - alpha1, type = 1))
  }
  if (any(in_conf, na.rm = TRUE)) max(t1K[in_conf], na.rm = TRUE) else max(t1K, na.rm = TRUE)
}

## Reject iff T_obs exceeds the calibrated threshold
t_decision <- function(t11_obs, t01_obs, n1, n0, epsilon, alpha_freq, calib,
                       zeta, n_sim = 100, clip = TRUE) {
  T_obs  <- t_stat(t11_obs, t01_obs, n1, n0)
  t_star <- if (calib == "worst_case")
    t_calibrate_worst_case(n1, n0, epsilon, alpha_freq, n_sim, clip)
  else
    t_calibrate_data_adaptive(n1, n0, epsilon, alpha_freq, zeta, t11_obs, t01_obs, n_sim, clip)
  if (T_obs > t_star) "reject" else "not_reject"
}

## Single randomization + DP-FRT + Freq decision (Psi)
sim_once_decision_freq <- function(
    Y1, Y0,
    n1,
    epsilon,
    prior,
    alpha_test,
    alpha_freq,
    calib,
    n_sim_calib = 2000,
    method_use   = "sampling",
    alpha_CI_use = 0.05,
    clip_tilde   = TRUE
){
  n  <- length(Y1)
  n0 <- n - n1

  Z <- c(rep(1L, n1), rep(0L, n0))
  Z <- sample(Z)

  Y_obs <- ifelse(Z == 1L, Y1, Y0)

  n11 <- sum(Z == 1L & Y_obs == 1L)
  n10 <- sum(Z == 1L & Y_obs == 0L)
  n01 <- sum(Z == 0L & Y_obs == 1L)
  n00 <- sum(Z == 0L & Y_obs == 0L)

  fit <- dp_frt(
    n1 = n1, n0 = n0,
    n11 = n11, n01 = n01,
    epsilon  = epsilon,
    prior    = prior,
    method   = method_use,
    alpha_CI = alpha_CI_use,
    clip     = clip_tilde
  )

  dec <- dp_frt_decision(
    frt_fit    = fit,
    alpha      = alpha_test,
    framework  = "freq",
    calib      = calib,
    alpha_freq = alpha_freq,
    zeta        = alpha_freq / 10,
    n_sim      = n_sim_calib
  )

  dec$decision
}

## Repeat randomization + DP-FRT + Freq decision S times (Psi)
sim_freq_decision_one_setting <- function(
    Y1, Y0,
    n1,
    epsilon,
    freq_grid,
    S = 1000,
    prior        = list(type = "uniform"),
    alpha_test   = 0.05,
    method_use   = "sampling",
    alpha_CI_use = 0.05,
    clip_tilde   = TRUE,
    n_sim_calib  = 2000
){
  n_freq <- nrow(freq_grid)
  n      <- length(Y1)
  n0     <- n - n1

  decisions <- matrix(NA_character_, nrow = S, ncol = n_freq)
  colnames(decisions) <- freq_grid$calib_id

  for (s in seq_len(S)) {
    Z <- c(rep(1L, n1), rep(0L, n0))
    Z <- sample(Z)

    Y_obs <- ifelse(Z == 1L, Y1, Y0)

    n11 <- sum(Z == 1L & Y_obs == 1L)
    n10 <- sum(Z == 1L & Y_obs == 0L)
    n01 <- sum(Z == 0L & Y_obs == 1L)
    n00 <- sum(Z == 0L & Y_obs == 0L)

    fit <- dp_frt(
      n1 = n1, n0 = n0,
      n11 = n11, n01 = n01,
      epsilon  = epsilon,
      prior    = prior,
      method   = method_use,
      alpha_CI = alpha_CI_use,
      clip     = clip_tilde
    )

    for (j in seq_len(n_freq)) {
      calib_j      <- freq_grid$calib[j]
      alpha_freq_j <- freq_grid$alpha_freq[j]

      dec <- dp_frt_decision(
        frt_fit    = fit,
        alpha      = alpha_test,
        framework  = "freq",
        calib      = calib_j,
        alpha_freq = alpha_freq_j,
        zeta        = alpha_freq_j / 10,
        n_sim      = n_sim_calib
      )

      decisions[s, j] <- dec$decision
    }
  }

  map_dfr(seq_len(n_freq), function(j) {
    dec_j <- decisions[, j]
    tibble(
      method      = "Psi",
      calib_id    = freq_grid$calib_id[j],
      calib       = freq_grid$calib[j],
      alpha_freq  = freq_grid$alpha_freq[j],
      prop_reject = mean(dec_j == "reject")
    )
  })
}

## Repeat randomization + one DP release + T-based Freq decision S times (T)
sim_freq_T_decision_one_setting <- function(
    Y1, Y0,
    n1,
    epsilon,
    freq_grid,
    S = 1000,
    zeta,
    n_sim_calib = 2000,
    clip_tilde  = TRUE
){
  n_freq <- nrow(freq_grid)
  n      <- length(Y1)
  n0     <- n - n1

  decisions <- matrix(NA_character_, nrow = S, ncol = n_freq)
  colnames(decisions) <- freq_grid$calib_id

  for (s in seq_len(S)) {
    Z <- sample(c(rep(1L, n1), rep(0L, n0)))
    Y_obs <- ifelse(Z == 1L, Y1, Y0)

    n11 <- sum(Z == 1L & Y_obs == 1L)
    n01 <- sum(Z == 0L & Y_obs == 1L)

    tt   <- privatize_counts(n11, n01, n1, n0, epsilon, clip = clip_tilde)
    t11o <- as.integer(tt["t11"]); t01o <- as.integer(tt["t01"])

    for (j in seq_len(n_freq)) {
      decisions[s, j] <- t_decision(
        t11_obs = t11o, t01_obs = t01o, n1 = n1, n0 = n0, epsilon = epsilon,
        alpha_freq = freq_grid$alpha_freq[j], calib = freq_grid$calib[j],
        zeta = zeta, n_sim = n_sim_calib, clip = clip_tilde
      )
    }
  }

  map_dfr(seq_len(n_freq), function(j) {
    dec_j <- decisions[, j]
    tibble(
      method      = "T",
      calib_id    = freq_grid$calib_id[j],
      calib       = freq_grid$calib[j],
      alpha_freq  = freq_grid$alpha_freq[j],
      prop_reject = mean(dec_j == "reject")
    )
  })
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

alpha_freq_grid <- c(0.10)

freq_grid <- tidyr::crossing(
  tibble::tibble(
    calib = c("worst_case","data_adaptive")
  ),
  tibble::tibble(
    alpha_freq = alpha_freq_grid
  )
) %>%
  mutate(
    calib_id = paste0(calib, "_af", alpha_freq)
  ) %>%
  select(calib_id, calib, alpha_freq)

################################################################################
## Simulation Settings
################################################################################

set.seed(2025)
eps_grid     <- c(0.2, 0.5, 1.0)
S_reps       <- 100
prior_use    <- list(type = "uniform")
alpha_use    <- 0.05
method_use   <- "sampling"
alpha_CI_use <- 0.05
n_sim_calib  <- 100
clip_tilde   <- TRUE
zeta_use     <- alpha_freq_grid[1] / 10

design_grid <- tidyr::crossing(
  scen_defs,
  sizes,
  epsilon = eps_grid
) %>%
  mutate(
    n1 = round(n_total * ratio / (1 + ratio)),
    n0 = n_total - n1
  ) %>%
  left_join(science_N, by = c("scenario", "n_total"))

table_df <- design_grid %>%
  rowwise() %>%
  mutate(
    science = list(make_science_from_N(
      N11 = N11,
      N10 = N10,
      N01 = N01,
      N00 = N00
    )),
    sim_res = list(bind_rows(
      sim_freq_decision_one_setting(
        Y1           = science$Y1,
        Y0           = science$Y0,
        n1           = n1,
        epsilon      = epsilon,
        freq_grid    = freq_grid,
        S            = S_reps,
        prior        = prior_use,
        alpha_test   = alpha_use,
        method_use   = method_use,
        alpha_CI_use = alpha_CI_use,
        clip_tilde   = clip_tilde,
        n_sim_calib  = n_sim_calib
      ),
      sim_freq_T_decision_one_setting(
        Y1           = science$Y1,
        Y0           = science$Y0,
        n1           = n1,
        epsilon      = epsilon,
        freq_grid    = freq_grid,
        S            = S_reps,
        zeta         = zeta_use,
        n_sim_calib  = n_sim_calib,
        clip_tilde   = clip_tilde
      )
    ))
  ) %>%
  ungroup() %>%
  select(-science, -N11, -N10, -N01, -N00) %>%
  unnest(sim_res) %>%
  mutate(
    Table = "FigS5",
    sample_size = factor(paste0("n=", n_total),
                         levels = paste0("n=", c(100, 500, 1000))),
    scenario = factor(
      scenario,
      levels = c("no_effect","small","medium","large"),
      labels = c("Case 5: No effect",
                 "Case 6: Small effect",
                 "Case 7: Medium effect",
                 "Case 8: Large effect")
    ),
    metric = if_else(scenario == "Case 5: No effect", "Type I error", "Power")
  ) %>%
  relocate(Table)

## ================================================================
## Plot Figure S5: Psi vs T overlay
##   color/shape = calibration ; linetype = statistic (Psi solid, T dashed)
## ================================================================

plot_df_all <- table_df %>%
  mutate(
    epsilon    = as.numeric(epsilon),
    method     = factor(method, levels = c("Psi","T")),
    calib      = factor(
      calib,
      levels = c("worst_case","data_adaptive"),
      labels = c("Worst-case","Data-adaptive")
    ),
    sample_size = factor(
      sample_size,
      levels = paste0("n=", c(100, 500, 1000))
    ),
    scenario = factor(
      scenario,
      levels = c("Case 5: No effect",
                 "Case 6: Small effect",
                 "Case 7: Medium effect",
                 "Case 8: Large effect")
    )
  )

eps_grid <- sort(unique(plot_df_all$epsilon))

theme_fig_common <- theme_bw(base_size = 15) +
  theme(
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.0, "cm"),
    legend.key.width  = unit(0.9, "cm"),
    strip.background  = element_rect(fill = "gray90"),
    strip.text.x      = element_text(size = 13, face = "bold"),
    strip.text.y      = element_text(size = 13, face = "bold"),
    panel.grid.minor  = element_blank(),
    panel.grid        = element_blank(),
    axis.title.x      = element_text(margin = margin(t = 15))
  )

calib_cols   <- c("Worst-case" = "#0072B2", "Data-adaptive" = "#D55E00")
calib_shapes <- c("Worst-case" = 16,        "Data-adaptive" = 17)
stat_lty     <- c("Psi" = "solid", "T" = "dashed")

df_af <- plot_df_all %>%
  filter(alpha_freq == 0.10)

p_af_010 <- ggplot(
  df_af,
  aes(
    x = epsilon,
    y = prop_reject,
    color    = calib,
    shape    = calib,
    linetype = method,
    group    = interaction(calib, method)
  )
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  geom_hline(
    data = subset(df_af, scenario == "Case 5: No effect"),
    aes(yintercept = 0.10),
    color = "black", linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = eps_grid
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = percent
  ) +
  scale_color_manual(
    values = calib_cols,
    name   = "Calibration"
  ) +
  scale_shape_manual(
    values = calib_shapes,
    name   = "Calibration"
  ) +
  scale_linetype_manual(
    values = stat_lty,
    breaks = c("Psi", "T"),
    labels = expression(Psi, italic(T)),
    name   = "Statistic"
  ) +
  facet_grid(
    rows = vars(sample_size),
    cols = vars(scenario)
  ) +
  labs(
    x = expression(paste("Privacy budget ", epsilon)),
    y = "Reject rate"
  ) +
  guides(
    color    = guide_legend(order = 1),
    shape    = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  theme_fig_common

print(p_af_010)

# ggsave("FigS5.png", p_af_010, width = 10, height = 6, dpi = 1000)
