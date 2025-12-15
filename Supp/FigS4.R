# ==============================================================================
# FigS4: Frequentist DP-FRT decision behavior (alpha_freq = 0.1)
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

## Single randomization + DP-FRT + Freq decision
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
    zeta        = alpha_freq / 2,
    n_sim      = n_sim_calib
  )
  
  dec$decision
}

## Repeat randomization + DP-FRT + Freq decision S times
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
        zeta        = alpha_freq_j / 2,
        n_sim      = n_sim_calib
      )
      
      decisions[s, j] <- dec$decision
    }
  }
  
  map_dfr(seq_len(n_freq), function(j) {
    dec_j <- decisions[, j]
    tibble(
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
    sim_res = list(
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
        clip_tilde   = TRUE,
        n_sim_calib  = n_sim_calib
      )
    )
  ) %>%
  ungroup() %>%
  select(-science, -N11, -N10, -N01, -N00) %>%
  unnest(sim_res) %>%
  mutate(
    Table = "FigS4",
    sample_size = factor("n=500", levels = "n=500"),
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
## Plot Figure S4
## ================================================================

plot_df_all <- table_df %>%
  mutate(
    epsilon    = as.numeric(epsilon),
    calib      = factor(
      calib,
      levels = c("worst_case","data_adaptive"),
      labels = c("Worst-case","Data-adaptive")
    ),
    sample_size = factor(
      sample_size,
      levels = c("n=500")
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
    strip.background  = element_rect(fill = "gray90"),
    strip.text.x      = element_text(size = 13, face = "bold"),
    strip.text.y      = element_text(size = 13, face = "bold"),
    panel.grid.minor  = element_blank(),
    panel.grid        = element_blank(),
    axis.title.x      = element_text(margin = margin(t = 15))
  )

calib_cols <- c(
  "Worst-case"    = "#0072B2",
  "Data-adaptive" = "#D55E00"
)

df_af <- plot_df_all %>%
  filter(alpha_freq == 0.10)

p_af_010 <- ggplot(
  df_af,
  aes(
    x = epsilon,
    y = prop_reject,
    color = calib,
    shape = calib,
    group = calib
  )
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
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
  facet_grid(
    rows = vars(sample_size),
    cols = vars(scenario)
  ) +
  labs(
    x = expression(paste("Privacy budget ", epsilon)),
    y = "Reject rate",
    title = bquote(alpha[Freq] == 0.10)
  ) +
  scale_shape_manual(
    values = c(
      "Worst-case"    = 16,
      "Data-adaptive" = 17
    ),
    name = "Calibration"
  ) +
  theme_fig_common

print(p_af_010)

# ggsave("FigS4.png", p_af_010, width = 10, height = 6, dpi = 1000)