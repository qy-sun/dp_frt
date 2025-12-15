# ==============================================================================
# Figure 4 & Figure 5: Bayesian DP-FRT decision behavior
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(readr)
library(scales)
library(ggpattern)

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

## Single randomization + DP-FRT + Bayes decision
sim_once_decision <- function(Y1, Y0,
                              n1,
                              epsilon,
                              prior,
                              alpha_test,
                              lambda_0, lambda_1, lambda_u,
                              method_use = "sampling",
                              alpha_CI_use = 0.05,
                              clip_tilde = TRUE) {
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
    frt_fit       = fit,
    alpha         = alpha_test,
    framework     = "bayes",
    bayes_abstain = TRUE,
    lambda_0      = lambda_0,
    lambda_1      = lambda_1,
    lambda_u      = lambda_u
  )
  
  dec$decision
}

## Repeat randomization + DP-FRT + Bayes decision S times

sim_decision_one_setting <- function(Y1, Y0,
                                     n1,
                                     epsilon,
                                     lambda_grid,
                                     S = 1000,
                                     prior       = list(type = "uniform"),
                                     alpha_test  = 0.05,
                                     method_use  = "sampling",
                                     alpha_CI_use = 0.05,
                                     clip_tilde  = TRUE) {
  n_lambda <- nrow(lambda_grid)
  n        <- length(Y1)
  n0       <- n - n1
  
  decisions <- matrix(NA_character_, nrow = S, ncol = n_lambda)
  colnames(decisions) <- lambda_grid$lambda_id
  
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
    
    for (j in seq_len(n_lambda)) {
      lam0 <- lambda_grid$lambda_0[j]
      lam1 <- lambda_grid$lambda_1[j]
      lamu <- lambda_grid$lambda_u[j]
      
      dec <- dp_frt_decision(
        frt_fit       = fit,
        alpha         = alpha_test,
        framework     = "bayes",
        bayes_abstain = TRUE,
        lambda_0      = lam0,
        lambda_1      = lam1,
        lambda_u      = lamu
      )
      decisions[s, j] <- dec$decision
    }
  }
  
  map_dfr(seq_len(n_lambda), function(j) {
    dec_j <- decisions[, j]
    tibble(
      lambda_id       = lambda_grid$lambda_id[j],
      lambda_0        = lambda_grid$lambda_0[j],
      lambda_1        = lambda_grid$lambda_1[j],
      lambda_u        = lambda_grid$lambda_u[j],
      prop_reject     = mean(dec_j == "reject"),
      prop_not_reject = mean(dec_j == "not_reject"),
      prop_abstain    = mean(dec_j == "abstain")
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

lambda_grid <- tibble::tibble(
  lambda_id       = c("sym_no_abstain",
                      "sym_wide_abstain",
                      "sym_medium_abstain",
                      "penalize_aggressive",
                      "penalize_conservative"),
  lambda_0        = c(1,   1,   1,   0.5, 2),
  lambda_1        = c(1,   1,   1,   2,   0.5),
  lambda_u_frac_H = c(0.6, 0.025, 0.10, 0.025, 0.025)
) %>%
  mutate(
    H        = 2 * lambda_0 * lambda_1 / (lambda_0 + lambda_1),
    lambda_u = lambda_u_frac_H * H
  )

################################################################################
################################################################################
################################################################################
## Simulation Settings

set.seed(2025)
eps_grid     <- c(0.2, 0.5, 1.0)
S_reps       <- 100
prior_use    <- list(type = "uniform")
alpha_use    <- 0.05
method_use   <- "sampling"
alpha_CI_use <- 0.05

design_grid <- tidyr::crossing(
  scen_defs,
  sizes,
  epsilon = eps_grid
) %>%
  mutate(
    n1 = round(n_total * ratio / (1 + ratio)),
    n0 = n_total - n1
  ) %>%
  left_join(science_N,
            by = c("scenario", "n_total"))

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
      sim_decision_one_setting(
        Y1          = science$Y1,
        Y0          = science$Y0,
        n1          = n1,
        epsilon     = epsilon,
        lambda_grid = lambda_grid,
        S           = S_reps,
        prior       = prior_use,
        alpha_test  = alpha_use,
        method_use  = method_use,
        alpha_CI_use = alpha_CI_use,
        clip_tilde  = TRUE
      )
    )
  ) %>%
  ungroup() %>%
  select(-science, -N11, -N10, -N01, -N00) %>%
  unnest(sim_res) %>%
  mutate(
    Table = "Fig4Fig5",
    sample_size = factor(sample_size,
                         levels = c("small","medium","large"),
                         labels = c("n=100","n=500","n=1000")),
    scenario = factor(scenario,
                      levels = c("no_effect","small","medium","large"),
                      labels = c("Case 5: No effect",
                                 "Case 6: Small effect",
                                 "Case 7: Medium effect",
                                 "Case 8: Large effect"))
  ) %>%
  relocate(Table)

plot_df_all <- table_df %>%
  mutate(
    sample_size = factor(sample_size,
                         levels = c("n=100","n=500","n=1000"),
                         labels = c("n=100","n=500","n=1000")),
    scenario = factor(scenario,
                      levels = c("Case 5: No effect",
                                 "Case 6: Small effect",
                                 "Case 7: Medium effect",
                                 "Case 8: Large effect")),
    epsilon   = as.numeric(epsilon),
    lambda_id = factor(lambda_id)
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

## ================================================================
## Figure 4: stacked bar plot
## ================================================================
lambda_id_for_stack <- levels(plot_df_all$lambda_id)[1]
message("Figure 4 uses stacked bars with lambda_id = ", lambda_id_for_stack)

df_stack <- plot_df_all %>%
  filter(lambda_id == lambda_id_for_stack) %>%
  select(scenario, sample_size, epsilon,
         prop_reject, prop_not_reject, prop_abstain) %>%
  pivot_longer(
    cols      = starts_with("prop_"),
    names_to  = "decision",
    values_to = "proportion"
  ) %>%
  mutate(
    decision = factor(
      decision,
      levels = c("prop_reject","prop_abstain","prop_not_reject"),
      labels = c("Reject","Abstain","Not Reject")
    ),
    epsilon = factor(epsilon,
                     levels = eps_grid,
                     labels = eps_grid)
  )

fig4 <- ggplot(
  df_stack,
  aes(
    x = epsilon,
    y = proportion,
    fill = decision,
    pattern = decision
  )
) +
  geom_col_pattern(
    width = 0.75,
    position = position_stack(),
    color = "black",
    linewidth = 0.2,
    pattern_fill = "black",
    pattern_color = NA,
    pattern_density = 0.2,
    pattern_spacing = 0.02
  ) +
  scale_fill_manual(
    values = c(
      "Reject"     = scales::alpha("black", 0.7),
      "Abstain"    = scales::alpha("black", 0.3),
      "Not Reject" = scales::alpha("black", 0)
    )
  ) +
  scale_pattern_manual(
    values = c(
      "Reject"     = "stripe",
      "Abstain"    = "circle",
      "Not Reject" = "none"
    )
  ) +
  facet_grid(rows = vars(sample_size), cols = vars(scenario)) +
  labs(
    x = expression(paste("Privacy budget ", epsilon)),
    y = "Decision proportion",
    fill    = "Decision",
    pattern = "Decision"
  ) +
  theme_fig_common +
  theme(
    legend.position = "bottom"
  )

## ================================================================
## Figure 5: line plot of reject rate
## ================================================================
lambda_info <- plot_df_all %>%
  distinct(lambda_id, lambda_0, lambda_1, lambda_u) %>%
  arrange(lambda_id)

lambda_label_expr <- with(lambda_info, {
  labs <- mapply(
    function(l0, l1, lu) {
      bquote(
        atop(
          lambda[0] == .(l0) * "," ~ lambda[1] == .(l1),
          lambda[u] == .(round(lu, 3))
        )
      )
    },
    lambda_0, lambda_1, lambda_u,
    SIMPLIFY = FALSE
  )
  names(labs) <- as.character(lambda_id)
  labs
})

df_line <- plot_df_all

lambda_cols <- setNames(
  hue_pal()(length(lambda_info$lambda_id)),
  as.character(lambda_info$lambda_id)
)

fig5 <- ggplot(df_line,
               aes(x = epsilon, y = prop_reject,
                   color = lambda_id, group = lambda_id)) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = eps_grid) +
  scale_color_manual(
    values = lambda_cols,
    breaks = lambda_info$lambda_id,
    labels = lambda_label_expr,
    name   = "Loss parameter"
  ) +
  facet_grid(rows = vars(sample_size), cols = vars(scenario)) +
  labs(
    x = expression(paste("Privacy budget ", epsilon)),
    y = "Reject rate"
  ) +
  theme_fig_common

print(fig4)
print(fig5)

# ggsave("Fig4.png", fig4, width = 10, height = 6, dpi = 1000)
# ggsave("Fig5.png", fig5,  width = 10, height = 6, dpi = 1000)