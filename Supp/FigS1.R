# ==============================================================================
# Figure S1: MSE comparison of DP-FRT-p vs DP-FRT-t vs DP-FRT-Bayes (Posterior Mean)
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Make sure DP-FRT.R defines: %||%, dp_frt, dp_frt_p, dp_frt_t
source("../DP-FRT.R")

################################################################################
## Utilities

# (1) Non-private FRT right-tail p under sharp null given (a,b) and margins
frt_p <- function(a, b, n1, n0) {
  K <- a + b
  phyper(q = a - 1, m = K, n = (n1 + n0) - K, k = n1, lower.tail = FALSE)
}

# (2) One Monte Carlo run: DP-FRT-Bayes (posterior mean) + DP-FRT-p + DP-FRT-t
sim_once_s1 <- function(n11, n10, n01, n00, epsilon,
                        prior = list(type="uniform"),
                        alpha_CI = 0.05,
                        method = "enumerate",
                        clip_tilde = FALSE) {
  n1 <- n11 + n10
  n0 <- n01 + n00
  p_true <- frt_p(a = n11, b = n01, n1 = n1, n0 = n0)
  
  # --- DP-FRT-Bayes: posterior mean ---
  out_bayes_mean <- tryCatch({
    fit <- dp_frt(
      n1 = n1, n0 = n0, n11 = n11, n01 = n01, epsilon = epsilon,
      prior = prior, method = method, alpha_CI = alpha_CI, clip = clip_tilde
    )
    as.numeric(fit$summaries$mean)
  }, error = function(e) NA_real_)
  
  # --- DP-FRT-p ---
  out_p <- tryCatch({
    fitp <- dp_frt_p(n1 = n1, n0 = n0, n11 = n11, n01 = n01, epsilon = epsilon)
    as.numeric(fitp$private$p_tilde)
  }, error = function(e) NA_real_)
  
  # --- DP-FRT-t ---
  out_t <- tryCatch({
    fitt <- dp_frt_t(n1 = n1, n0 = n0, n11 = n11, n01 = n01, epsilon = epsilon)
    as.numeric(fitt$private$p_tilde)
  }, error = function(e) NA_real_)
  
  c(
    p_true      = p_true,
    p_bayesMean = out_bayes_mean,
    p_dpftr_p   = out_p,
    p_dpftr_t   = out_t
  )
}

# (3) Replicate runs to compute MSE for three methods
sim_mse_s1 <- function(n11, n10, n01, n00, epsilon, S = 2000,
                       prior = list(type="uniform"),
                       method = "enumerate") {
  res <- vapply(
    X = seq_len(S),
    FUN = function(i) sim_once_s1(n11, n10, n01, n00, epsilon,
                                  prior = prior, method = method),
    FUN.VALUE = c(p_true=0, p_bayesMean=0, p_dpftr_p=0, p_dpftr_t=0)
  )
  p_true <- res["p_true", ]
  mse <- function(est) {
    ok <- is.finite(est) & is.finite(p_true)
    mean((est[ok] - p_true[ok])^2)
  }
  c(
    MSE_BayesMean = mse(res["p_bayesMean", ]),
    MSE_DPFRT_p   = mse(res["p_dpftr_p", ]),
    MSE_DPFRT_t   = mse(res["p_dpftr_t", ])
  )
}

################################################################################
## Scenario definitions (deterministic tables)

scen_defs <- tibble::tibble(
  scenario = factor(c("no_effect","small","medium","large"),
                    levels = c("no_effect","small","medium","large")),
  p_ctrl   = c(0.50, 0.50, 0.50, 0.50),
  p_treat  = c(0.50, 0.55, 0.65, 0.80),
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
## Simulation settings

set.seed(2025)
eps_grid <- c(0.2, 0.5, 1)
S_reps   <- 10
prior_cr <- list(type="uniform")
method_use <- "sampling"

################################################################################
## Run simulations

mse_df_s1 <- imap_dfr(grouped_scenarios, function(group_list, size_label) {
  imap_dfr(group_list, function(tab, sc_name) {
    n11 <- tab["n11"]; n10 <- tab["n10"]; n01 <- tab["n01"]; n00 <- tab["n00"]
    map_dfr(eps_grid, function(ep) {
      mses <- sim_mse_s1(n11, n10, n01, n00, epsilon = ep, S = S_reps,
                         prior = prior_cr, method = method_use)
      tibble(
        sample_size = size_label,
        scenario    = sc_name,
        n11 = n11, n10 = n10, n01 = n01, n00 = n00,
        epsilon = ep,
        MSE_BayesMean = mses["MSE_BayesMean"],
        MSE_DPFRT_p   = mses["MSE_DPFRT_p"],
        MSE_DPFRT_t   = mses["MSE_DPFRT_t"]
      )
    })
  })
})

################################################################################
## Plot (Figure S1)

plot_df_s1 <- mse_df_s1 |>
  pivot_longer(cols = starts_with("MSE_"),
               names_to = "method", values_to = "MSE") |>
  mutate(
    method = factor(
      method,
      levels = c("MSE_BayesMean","MSE_DPFRT_p","MSE_DPFRT_t"),
      labels = c("DP-FRT-Bayes (Posterior Mean)","DP-FRT-p","DP-FRT-t")
    ),
    sample_size = factor(sample_size, levels = c("small","medium","large"),
                         labels = c("n=100","n=500","n=1000")),
    scenario = factor(scenario, levels = c("no_effect","small","medium","large"),
                      labels = c("Case 1: No effect",
                                 "Case 2: Small effect",
                                 "Case 3: Medium effect",
                                 "Case 4: Large effect"))
  )

figS1_grid <- ggplot(plot_df_s1, 
                     aes(x = epsilon, y = MSE, 
                         color = method, shape = method, group = method)) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = eps_grid) +
  facet_grid(rows = vars(sample_size), cols = vars(scenario), scales = "free_y") +
  labs(
    x = expression(paste("Privacy budget ", epsilon)),
    y = "Mean Square Error (MSE)",
    color = "Method"
  ) +
  scale_shape_manual(
    values = c(
      "DP-FRT-Bayes (Posterior Mean)" = 16,
      "DP-FRT-p"                      = 17,
      "DP-FRT-t"                      = 15
    ),
    name = "Method"
  ) +
  scale_color_manual(
    values = c(
      "DP-FRT-Bayes (Posterior Mean)" = "#D55E00",
      "DP-FRT-p"                      = "#009E73",
      "DP-FRT-t"                      = "#0072B2"
    ),
    name = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text.x = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_text(margin = margin(t = 12))
  )

print(figS1_grid)

# ggsave("FigS1.png", figS1_grid, width = 10, height = 6, dpi = 1000)