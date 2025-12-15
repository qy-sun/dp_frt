# ==============================================================================
# Figure S2: MSE of the Bayesian denoising estimators
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
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
sim_once <- function(n11, n10, n01, n00, epsilon,
                     prior = list(type="uniform"),
                     alpha_CI = 0.05, method = "enumerate", clip_tilde = FALSE) {
  n1 <- n11 + n10
  n0 <- n01 + n00
  p_true <- frt_p(a = n11, b = n01, n1 = n1, n0 = n0)
  
  out <- tryCatch({
    fit <- dp_frt(
      n1 = n1, n0 = n0, n11 = n11, n01 = n01, epsilon = epsilon,
      prior = prior, method = method, alpha_CI = alpha_CI, clip = clip_tilde
    )
    t11 <- as.integer(round(fit$meta$tilde["t11"]))
    t01 <- as.integer(round(fit$meta$tilde["t01"]))
    t11 <- max(0L, min(n1, t11))
    t01 <- max(0L, min(n0, t01))
    
    c(
      p_true   = p_true,
      p_mean   = as.numeric(fit$summaries$mean),
      p_median = as.numeric(fit$summaries$median),
      p_map    = as.numeric(fit$summaries$MAP)
    )
  }, error = function(e) {
    c(p_true = p_true, p_mean = NA_real_,
      p_median = NA_real_, p_map = NA_real_)
  })
  
  out <- as.numeric(out)
  names(out) <- c("p_true","p_mean","p_median","p_map")
  out
}

# (3) Replicate runs to compute MSE
sim_mse <- function(n11, n10, n01, n00, epsilon, S = 2000,
                    prior = list(type="uniform"),
                    method = "enumerate") {
  res <- vapply(
    X = seq_len(S),
    FUN = function(i) sim_once(n11, n10, n01, n00, epsilon,
                               prior = prior, method = method),
    FUN.VALUE = c(p_true=0, p_mean=0, p_median=0, p_map=0)
  )
  p_true <- res["p_true", ]
  mse <- function(est) {
    ok <- is.finite(est) & is.finite(p_true)
    mean((est[ok] - p_true[ok])^2)
  }
  c(
    MSE_mean   = mse(res["p_mean", ]),
    MSE_median = mse(res["p_median", ]),
    MSE_map    = mse(res["p_map", ])
  )
}

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
################################################################################
################################################################################
## Simulation Settings

set.seed(2025)
eps_grid   <- c(0.2, 0.5, 1)
S_reps     <- 1000
prior_cr   <- list(type="uniform")
method_use <- "sampling"

mse_df_all <- imap_dfr(grouped_scenarios, function(group_list, size_label) {
  imap_dfr(group_list, function(tab, sc_name) {
    n11 <- tab["n11"]; n10 <- tab["n10"]; n01 <- tab["n01"]; n00 <- tab["n00"]
    map_dfr(eps_grid, function(ep) {
      mses <- sim_mse(n11, n10, n01, n00, epsilon = ep, S = S_reps,
                      prior = prior_cr, method = method_use)
      tibble(
        sample_size = size_label,
        scenario = sc_name,
        n11 = n11, n10 = n10, n01 = n01, n00 = n00,
        epsilon = ep,
        MSE_mean   = mses["MSE_mean"],
        MSE_median = mses["MSE_median"],
        MSE_map    = mses["MSE_map"]
      )
    })
  })
})

plot_df_all <- mse_df_all |>
  pivot_longer(cols = starts_with("MSE_"),
               names_to = "method", values_to = "MSE") |>
  mutate(method = factor(method,
                         levels = c("MSE_mean","MSE_median","MSE_map"),
                         labels = c("Posterior Mean","Posterior Median","Posterior MAP")),
         sample_size = factor(sample_size, levels = c("small","medium","large"),
                              labels = c("n=100","n=500", "n=1000")),
         scenario = factor(scenario, levels = c("no_effect","small","medium","large"),
                           labels = c("Case 1: No effect",
                                      "Case 2: Small effect",
                                      "Case 3: Medium effect",
                                      "Case 4: Large effect")))

figS2_grid <- ggplot(plot_df_all, aes(x = epsilon, y = MSE, color = method, group = method)) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = eps_grid) +
  scale_y_continuous() +
  scale_color_manual(
    values = c(
      "Posterior Mean" = scales::hue_pal()(3)[1],
      "Posterior Median" = scales::hue_pal()(3)[2],
      "Posterior MAP" = scales::hue_pal()(3)[3]
    )
  ) +
  facet_grid(rows = vars(sample_size), cols = vars(scenario), scales = "free_y") +
  labs(
    x = expression(paste("Privacy budget ", epsilon)),
    y = "Mean Square Error (MSE)",
    color = "Estimator"
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

print(figS2_grid)

# ggsave("FigS2.png", figS2_grid, width = 10, height = 6, dpi = 1000)