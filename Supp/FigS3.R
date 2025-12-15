# ==============================================================================
# Figure S3: Prior Sensitivity (posterior mean MSE) for two prior families:
# (A) type = "common_rate" and (B) type = "beta_binom"
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(readr)
library(tibble)
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
sim_once_mean <- function(n11, n10, n01, n00, epsilon,
                          prior, alpha_CI = 0.05, method = "enumerate", clip_tilde = FALSE) {
  n1 <- n11 + n10; n0 <- n01 + n00
  p_true <- frt_p(a = n11, b = n01, n1 = n1, n0 = n0)
  out <- tryCatch({
    fit <- dp_frt(n1=n1, n0=n0, n11=n11, n01=n01, epsilon=epsilon,
                  prior=prior, method=method, alpha_CI=alpha_CI, clip=clip_tilde)
    c(p_true = p_true, p_mean = as.numeric(fit$summaries$mean))
  }, error = function(e) c(p_true=p_true, p_mean=NA_real_))
  as.numeric(out)
}

# (3) Replicate runs to compute MSE
sim_mse_mean <- function(n11, n10, n01, n00, epsilon, S, prior, method) {
  res <- vapply(seq_len(S), function(i)
    sim_once_mean(n11, n10, n01, n00, epsilon, prior, method=method),
    FUN.VALUE = c(p_true=0, p_mean=0))
  p_true <- res["p_true", ]
  est    <- res["p_mean", ]
  ok <- is.finite(est) & is.finite(p_true)
  if (!any(ok)) return(NA_real_)
  mean((est[ok] - p_true[ok])^2)
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

design_grid <- tidyr::crossing(scen_defs, sizes) %>%
  rowwise() %>%
  mutate(
    n1  = round(n_total * ratio / (1 + ratio)),
    n0  = n_total - n1,
    n11 = round(n1 * p_treat),
    n01 = round(n0 * p_ctrl),
    n10 = n1 - n11,
    n00 = n0 - n01
  ) %>%
  ungroup()

# -- two prior families to compare ---------------------------------------------
# A) common_rate: Beta(alpha, beta) shared by arms
priors_common <- tibble(
  prior_label = c("Beta(1,1)","Beta(0.5,0.5)","Beta(2,2)","Beta(1,2)","Beta(2,1)"),
  prior_obj   = list(
    list(type="common_rate", params=list(alpha=1,   beta=1)),
    list(type="common_rate", params=list(alpha=0.5, beta=0.5)),
    list(type="common_rate", params=list(alpha=2,   beta=2)),
    list(type="common_rate", params=list(alpha=1,   beta=2)),
    list(type="common_rate", params=list(alpha=2,   beta=1))
  )
)

# B) beta_binom: independent Beta priors for each arm (control & treatment)
priors_bb <- tibble(
  prior_label = c("Beta(1,1)","Beta(0.5,0.5)","Beta(2,2)","Beta(1,2)","Beta(2,1)"),
  prior_obj   = list(
    list(type="beta_binom", params=list(alpha=1,   beta=1)),
    list(type="beta_binom", params=list(alpha=0.5, beta=0.5)),
    list(type="beta_binom", params=list(alpha=2,   beta=2)),
    list(type="beta_binom", params=list(alpha=1,   beta=2)),
    list(type="beta_binom", params=list(alpha=2,   beta=1))
  )
)

################################################################################
################################################################################
################################################################################
## Simulation Settings

set.seed(2025)
eps_grid   <- c(0.2, 0.5, 1)
S_reps     <- 1000
method_for <- function(n1, n0) if ((n1 + n0) >= 300) "sampling" else "enumerate"

run_family <- function(design_grid, priors_tbl, eps_grid, S_reps, fam_label) {
  expanded <- tidyr::crossing(design_grid, epsilon = eps_grid, priors_tbl)
  res <- expanded %>%
    rowwise() %>%
    mutate(
      method = method_for(n1, n0),
      MSE_mean = sim_mse_mean(n11, n10, n01, n00,
                              epsilon = epsilon, S = S_reps,
                              prior  = prior_obj, method = method)
    ) %>%
    ungroup() %>%
    select(sample_size, scenario, n11, n10, n01, n00, epsilon,
           prior = prior_label, MSE_mean) %>%
    mutate(family = fam_label)
  res
}

res_common <- run_family(design_grid, priors_common, eps_grid, S_reps, fam_label = "common_rate")
res_bb     <- run_family(design_grid, priors_bb,     eps_grid, S_reps, fam_label = "beta_binom")

write_csv(res_bb,     "FigS3a.csv")
write_csv(res_common, "FigS3b.csv")

################################################################################
################################################################################
################################################################################
## Plotting

make_plot <- function(df, title) {
  plot_df <- df %>%
    mutate(
      sample_size = factor(sample_size, levels=c("small","medium","large"),
                           labels=c("n=100","n=500", "n=1000")),
      scenario = factor(scenario, levels=c("no_effect","small","medium","large"),
                        labels=c("Case 1: No effect",
                                 "Case 2: Small effect",
                                 "Case 3: Medium effect",
                                 "Case 4: Large effect"))
    )
  
  ggplot(plot_df, aes(x = epsilon, y = MSE_mean, color = prior, group = prior)) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = eps_grid) +
    facet_grid(rows = vars(sample_size), cols = vars(scenario), scales = "free_y") +
    labs(
      x = expression(paste("Privacy budget ", epsilon)),
      y = "MSE of posterior mean",
      color = "Prior"
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
}

fig_bb     <- make_plot(res_bb)
fig_common <- make_plot(res_common)

print(fig_bb)
print(fig_common)

# ggsave("FigS3a.png", fig_bb,     width = 10, height = 6, dpi = 300)
# ggsave("FigS3b.png", fig_common, width = 10, height = 6, dpi = 300)