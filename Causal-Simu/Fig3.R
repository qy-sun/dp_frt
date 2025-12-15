# ==============================================================================
# Figure 3: Average posterior cumulative distribution function for FRT p-value
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source("../DP-FRT.R")
source("../DP-Decision.R")

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

eps_grid   <- c(0.2, 0.5, 1)
alpha_use  <- 0.05
prior_cr   <- list(type = "uniform")
method_use <- "enumerate"
clip_tilde <- TRUE

S_reps <- 100
set.seed(42)

cdf_all_raw <- purrr::map_dfr(seq_len(nrow(design_grid)), function(i) {
  scen_row <- design_grid[i, ]
  n11 <- scen_row$n11
  n10 <- scen_row$n10
  n01 <- scen_row$n01
  n00 <- scen_row$n00
  n1  <- n11 + n10
  n0  <- n01 + n00
  
  purrr::map_dfr(eps_grid, function(ep) {
    purrr::map_dfr(seq_len(S_reps), function(s_id) {
      fit <- dp_frt(
        n1 = n1, n0 = n0,
        n11 = n11, n01 = n01,
        epsilon = ep,
        prior   = prior_cr,
        method  = method_use,
        alpha   = alpha_use,
        clip    = clip_tilde
      )
      post <- fit$posterior_p
      psi  <- .dpfrt_posterior_prob_leq_alpha(fit, alpha_use)
      
      tibble(
        sim_id     = s_id,
        scenario   = scen_row$scenario,
        sample_size = scen_row$sample_size,
        epsilon    = ep,
        u          = post$u,
        cdf        = post$cdf,
        psi        = psi
      )
    })
  })
})

cdf_all_mean <- cdf_all_raw |>
  group_by(scenario, sample_size, epsilon, u) |>
  summarise(
    cdf = mean(cdf),
    psi = mean(psi),
    .groups = "drop"
  )

cdf_all_mean <- cdf_all_mean |>
  mutate(
    scenario = factor(
      scenario,
      levels = c("no_effect","small","medium","large"),
      labels = c("Case 1: No effect",
                 "Case 2: Small effect",
                 "Case 3: Medium effect",
                 "Case 4: Large effect")
    ),
    sample_size = factor(
      sample_size,
      levels = c("small","medium","large"),
      labels = c("n=100","n=500", "n=1000")
    ),
    epsilon_f = factor(
      epsilon,
      levels = eps_grid,
      labels = eps_grid
    )
  )

psi_points <- cdf_all_mean |>
  group_by(scenario, sample_size, epsilon_f, epsilon) |>
  summarise(
    psi = mean(psi),
    .groups = "drop"
  ) |>
  mutate(x_alpha = alpha_use)



cols <- scales::hue_pal()(length(eps_grid))
names(cols) <- levels(cdf_all_mean$epsilon_f)

fig3_mean <- ggplot(cdf_all_mean,
                    aes(x = u, y = cdf,
                        color = epsilon_f, group = epsilon_f)) +
  geom_step(linewidth = 0.75) +
  geom_point(data = psi_points,
             aes(x = x_alpha, y = psi, color = epsilon_f),
             size = 1.5,
             inherit.aes = FALSE) +
  geom_vline(xintercept = alpha_use, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = cols,
    labels = parse(text = paste0("epsilon == ", eps_grid)),
    name   = "Privacy Budget"
  ) +
  facet_grid(rows = vars(sample_size), cols = vars(scenario)) +
  labs(
    x = expression(p[FRT]),
    y = "Posterior CDF",
    color = expression(epsilon)
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size = 13),
    strip.background  = element_rect(fill = "gray90"),
    strip.text.x      = element_text(size = 13, face = "bold"),
    strip.text.y      = element_text(size = 13, face = "bold"),
    panel.grid        = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.title.x      = element_text(margin = margin(t = 15)),
    plot.title        = element_blank()
  )

print(fig3_mean)
# ggsave("Fig3.png", fig3_mean, width = 12, height = 6, dpi = 1000)