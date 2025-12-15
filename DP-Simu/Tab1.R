# ==============================================================================
# Table 1: Simulation settings for Cases 1--4 in the DP studies.
# ==============================================================================

library(dplyr)
library(tidyr)

frt_p <- function(a, b, n1, n0) {
  K <- a + b
  phyper(q = a - 1, m = K, n = (n1 + n0) - K, k = n1, lower.tail = FALSE)
}

scen_defs <- data.frame(
  case = paste0("Case ", 1:4),
  label = c("No effect", "Small effect", "Medium effect", "Large effect"),
  p_ctrl = 0.50,
  p_treat = c(0.50, 0.55, 0.60, 0.65),
  ratio = 1
)

sizes <- data.frame(
  sample_size = factor(c("n=100", "n=500", "n=1000"),
                       levels = c("n=100", "n=500", "n=1000"),
                       ordered = TRUE),
  n_total = c(100, 500, 1000)
)

design_grid <- tidyr::crossing(scen_defs, sizes) |>
  rowwise() |>
  mutate(
    n1  = round(n_total * ratio / (1 + ratio)),
    n0  = n_total - n1,
    n11 = round(n1 * p_treat),
    n01 = round(n0 * p_ctrl),
    n10 = n1 - n11,
    n00 = n0 - n01,
    p_nonprivate = frt_p(a = n11, b = n01, n1 = n1, n0 = n0)
  ) |>
  ungroup() |>
  select(case, label, sample_size, n_total, n11, n10, n01, n00, p_nonprivate)

print(design_grid, n = nrow(design_grid))