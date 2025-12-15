# ==============================================================================
# Simulation settings for Cases 5--8 in the causal studies
# ==============================================================================

library(dplyr)
library(tidyr)

N_grid <- c(100, 500, 1000)

scen_defs_ca <- tibble::tibble(
  case_id = 5:8,
  case_label = c("Case 5: No effect",
                 "Case 6: Small effect",
                 "Case 7: Medium effect",
                 "Case 8: Large effect"),
  p_ctrl   = 0.50,
  p_treat  = c(0.50, 0.55, 0.60, 0.65)
) %>%
  mutate(tau = p_treat - p_ctrl)

science_grid <- tidyr::crossing(
  scen_defs_ca,
  N = N_grid
) %>%
  rowwise() %>%
  mutate(
    N01 = 0L,
    N11 = round(N * p_ctrl),
    N10 = round(N * tau),
    N00 = N - N11 - N10 - N01
  ) %>%
  ungroup() %>%
  select(case_label, N, N11, N10, N01, N00, tau)

science_grid_table <- science_grid %>%
  mutate(
    N11_table = N11 + N10,      # sum Y(1)
    N01_table = N11 + N01,      # sum Y(0)
    N10_table = N - N11_table,    # sum (1 - Y(1))
    N00_table = N - N01_table     # sum (1 - Y(0))
  ) %>%
  select(
    case_label, N,
    N11 = N11_table,
    N10 = N10_table,
    N01 = N01_table,
    N00 = N00_table,
    tau
  )

print(science_grid_table, n = nrow(science_grid_table))