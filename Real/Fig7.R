# ==============================================================================
## Figure 7: Sequential DP-FRT with top-up release on Real Dataset
## ------------------------------------------------------------------------------
## Stage 1 uses a small initial budget eps_initial so that the Bayes decision
## lies in the abstention region. The analyst then obtains a lower bound on
## the required top-up budget eps_plus; Figure 7 shows stage-2 posterior CDFs
## at a sensitivity grid of eps_plus values around that lower bound, with
## line types distinguishing stage-2 decisions.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(matrixStats)

source("../DP-FRT.R")
source("../DP-Decision.R")
source("../DP-TopUp.R")

## ================================================================
## Dataset
## ================================================================
n1  <- 7536
n0  <- 7540
n11 <- 569
n01 <- 590

## ================================================================
## Stage-1 parameters and fit
## ================================================================
eps_initial  <- 0.06          # small on purpose -> stage-1 lies in abstain
xi_budget    <- 0.05
lambda_0_use <- 1
lambda_1_use <- 1
lambda_u_use <- 0.025
seed_topup   <- 343
alpha_test   <- 0.05

fit1 <- dp_frt(
  n1 = n1, n0 = n0, n11 = n11, n01 = n01,
  epsilon  = eps_initial,
  prior    = list(type = "uniform"),
  method   = "sampling",
  R        = 10000,
  alpha_CI = 0.05,
  seed     = seed_topup,
  clip     = TRUE
)

dec1 <- dp_frt_decision(
  frt_fit       = fit1,
  alpha         = alpha_test,
  framework     = "bayes",
  bayes_abstain = TRUE,
  lambda_0      = lambda_0_use,
  lambda_1      = lambda_1_use,
  lambda_u      = lambda_u_use
)

cat(sprintf("[Fig7] Stage 1: psi_obs = %.4f, decision = %s\n",
            dec1$psi_obs, dec1$decision))

## ================================================================
## Analyst-side lower bound for eps_plus
## ================================================================
lb_obj <- dp_frt_epsilon_plus_lower_bound(
  decision1 = dec1,
  frt_fit   = fit1,
  xi        = xi_budget,
  alpha     = alpha_test,
  seed      = seed_topup
)
eps_plus_lb <- lb_obj$epsilon_plus_lower_bound
cat(sprintf("[Fig7] analyst-side eps_plus lower bound = %.4f\n", eps_plus_lb))

## ================================================================
## Sensitivity grid: 4 values around the rounded lower bound
## ================================================================
lb_round      <- round(eps_plus_lb, 2)
eps_plus_grid <- lb_round + c(-0.04, -0.02, 0.00, 0.02)

stage2_list <- lapply(eps_plus_grid, function(ep) {
  top <- dp_frt_topup_release(
    n11 = n11, n01 = n01, n1 = n1, n0 = n0,
    epsilon_plus = ep,
    clip         = TRUE,
    seed         = seed_topup + 1000
  )
  res2 <- dp_frt_update_with_topup_and_decide(
    decision1    = dec1,
    tilde_plus   = top$tilde_plus,
    epsilon_plus = ep,
    frt_fit      = fit1,
    clip_plus    = TRUE,
    method       = "sampling",
    R            = 10000,
    alpha_CI     = 0.05,
    seed         = seed_topup + 2000
  )
  pp <- res2$frt_fit_plus$posterior_p
  cat(sprintf("  eps_plus = %.2f  psi2 = %.4f  decision = %s\n",
              ep, res2$decision2$psi_obs, res2$decision2$decision))
  data.frame(
    eps_plus = ep,
    u        = pp$u,
    cdf      = pp$cdf,
    decision = res2$decision2$decision
  )
})
stage2_df <- bind_rows(stage2_list)

## Non-private FRT p-value (dashed black reference)
p_nonpriv <- phyper(
  q = n11 - 1,
  m = n11 + n01,
  n = (n1 + n0) - (n11 + n01),
  k = n1,
  lower.tail = FALSE
)

## ================================================================
## Build the plotting data frame
## ================================================================
stage1_df <- data.frame(
  u   = fit1$posterior_p$u,
  cdf = fit1$posterior_p$cdf,
  grp = "Stage 1 (abstain)"
)

stage2_df$grp <- dplyr::case_when(
  abs(stage2_df$eps_plus - eps_plus_grid[1]) < 1e-6 ~
    sprintf("ε+ = %.2f (abstain)", eps_plus_grid[1]),
  abs(stage2_df$eps_plus - eps_plus_grid[2]) < 1e-6 ~
    sprintf("ε+ = %.2f (abstain)", eps_plus_grid[2]),
  abs(stage2_df$eps_plus - eps_plus_grid[3]) < 1e-6 ~
    sprintf("ε+ = %.2f (lower bound, not reject)", eps_plus_grid[3]),
  abs(stage2_df$eps_plus - eps_plus_grid[4]) < 1e-6 ~
    sprintf("ε+ = %.2f (not reject)", eps_plus_grid[4])
)

all_df <- bind_rows(
  stage1_df,
  stage2_df[, c("u", "cdf", "grp")]
)

grp_levels <- c(
  "Stage 1 (abstain)",
  sprintf("ε+ = %.2f (abstain)", eps_plus_grid[1]),
  sprintf("ε+ = %.2f (abstain)", eps_plus_grid[2]),
  sprintf("ε+ = %.2f (lower bound, not reject)", eps_plus_grid[3]),
  sprintf("ε+ = %.2f (not reject)", eps_plus_grid[4])
)
all_df$grp <- factor(all_df$grp, levels = grp_levels)

## ================================================================
## Plot: posterior CDFs at stage 2 under different top-up budgets
## ================================================================
col_vals <- setNames(
  c("grey40", "#F8766D", "#00BA38", "#000000", "#619CFF"),
  grp_levels
)
lty_vals <- setNames(
  c("dotted", "longdash", "longdash", "solid", "solid"),
  grp_levels
)
lw_vals <- setNames(
  c(1.3, 1.3, 1.3, 2.1, 1.3),
  grp_levels
)

fig7 <- ggplot(
  all_df,
  aes(x = u, y = cdf,
      color     = grp,
      linetype  = grp,
      linewidth = grp,
      group     = grp)
) +
  geom_line() +
  geom_vline(xintercept = alpha_test, linetype = "dotted",
             color = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = p_nonpriv, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  scale_color_manual(values = col_vals, name = NULL) +
  scale_linetype_manual(values = lty_vals, name = NULL, guide = "none") +
  scale_linewidth_manual(values = lw_vals, name = NULL, guide = "none") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(
    color = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(
        linetype  = unname(lty_vals),
        linewidth = unname(lw_vals)
      )
    )
  ) +
  labs(x = expression(p[FRT]), y = "Posterior CDF") +
  theme_bw(base_size = 18) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 15),
    legend.key.width = unit(40, "pt"),
    legend.spacing.x = unit(10, "pt"),
    panel.grid       = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x     = element_text(margin = margin(t = 18)),
    plot.title       = element_blank()
  )

print(fig7)

# ggsave(
#   "Fig7.png",
#   plot   = fig7,
#   width  = 11,
#   height = 6,
#   dpi    = 1000,
#   units  = "in"
# )
