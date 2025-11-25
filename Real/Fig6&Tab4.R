# ==============================================================================
## Figure 6 & Table 4: Real-data illustrations
# ==============================================================================

library(ggplot2)
library(patchwork)
library(matrixStats)

source("../DP-FRT.R")
source("../DP-Decision.R")

################################################################################
################################################################################
################################################################################
## Utilitys

# (1) non-private FRT p-value
compute_nonprivate_p <- function(n1, n0, n11, n01) {
  phyper(
    q = n11 - 1,
    m = n11 + n01,
    n = (n1 + n0) - (n11 + n01),
    k = n1,
    lower.tail = FALSE
  )
}

# (2) plotting function (posterior mass for p_FRT)
plot_post_bar <- function(res, p0, title = "") {
  df <- res$posterior_p
  ggplot(df, aes(x = u, y = mass)) +
    geom_col(width = 0.01, color = NA, alpha = 1, fill = "#A8C5FF") +
    geom_vline(xintercept = p0,
               linetype   = "dashed",
               color      = "#D7263D",
               linewidth  = 1) +
    labs(
      title = title,
      x     = expression(p[FRT]),
      y     = "Posterior mass"
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      panel.grid = element_blank()
    )
}

# (3) run DP-FRT + Bayes decision for a given dataset
run_dpfrt_for_dataset <- function(n1, n0, n11, n01,
                                  eps_list,
                                  seed       = 42,
                                  lambda_0   = 1,
                                  lambda_1   = 1,
                                  lambda_u   = 0.025,
                                  alpha_test = 0.05,
                                  prior      = list(type = "uniform"),
                                  method_use = "sampling",
                                  alpha_CI   = 0.05) {
  p_nonprivate <- compute_nonprivate_p(n1, n0, n11, n01)
  
  results   <- vector("list", length(eps_list))
  decisions <- vector("list", length(eps_list))
  plots     <- vector("list", length(eps_list))
  
  for (i in seq_along(eps_list)) {
    eps <- eps_list[i]
    
    res <- dp_frt(
      n1       = n1,
      n0       = n0,
      n11      = n11,
      n01      = n01,
      epsilon  = eps,
      prior    = prior,
      method   = method_use,
      alpha_CI = alpha_CI,
      seed     = seed,
      clip     = TRUE
    )
    
    res_bayes <- dp_frt_decision(
      frt_fit       = res,
      alpha         = alpha_test,
      framework     = "bayes",
      bayes_abstain = TRUE,
      lambda_0      = lambda_0,
      lambda_1      = lambda_1,
      lambda_u      = lambda_u
    )
    
    results[[i]]   <- res
    decisions[[i]] <- res_bayes
    plots[[i]]     <- plot_post_bar(
      res,
      p0    = p_nonprivate,
      title = bquote(epsilon == .(eps))
    )
  }
  
  list(
    p_nonprivate = p_nonprivate,
    results      = results,
    decisions    = decisions,
    plots        = plots
  )
}

## ================================================================
## Dataset 1
## ================================================================
n1  <- 7536
n0  <- 7540
n11 <- 569
n01 <- 590

eps_list <- c(0.1, 0.5, 1.0)

res_dataset1 <- run_dpfrt_for_dataset(
  n1       = n1,
  n0       = n0,
  n11      = n11,
  n01      = n01,
  eps_list = eps_list,
  seed     = 42
)

fig_real1 <- (res_dataset1$plots[[1]] |
                res_dataset1$plots[[2]] |
                res_dataset1$plots[[3]]) +
  plot_annotation(theme = theme())

print(fig_real1)
print(res_dataset1$p_nonprivate)
print(res_dataset1$decisions)

# ggsave(
#   "Fig6a.png",
#   plot   = fig_real1,
#   width  = 8,
#   height = 4,
#   dpi    = 1000,
#   units  = "in"
# )

## ================================================================
## Dataset 2
## ================================================================
n1  <- 7536
n0  <- 7540
n11 <- 44
n01 <- 53

res_dataset2 <- run_dpfrt_for_dataset(
  n1       = n1,
  n0       = n0,
  n11      = n11,
  n01      = n01,
  eps_list = eps_list,
  seed     = 42
)

fig_real2 <- (res_dataset2$plots[[1]] |
                res_dataset2$plots[[2]] |
                res_dataset2$plots[[3]]) +
  plot_annotation(theme = theme())

print(fig_real2)
print(res_dataset2$p_nonprivate)
print(res_dataset2$decisions)

# ggsave(
#   "Fig6b.png",
#   plot   = fig_real2,
#   width  = 8,
#   height = 4,
#   dpi    = 1000,
#   units  = "in"
# )