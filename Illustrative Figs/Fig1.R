# ==============================================================================
# Figure 1: Posterior distributions of Fisher's randomization p-value
# ==============================================================================

library(ggplot2)
library(patchwork)
library(matrixStats)
source("../DP-FRT.R")

n1  <- 500
n0  <- 500
n11 <- 260
n01 <- 250

p_nonprivate <- phyper(q = n11 - 1, 
                       m = n11 + n01, 
                       n = (n1 + n0) - (n11 + n01), 
                       k = n1, 
                       lower.tail = FALSE)

plot_post_bar <- function(res, p0, title = ""){
  df <- res$posterior_p
  ggplot(df, aes(x = u, y = mass)) +
    geom_col(width = 0.01, color = NA, alpha = 1, fill = "#A8C5FF") +
    geom_vline(xintercept = p0, linetype = "dashed", color = "#D7263D", linewidth = 1) +
    labs(
      title = title,
      x = expression(p[FRT]),
      y = "Posterior mass"
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      panel.grid = element_blank()
    )
}

eps_list <- c(0.1, 0.5, 1.0)
seed = 42

results <- vector("list", length(eps_list))
plots   <- vector("list", length(eps_list))

for (i in seq_along(eps_list)){
  eps <- eps_list[i]
  s   <- seed
  
  res <- dp_frt(
    n1 = n1, n0 = n0, n11 = n11, n01 = n01,
    epsilon = eps,
    prior   = list(type = "uniform"),
    method  = "enumerate",
    alpha   = 0.05,
    seed    = s,
    clip    = TRUE
  )
  
  results[[i]] <- res
  plots[[i]] <- plot_post_bar(res, p_nonprivate,
                              title = bquote(epsilon == .(eps)))
}

final_plot <- (plots[[1]] | plots[[2]] | plots[[3]]) +
  plot_annotation(theme = theme())

print(final_plot)

# ggsave("Fig1.png",
#        plot = last_plot(),  
#        width = 8,
#        height = 4,
#        dpi = 1000,
#        units = "in")