# ==============================================================================
# Figure 2: Illustration of the Bayes risk-optimal decision
# ==============================================================================

suppressWarnings({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("tidyr",    quietly = TRUE)) install.packages("tidyr")
})
library(ggplot2)
library(tidyr)
source("../DP-FRT.R")

decision_plot <- function(lambda0 = 1, lambda1 = 1,
                          b = NULL, conf = NULL, lambda_u = NULL,
                          show_labels = TRUE, clamp_b = TRUE,
                          base_size = 16,
                          palette = c(r0 = "#0072B2",
                                      r1 = "#E69F00",
                                      ru = "#009E73"),
                          legend_pos = c(0.85, 0.5)   
) {
  stopifnot(lambda0 > 0, lambda1 > 0)
  H  <- (lambda0 * lambda1) / (lambda0 + lambda1)
  w0 <- lambda0 / (lambda0 + lambda1)
  w1 <- 1 - w0
  
  if (!is.null(lambda_u)) {
    stopifnot(lambda_u >= 0)
    b_calc <- 1 - lambda_u / H
  } else if (!is.null(conf)) {
    stopifnot(conf > 0, conf < 1)
    b_calc <- conf
  } else if (!is.null(b)) {
    stopifnot(b >= 0, b <= 1)
    b_calc <- b
  } else {
    b_calc <- 0.60
  }
  
  if (clamp_b) {
    if (b_calc < 0 || b_calc > 1) {
      warning(sprintf("b=%.3f beyond [0,1], already being truncated.", b_calc))
      b_calc <- min(1, max(0, b_calc))
    }
  }
  
  lambda_u_eff <- H * (1 - b_calc)
  kappa        <- w0 * w1 * (1 - b_calc)
  
  t0 <- w0 * (1 - b_calc)
  t1 <- w0 + w1 * b_calc
  
  q  <- seq(0, 1, length.out = 1201)
  df <- data.frame(
    q  = q,
    r0 = w1 * q,              
    r1 = w0 * (1 - q),         
    ru = rep(kappa, length(q))
  )
  dfl <- pivot_longer(df, cols = c("r0", "r1", "ru"),
                      names_to = "risk", values_to = "value")
  dfl$risk <- factor(dfl$risk, levels = c("r0", "r1", "ru"))
  
  y_max <- max(dfl$value) * 1.05
  
  p <- ggplot() +
    geom_rect(aes(xmin = t0, xmax = t1, ymin = -Inf, ymax = Inf),
              fill = "grey85", alpha = 0.35, inherit.aes = FALSE) +
    geom_line(data = dfl,
              aes(x = q, y = value, color = risk, linetype = risk),
              linewidth = 1.3) +
    scale_color_manual(
      name = NULL,
      values = palette,
      labels = c(
        r0 = expression(R(delta==0 ~ "|" ~
                            tilde(bold(n)))),
        r1 = expression(R(delta==1 ~ "|" ~
                            tilde(bold(n)))),
        ru = expression(R(delta==u ~ "|" ~
                            tilde(bold(n))))
      )
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(r0 = "solid", r1 = "solid", ru = "twodash"),
      labels = c("", "", "")
    ) +
    geom_vline(xintercept = c(t0, t1), linetype = "dotted", linewidth = 0.9) +
    { if (show_labels)
      annotate("label", x = (0 + t0)/2,  y = y_max*0.98, label = "Decision: 0",
               size = 4, alpha = 0.9, hjust = 0.5) } +
    { if (show_labels && t1 > t0)
      annotate("label", x = (t0 + t1)/2, y = y_max*0.98, label = "Decision: u",
               size = 4, alpha = 0.9, hjust = 0.5) } +
    { if (show_labels)
      annotate("label", x = (t1 + 1)/2,  y = y_max*0.98, label = "Decision: 1",
               size = 4, alpha = 0.9, hjust = 0.5) } +
    labs(
      x = expression(Pr(p[FRT] <= alpha ~ "|" ~ tilde(bold(n)))),
      y = "Normalized posterior risk",
      caption = NULL
    ) +
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(linetype = c("solid","solid","twodash"), linewidth = 1.3)
      ),
      linetype = "none"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold"),
      legend.position   = legend_pos,
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.85),
                                       color = NA),
      legend.key.height = unit(12, "pt"),
      legend.key.width  = unit(24, "pt"),
      legend.text       = element_text(size = base_size * 0.9),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x= element_blank(),
      panel.grid = element_blank(),
      axis.title.x      = element_text(margin = margin(t = 8)),
      axis.title.y      = element_text(margin = margin(r = 8))
    )
  
  print(p)
  invisible(list(
    w0 = w0, w1 = w1, b = b_calc, t0 = t0, t1 = t1,
    kappa = kappa, lambda_u = lambda_u_eff, H = H
  ))
}

decision_plot_by_conf <- function(lambda0 = 1, lambda1 = 1, conf = 0.95, ...) {
  stopifnot(conf > 0, conf < 1)
  decision_plot(lambda0 = lambda0, lambda1 = lambda1, conf = conf, ...)
}

decision_plot(lambda0 = 0.2, lambda1 = 0.5, lambda_u = 0.025*0.2857,
              legend_pos = c(0.9, 0.53))
res <- decision_plot(lambda0 = 0.2, lambda1 = 0.5, lambda_u = 0.1,
                     legend_pos = c(0.9, 0.53))

# ggsave("Fig2.png",
#        plot = last_plot(),
#        width = 8,
#        height = 5,
#        dpi = 1000,
#        units = "in")