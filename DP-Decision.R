# ==============================================================================
# DP-Decision: Decision-making under DP-FRT
# ------------------------------------------------------------------------------
# - Bayes Risk-optimal Decision Framework
# - Frequentist-calibrated Decision Framework
# ==============================================================================

if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
current_dir <- dirname(this.path::this.path())
source(file.path(current_dir, "DP-FRT.R"))

# Compute Psi from dp_frt() output
.dpfrt_posterior_prob_leq_alpha <- function(frt_fit, alpha){
  stopifnot(!is.null(frt_fit$posterior_p))
  u  <- frt_fit$posterior_p$u
  w  <- frt_fit$posterior_p$mass
  sum(w[u <= alpha])
}

# ------------------------------------------------------------------------------
# Frequentist calibration: Worst-case Calibration
# ------------------------------------------------------------------------------
.dpfrt_calibrate_freq_worst_case <- function(
    n1, n0, epsilon, alpha, alpha_freq,
    prior, method = c("enumerate","sampling"),
    n_sim = 2000, clip = NULL, seed = NULL
){
  method <- match.arg(method)
  n  <- n1 + n0
  if (!is.null(seed)) set.seed(seed)
  
  K_vals <- 0:n
  tK_vec <- numeric(length(K_vals))
  
  for (idx in seq_along(K_vals)){
    K <- K_vals[idx]
    if (K == 0 || K == n){
      n_sim_K <- max(20, ceiling(n_sim / 20))
    } else {
      n_sim_K <- n_sim
    }
    
    psi_sim <- numeric(n_sim_K)
    
    for (r in seq_len(n_sim_K)){
      A_r <- rhyper(1, m = K, n = n - K, k = n1)
      n11_r <- A_r
      n01_r <- K - A_r
      
      fit_r <- dp_frt(
        n1 = n1, n0 = n0,
        n11 = n11_r, n01 = n01_r,
        epsilon = epsilon,
        prior   = prior,
        method  = method,
        clip    = clip
      )
      
      psi_sim[r] <- .dpfrt_posterior_prob_leq_alpha(fit_r, alpha = alpha)
    }
    
    tK_vec[idx] <- as.numeric(stats::quantile(
      psi_sim,
      probs = 1 - alpha_freq,
      type  = 1
    ))
  }
  
  t_star_LFC <- max(tK_vec, na.rm = TRUE)
  
  list(
    K_vals     = K_vals,
    tK         = tK_vec,
    t_star     = t_star_LFC,
    calib_type = "worst_case"
  )
}

# ------------------------------------------------------------------------------
# Frequentist calibration: Data-adaptive Calibration with Confidence Sets
# ------------------------------------------------------------------------------
.dpfrt_calibrate_freq_data_adaptive <- function(
    n1, n0, epsilon, alpha, alpha_freq,
    prior, method = c("enumerate","sampling"),
    n_sim = 2000, zeta = 0.05,
    t11_obs, t01_obs,
    clip = NULL, seed = NULL
){
  method <- match.arg(method)
  n  <- n1 + n0
  if (!is.null(seed)) set.seed(seed)
  
  K_vals <- 0:n
  t1K_vec     <- rep(NA_real_, length(K_vals))
  in_conf_set <- rep(FALSE,    length(K_vals))
  
  alpha1 <- alpha_freq - zeta
  if (alpha1 <= 0) 
    stop("alpha_freq must be larger than zeta in data-adaptive calibration.")
  
  for (idx in seq_along(K_vals)){
    K <- K_vals[idx]
    n_sim_K <- n_sim
    
    psi_sim <- numeric(n_sim_K)
    tilde_mat <- matrix(NA_integer_, nrow = n_sim_K, ncol = 2)
    
    for (r in seq_len(n_sim_K)){
      A_r <- rhyper(1, m = K, n = n - K, k = n1)
      n11_r <- A_r
      n01_r <- K - A_r
      
      fit_r <- dp_frt(
        n1 = n1, n0 = n0,
        n11 = n11_r, n01 = n01_r,
        epsilon = epsilon,
        prior   = prior,
        method  = method,
        clip    = clip
      )
      
      psi_sim[r]      <- .dpfrt_posterior_prob_leq_alpha(fit_r, alpha = alpha)
      tilde_mat[r, ]  <- as.integer(unname(fit_r$meta$tilde[c("t11","t01")]))
    }
    
    tilde_df <- data.frame(t11 = tilde_mat[,1], t01 = tilde_mat[,2])
    freq_tab <- as.data.frame(table(tilde_df$t11, tilde_df$t01))
    names(freq_tab) <- c("t11","t01","freq")
    freq_tab$t11  <- as.integer(as.character(freq_tab$t11))
    freq_tab$t01  <- as.integer(as.character(freq_tab$t01))
    freq_tab$prob <- freq_tab$freq / sum(freq_tab$freq)
    
    freq_tab <- freq_tab[order(-freq_tab$prob), , drop = FALSE]
    freq_tab$cum_prob <- cumsum(freq_tab$prob)
    
    i_star <- which(freq_tab$cum_prob >= (1 - zeta))[1]
    
    idx_obs <- which(freq_tab$t11 == t11_obs & freq_tab$t01 == t01_obs)
    
    if (length(idx_obs) >= 1L){
      in_conf_set[idx] <- any(idx_obs <= i_star)
    } else {
      in_conf_set[idx] <- FALSE
    }
    
    t1K_vec[idx] <- as.numeric(stats::quantile(
      psi_sim,
      probs = 1 - alpha1,
      type  = 1
    ))
  }
  
  K_conf <- K_vals[in_conf_set]
  
  if (length(K_conf) == 0L){
    message("Data-adaptive calibration: C_{1-zeta}(tilde n) is empty based on simulations; fallback to all K.")
    K_conf <- K_vals
  }
  
  t_star_Neyman <- if (any(in_conf_set, na.rm = TRUE))
    max(t1K_vec[in_conf_set], na.rm = TRUE)
  else
    max(t1K_vec, na.rm = TRUE)
  
  
  list(
    K_vals      = K_vals,
    t1K         = t1K_vec,
    in_conf_set = in_conf_set,
    K_conf      = K_conf,
    t_star      = t_star_Neyman,
    calib_type  = "data_adaptive",
    zeta         = zeta
  )
}

# ================================================================
# Main function: decision-making
# ================================================================
dp_frt_decision <- function(
    frt_fit,
    alpha = 0.05,
    framework = c("bayes","freq"),
    # ---------- Bayesian decision parameters ----------
    bayes_abstain = TRUE,
    lambda_0 = 1,
    lambda_1 = 1,
    lambda_u = NULL,
    # ---------- Frequentist decision parameters ----------
    alpha_freq = 0.05,
    calib = c("worst_case","data_adaptive"),
    zeta   = 0.05,
    n_sim = 2000,
    seed  = NULL
){
  framework <- match.arg(framework)
  calib     <- match.arg(calib)
  
  meta <- frt_fit$meta
  n1   <- meta$n1
  n0   <- meta$n0
  n    <- n1 + n0
  epsilon <- meta$epsilon
  prior   <- meta$prior
  clip    <- meta$clip
  method  <- frt_fit$method
  
  psi_obs <- .dpfrt_posterior_prob_leq_alpha(frt_fit, alpha = alpha)
  
  # --------------------------------------------------------------
  # 1. Bayes Risk-optimal Decision Framework (Sec 4.1)
  # --------------------------------------------------------------
  if (framework == "bayes"){
    thresh_bayes <- lambda_0 / (lambda_0 + lambda_1)
    
    if (!bayes_abstain){
      delta_num <- as.integer(psi_obs > thresh_bayes)
      decision  <- if (delta_num == 1L) "reject" else "not_reject"
      
      return(list(
        framework   = "bayes",
        abstain     = FALSE,
        decision    = decision,
        delta       = delta_num,
        psi_obs     = psi_obs,
        threshold   = thresh_bayes,
        lambda_0    = lambda_0,
        lambda_1    = lambda_1,
        lambda_u    = NULL
      ))
    }
    
    H <- 2 * lambda_0 * lambda_1 / (lambda_0 + lambda_1)
    lambda_u <- lambda_u %||% (0.025 * H)
    
    if (lambda_u >= H / 2){
      delta_num <- as.integer(psi_obs > thresh_bayes)
      decision  <- if (delta_num == 1L) "reject" else "not_reject"
      
      return(list(
        framework   = "bayes",
        abstain     = FALSE,
        decision    = decision,
        delta       = delta_num,
        psi_obs     = psi_obs,
        threshold   = thresh_bayes,
        lambda_0    = lambda_0,
        lambda_1    = lambda_1,
        lambda_u    = lambda_u
      ))
    }
    
    t_center <- thresh_bayes
    t_low  <- min(t_center, lambda_u / lambda_1)
    t_high <- max(t_center, 1 - lambda_u / lambda_0)
    
    if (psi_obs > t_high){
      delta_num <- 1L
      decision  <- "reject"
    } else if (psi_obs < t_low){
      delta_num <- 0L
      decision  <- "not_reject"
    } else {
      delta_num <- NA_integer_
      decision  <- "abstain"
    }
    
    return(list(
      framework   = "bayes",
      abstain     = TRUE,
      alpha       = alpha,
      decision    = decision,
      delta       = delta_num,
      psi_obs     = psi_obs,
      threshold_center = t_center,
      threshold_low    = t_low,
      threshold_high   = t_high,
      lambda_0    = lambda_0,
      lambda_1    = lambda_1,
      lambda_u    = lambda_u
    ))
  }
  
  # --------------------------------------------------------------
  # 2. Frequentist-calibrated Decision Framework (Sec 4.2)
  # --------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  t11_obs <- as.integer(unname(meta$tilde["t11"]))
  t01_obs <- as.integer(unname(meta$tilde["t01"]))
  
  if (calib == "worst_case"){
    calib_res <- .dpfrt_calibrate_freq_worst_case(
      n1       = n1,
      n0       = n0,
      epsilon  = epsilon,
      alpha    = alpha,
      alpha_freq = alpha_freq,
      prior    = prior,
      method   = method,
      n_sim    = n_sim,
      clip     = clip,
      seed     = seed
    )
  } else { # data_adaptive
    calib_res <- .dpfrt_calibrate_freq_data_adaptive(
      n1       = n1,
      n0       = n0,
      epsilon  = epsilon,
      alpha    = alpha,
      alpha_freq = alpha_freq,
      prior    = prior,
      method   = method,
      n_sim    = n_sim,
      zeta      = zeta,
      t11_obs  = t11_obs,
      t01_obs  = t01_obs,
      clip     = clip,
      seed     = seed
    )
  }
  
  t_star <- calib_res$t_star
  delta_num <- as.integer(psi_obs > t_star)
  decision  <- if (delta_num == 1L) "reject" else "not_reject"
  
  return(list(
    framework   = "freq",
    calib       = calib_res$calib_type,
    decision    = decision,
    delta       = delta_num,
    psi_obs     = psi_obs,
    threshold   = t_star,
    alpha       = alpha,
    alpha_freq  = alpha_freq,
    zeta         = if (calib_res$calib_type == "data_adaptive") calib_res$zeta else NULL,
    meta        = list(
      n1      = n1,
      n0      = n0,
      n       = n,
      epsilon = epsilon,
      t11_obs = t11_obs,
      t01_obs = t01_obs
    ),
    calib_raw   = calib_res
  ))
}