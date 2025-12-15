# ==============================================================================
# DP-FRT-Bayes: Mechanism-aware Bayesian Denoising Method
# ==============================================================================

# Utility: null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Main Function
dp_frt <- function(
    n1, n0, n11, n01, epsilon,
    prior = list(
      type   = c("uniform", "beta_binom", "common_rate")[1],
      params = list()
    ),
    method  = c("sampling", "enumerate")[1],
    R       = 10000,        # Sample size for "sampling" (ignored for enumeration)
    alpha_CI   = 0.05,      # Confidence level 1 - alpha_CI
    seed    = NULL,         # Random seed (used in sampling)
    clip    = FALSE         # Whether to clip noisy counts (posterior unaffected per lemma)
){
  # ---------- Input checks ----------
  stopifnot(n11 <= n1, n01 <= n0, n1 > 0, n0 > 0)
  stopifnot(is.finite(epsilon), epsilon > 0)
  n  <- n1 + n0
  rho <- exp(-epsilon)
  log_rho <- -epsilon
  log_const <- log1p(-rho) - log1p(rho)
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---------- Two-sided geometric mechanism ----------
  rgeom_two_sided <- function(size, rho) {
    stopifnot(rho > 0, rho < 1)
    rgeom(size, prob = 1 - rho) - rgeom(size, prob = 1 - rho)
  }
  
  eta11 <- rgeom_two_sided(1, rho)
  eta01 <- rgeom_two_sided(1, rho)
  t11   <- n11 + eta11
  t01   <- n01 + eta01
  
  if (clip){
    t11 <- max(0L, min(n1, t11))
    t01 <- max(0L, min(n0, t01))
  }
  
  log_kappa <- function(h_abs) log_const + h_abs * log_rho
  
  # ---------- Prior log weight for (a,b) ----------
  # prior$type ∈ {"uniform", "beta_binom", "common_rate"}
  # "uniform": no parameters
  # "beta_binom": params = list(alpha1, beta1, alpha0, beta0)
  # "common_rate": params = list(alpha, beta)
  get_logprior <- function(a, b){
    type <- prior$type
    par  <- prior$params %||% list()
    if (type == "uniform"){
      return(0)
    } else if (type == "beta_binom"){
      alpha1 <- par$alpha1 %||% 1
      beta1  <- par$beta1  %||% 1
      alpha0 <- par$alpha0 %||% 1
      beta0  <- par$beta0  %||% 1
      lp1 <- lchoose(n1, a) + lbeta(a + alpha1, n1 - a + beta1) - lbeta(alpha1, beta1)
      lp0 <- lchoose(n0, b) + lbeta(b + alpha0, n0 - b + beta0) - lbeta(alpha0, beta0)
      return(lp1 + lp0)
    } else if (type == "common_rate"){
      alpha <- par$alpha %||% 1
      beta  <- par$beta  %||% 1
      K     <- a + b
      return(lchoose(n1, a) + lchoose(n0, b) +
               lbeta(K + alpha, n - K + beta) - lbeta(alpha, beta))
    } else {
      stop("Unknown prior$type. Use 'uniform', 'beta_binom', or 'common_rate'.")
    }
  }
  
  # ---------- p_FRT(a,b): hypergeometric right tail ----------
  p_frta <- function(a, b){
    K <- a + b
    phyper(q = a - 1, m = K, n = n - K, k = n1, lower.tail = FALSE)
  }
  
  # ---------- Compute posterior (or sample) ----------
  if (method == "enumerate"){
    A <- 0:n1
    B <- 0:n0
    
    logw_mat <- matrix(-Inf, nrow = n1 + 1, ncol = n0 + 1)
    p_mat    <- matrix(NA_real_, nrow = n1 + 1, ncol = n0 + 1)
    
    for (a in A){
      for (b in B){
        lp <- get_logprior(a,b) +
          log_kappa(abs(t11 - a)) + log_kappa(abs(t01 - b))
        logw_mat[a + 1, b + 1] <- lp
        p_mat[a + 1, b + 1]    <- p_frta(a,b)
      }
    }
    
    lse <- matrixStats::logSumExp(as.numeric(logw_mat))
    w_mat <- exp(logw_mat - lse)
    
    grid <- expand.grid(a = A, b = B)
    grid$p      <- as.vector(p_mat)
    grid$weight <- as.vector(w_mat)
    
    p_key <- round(grid$p, 12)
    agg <- aggregate(grid$weight, by = list(u = p_key), FUN = sum)
    names(agg) <- c("u","mass")
    agg <- agg[order(agg$u), , drop = FALSE]
    agg$cdf <- cumsum(agg$mass)
    
    post_mean   <- sum(agg$u * agg$mass)
    idx_med <- which(agg$cdf >= 0.5)
    med_u   <- agg$u[min(idx_med)]
    idx_map <- which.max(agg$mass)
    map_u   <- agg$u[idx_map]
    
    L <- agg$u[which(agg$cdf >= alpha_CI/2)[1]]
    U <- agg$u[which(agg$cdf >= 1 - alpha_CI/2)[1]]
    cred_set <- agg$u[agg$u >= L & agg$u <= U]
    
    ord <- order(agg$mass, decreasing = TRUE)
    mass_cum <- cumsum(agg$mass[ord])
    k <- which(mass_cum >= (1 - alpha_CI))[1]
    hpd_set <- sort(agg$u[ord][seq_len(k)])
    
    tau_vec  <- grid$a / n1 - grid$b / n0
    
    a <- grid$a; b <- grid$b
    c_ <- n1 - a; d_ <- n0 - b
    
    need_cc_rr <- (a == 0L) | (b == 0L)
    risk_t <- ifelse(need_cc_rr, (a + 0.5) / (n1 + 1), a / n1)
    risk_c <- ifelse(need_cc_rr, (b + 0.5) / (n0 + 1), b / n0)
    rr_vec <- risk_t / risk_c
    
    need_cc_or <- (a == 0L) | (b == 0L) | (c_ == 0L) | (d_ == 0L)
    a_cc <- ifelse(need_cc_or, a  + 0.5, a)
    b_cc <- ifelse(need_cc_or, b  + 0.5, b)
    c_cc <- ifelse(need_cc_or, c_ + 0.5, c_)
    d_cc <- ifelse(need_cc_or, d_ + 0.5, d_)
    or_vec <- (a_cc * d_cc) / (c_cc * b_cc)
    
    w <- grid$weight
    summarize_weighted <- function(x, w, alpha_CI){
      mu <- sum(x * w)
      o  <- order(x)
      xs <- x[o]; ws <- w[o]; cws <- cumsum(ws)
      qL <- xs[which(cws >= alpha_CI/2)[1]]
      qU <- xs[which(cws >= 1 - alpha_CI/2)[1]]
      list(mean = mu, qL = qL, qU = qU)
    }
    tau_sum <- summarize_weighted(tau_vec, w, alpha_CI)
    rr_sum  <- summarize_weighted(rr_vec,  w, alpha_CI)
    or_sum  <- summarize_weighted(or_vec,  w, alpha_CI)
    
    return(list(
      method   = "enumerate",
      meta     = list(n1 = n1, n0 = n0, n11 = n11, n01 = n01,
                      clip = clip, epsilon = epsilon, rho = rho,
                      tilde = c(t11 = t11, t01 = t01),
                      prior  = prior, alpha_CI = alpha_CI),
      posterior_p = agg[, c("u","mass","cdf")],
      summaries = list(
        mean   = post_mean,
        median = med_u,
        MAP    = map_u,
        credible_set = cred_set,
        HPD_set      = hpd_set,
        credible_interval = c(L = L, U = U),
        HPD_interval = range(hpd_set)
      ),
      counts_posterior = grid[, c("a","b","p","weight")],
      synthetic_effects = list(
        tau = list(mean = tau_sum$mean, CI = c(tau_sum$qL, tau_sum$qU)),
        RR  = list(mean = rr_sum$mean,  CI = c(rr_sum$qL,  rr_sum$qU)),
        OR  = list(mean = or_sum$mean,  CI = c(or_sum$qL,  or_sum$qU))
      )
    ))
    
  } else if (method == "sampling"){
    type <- prior$type
    
    build_1d_gamma <- function(N, tN, logprior_fun_for_axis){
      idx <- 0:N
      lp <- vapply(idx, function(z) logprior_fun_for_axis(z) + log_kappa(abs(tN - z)), 0.0)
      lse <- matrixStats::logSumExp(lp)
      prob <- exp(lp - lse)
      list(idx = idx, prob = prob)
    }
    
    sample_cat <- function(values, prob, R){
      prob <- pmax(prob, 0)
      prob <- prob / sum(prob)
      sample(values, size = R, replace = TRUE, prob = prob)
    }
    
    if (type %in% c("uniform", "beta_binom")){
      if (type == "uniform"){
        lp_a <- function(a) 0
        lp_b <- function(b) 0
      } else {
        par <- prior$params %||% list()
        alpha1 <- par$alpha1 %||% 1
        beta1  <- par$beta1  %||% 1
        alpha0 <- par$alpha0 %||% 1
        beta0  <- par$beta0  %||% 1
        lp_a <- function(a) lchoose(n1, a) + lbeta(a + alpha1, n1 - a + beta1) - lbeta(alpha1, beta1)
        lp_b <- function(b) lchoose(n0, b) + lbeta(b + alpha0, n0 - b + beta0) - lbeta(alpha0, beta0)
      }
      
      g11 <- build_1d_gamma(n1, t11, lp_a)
      g01 <- build_1d_gamma(n0, t01, lp_b)
      
      a_s <- sample_cat(g11$idx, g11$prob, R)
      b_s <- sample_cat(g01$idx, g01$prob, R)
      
      u_s <- mapply(p_frta, a_s, b_s)
      
    } else if (type == "common_rate"){
      geom_axis_prob <- function(N, tN){
        z  <- 0:N
        lp <- log_kappa(abs(tN - z))
        lse <- matrixStats::logSumExp(lp)
        exp(lp - lse)
      }
      qA_prob <- geom_axis_prob(n1, t11)
      qB_prob <- geom_axis_prob(n0, t01)
      
      a_cur <- sample(0:n1, 1, prob = qA_prob)
      b_cur <- sample(0:n0, 1, prob = qB_prob)
      lp_cur <- get_logprior(a_cur, b_cur)
      
      burn <- max(1000, ceiling(0.1 * R))
      thin <- 1
      out_a <- integer(R); out_b <- integer(R); out_u <- numeric(R)
      
      i <- 0; kept <- 0; iters <- 0
      while (kept < R){
        iters <- iters + 1
        a_new <- sample(0:n1, 1, prob = qA_prob)
        b_new <- sample(0:n0, 1, prob = qB_prob)
        lp_new <- get_logprior(a_new, b_new)
        
        if (log(runif(1)) < (lp_new - lp_cur)){
          a_cur <- a_new; b_cur <- b_new; lp_cur <- lp_new
        }
        if (iters > burn && ((iters - burn) %% thin == 0)){
          i <- i + 1
          out_a[i] <- a_cur
          out_b[i] <- b_cur
          out_u[i] <- p_frta(a_cur, b_cur)
          kept <- i
        }
        if (iters > 1e8) stop("MH did not finish; reduce R or check prior.")
      }
      a_s <- out_a; b_s <- out_b; u_s <- out_u
    } else {
      stop("Unknown prior$type. Use 'uniform', 'beta_binom', or 'common_rate'.")
    }

    u_key <- round(u_s, 12)
    agg <- aggregate(rep(1/length(u_key), length(u_key)),
                     by = list(u = u_key), FUN = sum)
    names(agg) <- c("u","mass")
    agg <- agg[order(agg$u), , drop = FALSE]
    agg$cdf <- cumsum(agg$mass)
    
    post_mean <- mean(u_s)
    med_u <- unname(stats::quantile(u_s, probs = 0.5, type = 1))
    idx_map <- which.max(agg$mass)
    map_u   <- agg$u[idx_map]
    
    L <- agg$u[which(agg$cdf >= alpha_CI/2)[1]]
    U <- agg$u[which(agg$cdf >= 1 - alpha_CI/2)[1]]
    cred_set <- agg$u[agg$u >= L & agg$u <= U]
    
    ord <- order(agg$mass, decreasing = TRUE)
    mass_cum <- cumsum(agg$mass[ord])
    k <- which(mass_cum >= (1 - alpha_CI))[1]
    hpd_set <- sort(agg$u[ord][seq_len(k)])
    
    tau_s <- a_s / n1 - b_s / n0
    
    need_cc_rr <- (a_s == 0L) | (b_s == 0L)
    risk_t <- ifelse(need_cc_rr, (a_s + 0.5) / (n1 + 1), a_s / n1)
    risk_c <- ifelse(need_cc_rr, (b_s + 0.5) / (n0 + 1), b_s / n0)
    rr_s <- risk_t / risk_c
    
    need_cc_or <- (a_s == 0L) | (b_s == 0L) | ((n1 - a_s) == 0L) | ((n0 - b_s) == 0L)
    a_cc <- ifelse(need_cc_or, a_s + 0.5, a_s)
    b_cc <- ifelse(need_cc_or, b_s + 0.5, b_s)
    c_cc <- ifelse(need_cc_or, (n1 - a_s) + 0.5, (n1 - a_s))
    d_cc <- ifelse(need_cc_or, (n0 - b_s) + 0.5, (n0 - b_s))
    or_s <- (a_cc * d_cc) / (c_cc * b_cc)
    
    qfun <- function(x) unname(stats::quantile(x, probs = c(alpha_CI/2, 1 - alpha_CI/2), type = 1))
    
    return(list(
      method   = "sampling",
      meta     = list(n1 = n1, n0 = n0, n11 = n11, n01 = n01,
                      clip = clip, epsilon = epsilon, rho = rho,
                      tilde = c(t11 = t11, t01 = t01),
                      prior  = prior, alpha_CI = alpha_CI, R = R, seed = seed),
      posterior_p = agg[, c("u","mass","cdf")],
      summaries = list(
        mean   = post_mean,
        median = med_u,
        MAP    = map_u,
        credible_set = cred_set,
        HPD_set      = hpd_set,
        credible_interval = c(L = L, U = U),
        HPD_interval = range(hpd_set)
      ),
      samples = data.frame(a = a_s, b = b_s, p = u_s,
                           tau = tau_s, RR = rr_s, OR = or_s),
      synthetic_effects = list(
        tau = list(mean = mean(tau_s),  CI = qfun(tau_s)),
        RR  = list(mean = mean(rr_s),   CI = qfun(rr_s)),
        OR  = list(mean = mean(or_s),   CI = qfun(or_s))
      )
    ))
  } else {
    stop("method must be 'enumerate' or 'sampling'.")
  }
}



# ==============================================================================
# DP-FRT-p: Direct Perturbation of the Exact FRT p-value (Laplace)
# ==============================================================================

dp_frt_p <- function(
    n1, n0, n11, n01, epsilon,
    seed  = NULL
){
  # ---------- Input checks ----------
  stopifnot(n11 <= n1, n01 <= n0, n1 > 0, n0 > 0)
  stopifnot(is.finite(epsilon), epsilon > 0)
  n <- n1 + n0
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---------- Non-private one-sided FRT p-value ----------
  # Under sharp null, treated successes ~ Hypergeometric(n, K, n1), K = n+1 = n11 + n01
  K <- n11 + n01
  p_frt <- phyper(q = n11 - 1, m = K, n = n - K, k = n1, lower.tail = FALSE)
  
  # ---------- Sensitivity (Lemma 3.1): Δp = max{n1/n, n0/n} ----------
  delta_p <- max(n1 / n, n0 / n)
  
  # ---------- Laplace mechanism ----------
  rlaplace <- function(m, scale){
    u <- runif(m, min = -0.5, max = 0.5)
    -scale * sign(u) * log(1 - 2 * abs(u))
  }
  eta <- rlaplace(1, scale = delta_p / epsilon)
  
  # ---------- Clip to feasible p-value range [|Z|^{-1}, 1] ----------
  L <- 1 / choose(n, n1)
  U <- 1
  p_tilde_raw <- p_frt + eta
  p_tilde <- min(U, max(L, p_tilde_raw))
  
  return(list(
    method = "DP-FRT-p",
    meta = list(
      n1 = n1, n0 = n0, n11 = n11, n01 = n01, n = n, K = K,
      epsilon = epsilon,
      sensitivity = delta_p,
      clip_range = c(L = L, U = U),
      noise = list(type = "Laplace", value = unname(eta), scale = delta_p / epsilon),
      seed = seed
    ),
    non_private = list(p_frt = unname(p_frt)),
    private = list(p_tilde = unname(p_tilde), p_tilde_raw = unname(p_tilde_raw))
  ))
}



# ==============================================================================
# DP-FRT-t: Perturb Test Statistic + Privatize Randomization Distribution
#          (Laplace on tau_hat, Geometric on n_{+1})
# ==============================================================================

dp_frt_t <- function(
    n1, n0, n11, n01, epsilon,
    epsilon_obs = NULL,   # default epsilon/2
    epsilon_ref = NULL,   # default epsilon/2
    seed = NULL
){
  # ---------- Input checks ----------
  stopifnot(n11 <= n1, n01 <= n0, n1 > 0, n0 > 0)
  stopifnot(is.finite(epsilon), epsilon > 0)
  n <- n1 + n0
  
  # split budget
  epsilon_obs <- epsilon_obs %||% (epsilon / 2)
  epsilon_ref <- epsilon_ref %||% (epsilon / 2)
  stopifnot(epsilon_obs > 0, epsilon_ref > 0)
  stopifnot(abs((epsilon_obs + epsilon_ref) - epsilon) < 1e-12)
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---------- Utilities ----------
  rlaplace <- function(m, scale){
    u <- runif(m, min = -0.5, max = 0.5)
    -scale * sign(u) * log(1 - 2 * abs(u))
  }
  rgeom_two_sided <- function(size, rho){
    stopifnot(rho > 0, rho < 1)
    rgeom(size, prob = 1 - rho) - rgeom(size, prob = 1 - rho)
  }
  
  # ---------- Non-private statistic ----------
  tau_hat <- n11 / n1 - n01 / n0
  
  # ---------- Sensitivity (Lemma 3.3): Δtau = max{1/n1, 1/n0} ----------
  delta_tau <- max(1 / n1, 1 / n0)
  
  # ---------- Perturb observed statistic (Laplace), clip to [-1,1] ----------
  eta_obs <- rlaplace(1, scale = delta_tau / epsilon_obs)
  T_tilde_raw <- tau_hat + eta_obs
  T_tilde <- min(1, max(-1, T_tilde_raw))
  
  # ---------- Privatize n_{+1} with two-sided geometric (Δ=1) ----------
  K <- n11 + n01
  rho_ref <- exp(-epsilon_ref)
  eta_ref <- rgeom_two_sided(1, rho_ref)
  K_tilde_raw <- K + eta_ref
  K_tilde <- max(0L, min(n, as.integer(round(K_tilde_raw))))
  
  # ---------- Compute privatized Fisher tail probability ----------
  # t ~ Hypergeometric(n, K_tilde, n1), reject region: t/n1 - (K_tilde - t)/n0 >= T_tilde
  a <- max(0L, K_tilde - n0)
  b <- min(n1, K_tilde)
  
  # monotone in t, solve threshold
  denom <- (1 / n1 + 1 / n0)
  tcrit <- ceiling((T_tilde + K_tilde / n0) / denom)
  
  if (tcrit <= a){
    p_tilde <- 1
  } else if (tcrit > b){
    p_tilde <- 0
  } else {
    p_tilde <- phyper(q = tcrit - 1, m = K_tilde, n = n - K_tilde, k = n1, lower.tail = FALSE)
  }
  
  return(list(
    method = "DP-FRT-t",
    meta = list(
      n1 = n1, n0 = n0, n11 = n11, n01 = n01, n = n, K = K,
      epsilon = epsilon, epsilon_obs = epsilon_obs, epsilon_ref = epsilon_ref,
      sensitivities = list(delta_tau = delta_tau, delta_K = 1),
      noisy = list(
        T_tilde = unname(T_tilde), T_tilde_raw = unname(T_tilde_raw),
        K_tilde = as.integer(K_tilde), K_tilde_raw = unname(K_tilde_raw)
      ),
      noise = list(
        obs = list(type = "Laplace", value = unname(eta_obs), scale = delta_tau / epsilon_obs),
        ref = list(type = "TwoSidedGeometric", value = unname(eta_ref), rho = rho_ref)
      ),
      seed = seed
    ),
    non_private = list(
      tau_hat = unname(tau_hat),
      K = as.integer(K)
    ),
    private = list(
      p_tilde = unname(p_tilde),
      threshold = list(a = a, b = b, tcrit = tcrit)
    )
  ))
}