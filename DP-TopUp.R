# ==============================================================================
# DP-TopUp: Sequential Decision under Additional Privacy Budget
# ------------------------------------------------------------------------------
# Separation of roles:
#   (A) Data-holder side: generate second DP release with additional epsilon_plus
#   (B) Analyst side: update posterior using (tilde, tilde_plus) and decide
#
# Implements:
#   (1) epsilon_plus lower bound via a data-adaptive Delta
#   (2) Posterior update with two independent releases
#   (3) Bayes-optimal decision with abstention using dp_frt_decision()
# ==============================================================================

if (!requireNamespace("this.path", quietly = TRUE)) install.packages("this.path")
current_dir <- dirname(this.path::this.path())

source(file.path(current_dir, "DP-FRT.R"))
source(file.path(current_dir, "DP-Decision.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Two-sided geometric noise: Geom(prob=1-rho) - Geom(prob=1-rho)
.dpfrt_rgeom_two_sided <- function(size, rho) {
  stopifnot(rho > 0, rho < 1)
  stats::rgeom(size, prob = 1 - rho) - stats::rgeom(size, prob = 1 - rho)
}

# log kappa for two-sided geometric kernel at |h| with epsilon
.dpfrt_log_kappa_abs <- function(h_abs, epsilon) {
  rho <- exp(-epsilon)
  log_rho   <- -epsilon
  log_const <- log1p(-rho) - log1p(rho)
  log_const + h_abs * log_rho
}

# log P(T = t | true = z) for two-sided geometric + optional clipping to [0, N]
.dpfrt_log_kappa_clipped <- function(t, z, N, epsilon) {
  rho <- exp(-epsilon)
  
  if (t <= 0L) {
    return(z * log(rho) - log1p(rho))
  }
  if (t >= N) {
    return((N - z) * log(rho) - log1p(rho))
  }
  
  .dpfrt_log_kappa_abs(abs(t - z), epsilon)
}

# One-sided Fisher randomization p-value under sharp null for candidate counts (a,b)
.dpfrt_p_frta <- function(a, b, n1, n0) {
  n <- n1 + n0
  K <- a + b
  stats::phyper(q = a - 1, m = K, n = n - K, k = n1, lower.tail = FALSE)
}

# Binary search: for fixed K, find smallest a in [a_min, a_max] s.t. p(a) <= alpha.
# If none, return Inf.
.dpfrt_find_a_crit <- function(K, n1, n0, alpha) {
  n <- n1 + n0
  a_min <- max(0L, K - n0)
  a_max <- min(n1, K)
  if (a_min > a_max) return(Inf)
  
  p_at_max <- stats::phyper(q = a_max - 1, m = K, n = n - K, k = n1, lower.tail = FALSE)
  if (p_at_max > alpha) return(Inf)
  
  p_at_min <- stats::phyper(q = a_min - 1, m = K, n = n - K, k = n1, lower.tail = FALSE)
  if (p_at_min <= alpha) return(a_min)
  
  lo <- a_min
  hi <- a_max
  while (lo < hi) {
    mid <- as.integer(floor((lo + hi) / 2))
    p_mid <- stats::phyper(q = mid - 1, m = K, n = n - K, k = n1, lower.tail = FALSE)
    if (p_mid <= alpha) {
      hi <- mid
    } else {
      lo <- mid + 1L
    }
  }
  lo
}

# s(x)=tanh(x/2)
.dpfrt_s_tanh <- function(x) tanh(x / 2)

# Data-adaptive Delta(tilde n; eps_plus):
#   Delta = E_{X~mu1, Y~mu0}[ s(eps_plus * d(X,Y)) ],
# where mu1 is posterior over (a,b) restricted to S1 (p<=alpha),
# and mu0 is posterior restricted to S0 (p>alpha).
#
# Works with:
#   - enumerate output: frt_fit$counts_posterior has a,b,p,weight
#   - sampling output:  frt_fit$samples has a,b,p (equal weights)
.dpfrt_compute_Delta <- function(
    frt_fit,
    alpha,
    epsilon_plus,
    M = 20000,
    seed = NULL
){
  if (!is.null(seed)) set.seed(seed)
  
  # Enumerate case: weighted grid
  if (!is.null(frt_fit$counts_posterior)) {
    df <- frt_fit$counts_posterior
    a <- df$a; b <- df$b; p <- df$p; w <- df$weight
    w <- pmax(w, 0); w <- w / sum(w)
    
    idx1 <- which(p <= alpha)
    idx0 <- which(p >  alpha)
    if (length(idx1) == 0 || length(idx0) == 0) {
      return(list(Delta = 0, details = list(mode="weighted_MC", M=M, n_S1=length(idx1), n_S0=length(idx0))))
    }
    
    mu1 <- w[idx1]; mu1 <- mu1 / sum(mu1)
    mu0 <- w[idx0]; mu0 <- mu0 / sum(mu0)
    
    i <- sample(idx1, size = M, replace = TRUE, prob = mu1)
    j <- sample(idx0, size = M, replace = TRUE, prob = mu0)
    
    d <- abs(a[i] - a[j]) + abs(b[i] - b[j])
    Delta <- mean(.dpfrt_s_tanh(epsilon_plus * d))
    
    return(list(
      Delta = Delta,
      details = list(mode = "weighted_MC", M = M,
                     mass_S1 = sum(w[idx1]), mass_S0 = sum(w[idx0]),
                     n_S1 = length(idx1), n_S0 = length(idx0))
    ))
  }
  
  # Sampling case: equal-weight draws
  if (!is.null(frt_fit$samples)) {
    s <- frt_fit$samples
    idx1 <- which(s$p <= alpha)
    idx0 <- which(s$p >  alpha)
    if (length(idx1) == 0 || length(idx0) == 0) {
      return(list(Delta = 0, details = list(mode="samples_MC", M=M, n_S1=length(idx1), n_S0=length(idx0))))
    }
    
    i <- sample(idx1, size = M, replace = TRUE)
    j <- sample(idx0, size = M, replace = TRUE)
    
    d <- abs(s$a[i] - s$a[j]) + abs(s$b[i] - s$b[j])
    Delta <- mean(.dpfrt_s_tanh(epsilon_plus * d))
    
    return(list(
      Delta = Delta,
      details = list(mode="samples_MC", M=M, n_S1=length(idx1), n_S0=length(idx0))
    ))
  }
  
  stop("frt_fit must contain either $counts_posterior (enumerate) or $samples (sampling).")
}

# ------------------------------------------------------------------------------
# Analyst-side posterior: two independent DP releases (NO true counts needed)
# ------------------------------------------------------------------------------

.dpfrt_posterior_two_release_public <- function(
    n1, n0,
    t11, t01, epsilon,
    t11_plus, t01_plus, epsilon_plus,
    prior,
    clip = FALSE,
    clip_plus = FALSE,
    method = c("sampling","enumerate"),
    R = 10000,
    alpha_CI = 0.05,
    seed = NULL
) {
  method <- match.arg(method, choices = c("sampling","enumerate"))
  stopifnot(is.finite(epsilon) && epsilon > 0,
            is.finite(epsilon_plus) && epsilon_plus > 0,
            n1 > 0, n0 > 0)
  
  if (!is.null(seed)) set.seed(seed)
  
  # prior helper
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
      n  <- n1 + n0
      K  <- a + b
      return(lchoose(n1, a) + lchoose(n0, b) +
               lbeta(K + alpha, n - K + beta) - lbeta(alpha, beta))
      
    } else {
      stop("Unknown prior$type. Use 'uniform', 'beta_binom', or 'common_rate'.")
    }
  }
  
  p_frta <- function(a, b) .dpfrt_p_frta(a, b, n1, n0)
  
  loglik_ab <- function(a, b){
    ll11 <- if (clip)      .dpfrt_log_kappa_clipped(t11, a, n1, epsilon)
            else           .dpfrt_log_kappa_abs(abs(t11 - a), epsilon)
    
    ll01 <- if (clip)      .dpfrt_log_kappa_clipped(t01, b, n0, epsilon)
            else           .dpfrt_log_kappa_abs(abs(t01 - b), epsilon)
    
    ll11p <- if (clip_plus) .dpfrt_log_kappa_clipped(t11_plus, a, n1, epsilon_plus)
             else           .dpfrt_log_kappa_abs(abs(t11_plus - a), epsilon_plus)
    
    ll01p <- if (clip_plus) .dpfrt_log_kappa_clipped(t01_plus, b, n0, epsilon_plus)
             else           .dpfrt_log_kappa_abs(abs(t01_plus - b), epsilon_plus)
    
    ll11 + ll01 + ll11p + ll01p
  }
  
  
  if (method == "enumerate"){
    A <- 0:n1
    B <- 0:n0
    
    logw_mat <- matrix(-Inf, nrow = n1 + 1, ncol = n0 + 1)
    p_mat    <- matrix(NA_real_, nrow = n1 + 1, ncol = n0 + 1)
    
    for (a in A){
      for (b in B){
        logw_mat[a + 1, b + 1] <- get_logprior(a,b) + loglik_ab(a,b)
        p_mat[a + 1, b + 1]    <- p_frta(a,b)
      }
    }
    
    if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
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
    
    post_mean <- sum(agg$u * agg$mass)
    med_u <- agg$u[which(agg$cdf >= 0.5)[1]]
    map_u <- agg$u[which.max(agg$mass)]
    
    L <- agg$u[which(agg$cdf >= alpha_CI/2)[1]]
    U <- agg$u[which(agg$cdf >= 1 - alpha_CI/2)[1]]
    
    ord <- order(agg$mass, decreasing = TRUE)
    mass_cum <- cumsum(agg$mass[ord])
    k <- which(mass_cum >= (1 - alpha_CI))[1]
    hpd_set <- sort(agg$u[ord][seq_len(k)])
    
    return(list(
      method = "enumerate",
      meta = list(
        n1 = n1, n0 = n0,
        epsilon_total   = epsilon + epsilon_plus,
        epsilon_initial = epsilon,
        epsilon_plus    = epsilon_plus,
        tilde      = c(t11 = t11, t01 = t01),
        tilde_plus = c(t11_plus = t11_plus, t01_plus = t01_plus),
        prior = prior,
        alpha_CI = alpha_CI
      ),
      posterior_p = agg[, c("u","mass","cdf")],
      summaries = list(
        mean   = post_mean,
        median = med_u,
        MAP    = map_u,
        credible_interval = c(L = L, U = U),
        HPD_interval = range(hpd_set),
        HPD_set = hpd_set
      ),
      counts_posterior = grid[, c("a","b","p","weight")]
    ))
  }
  
  # ---------------- sampling ----------------
  type <- prior$type
  
  sample_cat <- function(values, prob, R){
    prob <- pmax(prob, 0)
    prob <- prob / sum(prob)
    sample(values, size = R, replace = TRUE, prob = prob)
  }
  
  build_1d_gamma_two <- function(N, tN, eps1, tN2, eps2, logprior_axis){
    idx <- 0:N
    lp <- vapply(
      idx,
      function(z){
        logprior_axis(z) +
          .dpfrt_log_kappa_abs(abs(tN  - z), eps1) +
          .dpfrt_log_kappa_abs(abs(tN2 - z), eps2)
      },
      0.0
    )
    if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
    lse <- matrixStats::logSumExp(lp)
    prob <- exp(lp - lse)
    list(idx = idx, prob = prob)
  }
  
  if (type %in% c("uniform","beta_binom")){
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
    
    g11 <- build_1d_gamma_two(n1, t11, epsilon, t11_plus, epsilon_plus, lp_a)
    g01 <- build_1d_gamma_two(n0, t01, epsilon, t01_plus, epsilon_plus, lp_b)
    
    a_s <- sample_cat(g11$idx, g11$prob, R)
    b_s <- sample_cat(g01$idx, g01$prob, R)
    u_s <- mapply(p_frta, a_s, b_s)
    
  } else if (type == "common_rate") {
    
    geom_axis_prob_two <- function(N, tN, eps1, tN2, eps2){
      z <- 0:N
      lp <- .dpfrt_log_kappa_abs(abs(tN  - z), eps1) +
        .dpfrt_log_kappa_abs(abs(tN2 - z), eps2)
      if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
      lse <- matrixStats::logSumExp(lp)
      exp(lp - lse)
    }
    
    qA <- geom_axis_prob_two(n1, t11, epsilon, t11_plus, epsilon_plus)
    qB <- geom_axis_prob_two(n0, t01, epsilon, t01_plus, epsilon_plus)
    
    a_cur <- sample(0:n1, 1, prob = qA)
    b_cur <- sample(0:n0, 1, prob = qB)
    lp_cur <- get_logprior(a_cur, b_cur)
    
    burn <- max(1000, ceiling(0.1 * R))
    out_a <- integer(R); out_b <- integer(R); out_u <- numeric(R)
    
    kept <- 0; iters <- 0
    while (kept < R){
      iters <- iters + 1
      a_new <- sample(0:n1, 1, prob = qA)
      b_new <- sample(0:n0, 1, prob = qB)
      lp_new <- get_logprior(a_new, b_new)
      
      if (log(runif(1)) < (lp_new - lp_cur)){
        a_cur <- a_new; b_cur <- b_new; lp_cur <- lp_new
      }
      if (iters > burn){
        kept <- kept + 1
        out_a[kept] <- a_cur
        out_b[kept] <- b_cur
        out_u[kept] <- p_frta(a_cur, b_cur)
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
  map_u <- agg$u[which.max(agg$mass)]
  
  L <- agg$u[which(agg$cdf >= alpha_CI/2)[1]]
  U <- agg$u[which(agg$cdf >= 1 - alpha_CI/2)[1]]
  
  return(list(
    method = "sampling",
    meta = list(
      n1 = n1, n0 = n0,
      epsilon_total   = epsilon + epsilon_plus,
      epsilon_initial = epsilon,
      epsilon_plus    = epsilon_plus,
      tilde      = c(t11 = t11, t01 = t01),
      tilde_plus = c(t11_plus = t11_plus, t01_plus = t01_plus),
      prior = prior,
      alpha_CI = alpha_CI,
      R = R,
      seed = seed
    ),
    posterior_p = agg[, c("u","mass","cdf")],
    summaries = list(
      mean   = post_mean,
      median = med_u,
      MAP    = map_u,
      credible_interval = c(L = L, U = U)
    ),
    samples = data.frame(a = a_s, b = b_s, p = u_s)
  ))
}

# ------------------------------------------------------------------------------
# Main 1 (Analyst): epsilon_plus lower bound using data-adaptive Delta
# ------------------------------------------------------------------------------
dp_frt_epsilon_plus_lower_bound <- function(
    decision1,
    frt_fit = NULL,
    xi = 0.05,
    alpha = NULL,
    M = 20000,
    eps_upper_init = 1,
    seed = NULL
){
  frt_fit <- frt_fit %||% decision1$frt_fit %||% decision1$fit %||% NULL
  if (is.null(frt_fit)) stop("Need frt_fit (or store it in decision1$frt_fit).")
  
  # psi + abstention interval from decision1
  Psi <- decision1$psi_obs %||% stop("decision1$psi_obs not found.")
  tlow  <- decision1$threshold_low  %||% NA_real_
  thigh <- decision1$threshold_high %||% NA_real_
  if (!is.finite(tlow) || !is.finite(thigh)) {
    stop("Need decision1$threshold_low and decision1$threshold_high (Bayes abstain output).")
  }
  
  alpha <- alpha %||% decision1$alpha %||% 0.05
  
  # If already outside abstention region, no top-up needed
  if (Psi <= tlow || Psi >= thigh) {
    return(list(
      epsilon_plus_lower_bound = 0,
      inputs = list(xi = xi, alpha = alpha, Psi = Psi, tlow = tlow, thigh = thigh, rPsi = 0),
      Delta_details = NULL
    ))
  }
  
  rPsi <- min(Psi - tlow, thigh - Psi)
  if (!is.finite(rPsi) || rPsi <= .Machine$double.eps) {
    return(list(
      epsilon_plus_lower_bound = Inf,
      inputs = list(xi = xi, alpha = alpha, Psi = Psi, tlow = tlow, thigh = thigh, rPsi = rPsi),
      Delta_details = NULL
    ))
  }
  
  target <- 1 - xi
  coef <- 2 * Psi * (1 - Psi) / rPsi
  
  if (coef < target) {
    return(list(
      epsilon_plus_lower_bound = Inf,
      inputs = list(xi = xi, alpha = alpha, Psi = Psi, tlow = tlow, thigh = thigh, rPsi = rPsi,
                    coef = coef, target = target),
      Delta_details = list(note = "coef < 1-xi so even Delta<=1 cannot reach target; no finite lower bound.")
    ))
  }
  
  # Make g(eps) deterministic & monotone by using a fixed seed for Delta sampling
  seed_eff <- seed %||% 2025L
  
  g_with_details <- function(eps){
    res <- .dpfrt_compute_Delta(
      frt_fit = frt_fit,
      alpha = alpha,
      epsilon_plus = eps,
      M = M,
      seed = seed_eff
    )
    list(val = coef * res$Delta - target, details = res$details, Delta = res$Delta)
  }
  
  # Bracket smallest eps with g(eps) >= 0
  lo <- 0
  hi <- eps_upper_init
  out_hi <- g_with_details(hi)
  
  iter <- 0
  while (out_hi$val < 0 && iter < 60) {
    hi <- hi * 2
    out_hi <- g_with_details(hi)
    iter <- iter + 1
  }
  
  if (out_hi$val < 0) {
    return(list(
      epsilon_plus_lower_bound = Inf,
      inputs = list(xi = xi, alpha = alpha, Psi = Psi, tlow = tlow, thigh = thigh, rPsi = rPsi,
                    coef = coef, target = target, M = M, eps_upper_reached = hi),
      Delta_details = list(note = "Failed to bracket; increase eps_upper_init or M (Monte Carlo noise).",
                           last = out_hi)
    ))
  }
  
  # uniroot on deterministic, monotone g(eps)
  root <- stats::uniroot(
    f = function(eps) g_with_details(eps)$val,
    lower = lo,
    upper = hi
  )$root
  
  out_root <- g_with_details(root)
  
  list(
    epsilon_plus_lower_bound = root,
    inputs = list(
      xi = xi,
      alpha = alpha,
      Psi = Psi,
      tlow = tlow,
      thigh = thigh,
      rPsi = rPsi,
      coef = coef,
      target = target,
      M = M,
      seed = seed_eff
    ),
    Delta_at_root = out_root$Delta,
    Delta_details = out_root$details
  )
}

# ------------------------------------------------------------------------------
# Main 2 (Data-holder): Generate top-up DP release with epsilon_plus
# ------------------------------------------------------------------------------

dp_frt_topup_release <- function(
    n11, n01, n1, n0,
    epsilon_plus,
    clip = FALSE,
    seed = NULL
){
  stopifnot(is.finite(epsilon_plus) && epsilon_plus > 0)
  stopifnot(n11 >= 0, n11 <= n1, n01 >= 0, n01 <= n0)
  
  # NOTE: for real deployments, do not let analyst control RNG seed.
  if (!is.null(seed)) set.seed(seed)
  
  rho_plus <- exp(-epsilon_plus)
  eta11_plus <- .dpfrt_rgeom_two_sided(1, rho_plus)
  eta01_plus <- .dpfrt_rgeom_two_sided(1, rho_plus)
  
  t11_plus <- as.integer(n11 + eta11_plus)
  t01_plus <- as.integer(n01 + eta01_plus)
  
  if (clip){
    t11_plus <- max(0L, min(n1, t11_plus))
    t01_plus <- max(0L, min(n0, t01_plus))
  }
  
  list(
    tilde_plus = c(t11_plus = t11_plus, t01_plus = t01_plus),
    epsilon_plus = epsilon_plus,
    clip = clip
  )
}

# ------------------------------------------------------------------------------
# Main 3 (Analyst): Update posterior with top-up release and decide
# ------------------------------------------------------------------------------

dp_frt_update_with_topup_and_decide <- function(
    decision1,
    tilde_plus,
    epsilon_plus,
    frt_fit = NULL,
    alpha = NULL,
    bayes_abstain = TRUE,
    clip_plus = FALSE,
    lambda_0 = NULL,
    lambda_1 = NULL,
    lambda_u = NULL,
    method = NULL,
    R = 10000,
    alpha_CI = 0.05,
    seed = NULL
){
  frt_fit <- frt_fit %||% decision1$frt_fit %||% decision1$fit %||% NULL
  if (is.null(frt_fit)) stop("Need frt_fit (or store it in decision1$frt_fit).")
  
  meta <- frt_fit$meta
  n1 <- meta$n1; n0 <- meta$n0
  eps0 <- meta$epsilon
  prior <- meta$prior
  
  # first release (public)
  t11 <- as.integer(unname(meta$tilde["t11"]))
  t01 <- as.integer(unname(meta$tilde["t01"]))
  
  # second release (public, provided)
  t11_plus <- as.integer(unname(tilde_plus["t11_plus"]))
  t01_plus <- as.integer(unname(tilde_plus["t01_plus"]))
  if (anyNA(c(t11_plus, t01_plus))) {
    stop("tilde_plus must contain names 't11_plus' and 't01_plus'.")
  }
  
  if (is.null(method) || length(method) == 0L) {
    method <- frt_fit$method %||% "sampling"
  }
  method <- match.arg(method, choices = c("sampling","enumerate"))
  
  clip1 <- frt_fit$meta$clip %||% FALSE
  
  frt_fit_plus <- .dpfrt_posterior_two_release_public(
    n1 = n1, n0 = n0,
    t11 = t11, t01 = t01, epsilon = eps0,
    t11_plus = t11_plus, t01_plus = t01_plus, epsilon_plus = epsilon_plus,
    prior = prior,
    clip = clip1,
    clip_plus = clip_plus,
    method = method,
    R = R,
    alpha_CI = alpha_CI,
    seed = seed
  )
  
  alpha <- alpha %||% decision1$alpha %||% 0.05
  lambda_0 <- lambda_0 %||% decision1$lambda_0 %||% 1
  lambda_1 <- lambda_1 %||% decision1$lambda_1 %||% 1
  lambda_u <- lambda_u %||% decision1$lambda_u %||% NULL
  
  decision2 <- dp_frt_decision(
    frt_fit = frt_fit_plus,
    alpha = alpha,
    framework = "bayes",
    bayes_abstain = bayes_abstain,
    lambda_0 = lambda_0,
    lambda_1 = lambda_1,
    lambda_u = lambda_u
  )
  
  list(
    epsilon_plus = epsilon_plus,
    epsilon_total = eps0 + epsilon_plus,
    tilde_plus = c(t11_plus = t11_plus, t01_plus = t01_plus),
    frt_fit_plus = frt_fit_plus,
    decision2 = decision2
  )
}