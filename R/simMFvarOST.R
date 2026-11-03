#' Simulate a mixed-frequency VAR with centered orthogonal skew-t innovations
#'
#' Structural innovations are A^{-1}((w-mu)*gamma + sqrt(w)*epsilon),
#' where epsilon ~ N(0, diag(exp(h))), w_i ~ IG(nu_i/2, nu_i/2),
#' and mu_i = nu_i/(nu_i-2). This preserves the original OT recursion.
#'
#' @param K Number of variables.
#' @param p Number of lags (currently >= 1).
#' @param t_max Number of retained monthly periods.
#' @param b0 VAR coefficients, as in sim.MFVAR.OT.novol().
#' @param a0 Coefficients of the unit-lower-triangular A matrix, as in OT.
#' @param h Log baseline variances; sigma_i = exp(h_i/2).
#' @param y0 Initial p by K observations.
#' @param nu Equation-specific degrees of freedom, each > 4.
#' @param gamma Equation-specific skewness coefficients.
#' @param idq Indices of variables observed at quarterly frequency.
#' @param freq Spacing of observed quarterly data.
#' @param aggregation One of identity, average, or triangular.
#' @param seednum Random seed; NULL leaves the global generator unchanged.
#' @param burn_in Number of simulated periods discarded before the sample.
#' @param measurement_variance Variance of added noise in observed quarterly
#'   aggregates, for average/triangular only. Set equal to
#'   prior$sdeviation to match the estimator; default 0.001.
#' @return A list with observed y, latent y_true, aggregated y_atrue,
#'   true B0, A0, sigma, sigma2, gamma, nu, w, and sample initial y0.
#' @export
sim.MFVAR.OST.novol <- function(K = 5, p = 2, t_max = 1000,
                               b0 = 0.5, a0 = 0.5, h = 0,
                               y0 = matrix(0, nrow = p, ncol = K),
                               nu = 8, gamma = 0,
                               idq = integer(0), freq = 3,
                               aggregation = "identity",
                               seednum = 0, burn_in = 0,
                               measurement_variance = 0.001) {
  if (!is.null(seednum)) set.seed(seednum)
  if (length(K) != 1 || !is.finite(K) || K < 1 || K != as.integer(K) ||
      length(p) != 1 || !is.finite(p) || p < 1 || p != as.integer(p) ||
      length(t_max) != 1 || !is.finite(t_max) || t_max < 1 ||
      t_max != as.integer(t_max) || length(burn_in) != 1 ||
      !is.finite(burn_in) || burn_in < 0 || burn_in != as.integer(burn_in)) {
    stop("K and p must be positive integers; t_max > 0 and burn_in >= 0.")
  }
  if (!aggregation %in% c("identity", "average", "triangular"))
    stop("Invalid aggregation.")
  if (length(freq) != 1 || !is.finite(freq) || freq < 1 ||
      freq != as.integer(freq)) stop("freq must be a positive integer.")
  # Legacy simulation wrapper uses rep(0, K) to mean no quarterly variables.
  if (length(idq) && all(idq == 0)) idq <- integer(0)
  if (anyNA(idq) || any(!idq %in% seq_len(K)) || anyDuplicated(idq))
    stop("idq must contain unique indices between 1 and K.")
  if (length(h) == 1) h <- rep(h, K)
  if (length(nu) == 1) nu <- rep(nu, K)
  if (length(gamma) == 1) gamma <- rep(gamma, K)
  if (length(h) != K || length(nu) != K || length(gamma) != K ||
      any(!is.finite(c(h, nu, gamma))) || any(nu <= 4))
    stop("h, nu, gamma must be finite length-K vectors and nu > 4.")
  if (!is.matrix(y0) || !identical(dim(y0), c(as.integer(p), as.integer(K))) ||
      any(!is.finite(y0))) stop("y0 must be a finite p by K matrix.")
  if (length(measurement_variance) != 1 ||
      !is.finite(measurement_variance) || measurement_variance < 0)
    stop("measurement_variance must be a nonnegative finite scalar.")

  if (length(b0) == 1) {
    B0 <- cbind(rep(0, K), do.call(cbind, lapply(seq_len(p), function(j)
      b0^j * diag(K))))
  } else {
    B0 <- matrix(b0, nrow = K)
  }
  if (!identical(dim(B0), c(as.integer(K), as.integer(1 + K*p))) ||
      any(!is.finite(B0)))
    stop("b0 must define a K by (1+K*p) coefficient matrix.")

  A0 <- diag(K)
  if (length(a0) == 1) {
    A0[lower.tri(A0)] <- a0
  } else if (length(a0) == K*(K-1)/2) {
    A0[lower.tri(A0)] <- a0
  } else {
    stop("a0 must be scalar or a vector of K*(K-1)/2 coefficients.")
  }
  if (any(!is.finite(A0))) stop("a0 must be finite.")

  total <- t_max + burn_in
  w_t <- matrix(NA_real_, nrow = total, ncol = K)
  for (k in seq_len(K)) {
    w_t[, k] <- invgamma::rinvgamma(total, shape = nu[k]/2,
                                    rate = nu[k]/2)
  }
  eps <- matrix(rnorm(total*K), nrow = total, ncol = K)
  mu <- nu/(nu-2)
  sigma <- exp(h/2)
  inv_A0 <- solve(A0)
  ystar <- y0
  y_mean <- matrix(NA_real_, nrow = total, ncol = K)
  volatility <- matrix(NA_real_, nrow = total, ncol = K)
  q_t <- sweep(w_t, 2, mu, "-")
  q_t <- sweep(q_t, 2, gamma, "*")

  for (tt in seq_len(total)) {
    lag_rows <- (p + tt - 1):tt
    xt <- c(1, as.vector(t(ystar[lag_rows, , drop = FALSE])))
    y_mean[tt, ] <- as.vector(B0 %*% xt)
    S_t <- inv_A0 %*% diag(sqrt(w_t[tt, ]) * sigma, K)
    volatility[tt, ] <- diag(S_t %*% t(S_t))
    ysim <- as.vector(B0 %*% xt + inv_A0 %*% q_t[tt, ] +
                        S_t %*% eps[tt, ])
    ystar <- rbind(ystar, ysim)
  }

  retained <- p + burn_in + seq_len(t_max)
  kept_draws <- burn_in + seq_len(t_max)
  y_true <- ystar[retained, , drop = FALSE]
  y0_sample <- ystar[burn_in + seq_len(p), , drop = FALSE]
  y_atrue <- y_true
  if (aggregation %in% c("average", "triangular") && length(idq)) {
    aggregate_all <- as.matrix(Make_bandSparse(nrow(ystar), aggregation) %*%
                                 ystar[, idq, drop = FALSE])
    y_atrue[, idq] <- aggregate_all[retained, , drop = FALSE]
  }
  y <- y_atrue
  observed <- seq(1L, t_max, by = freq)
  if (length(idq)) {
    if (aggregation %in% c("average", "triangular") &&
        measurement_variance > 0) {
      y[observed, idq] <- y[observed, idq, drop = FALSE] +
        matrix(rnorm(length(observed)*length(idq),
                     sd = sqrt(measurement_variance)),
               nrow = length(observed), ncol = length(idq))
    }
    not_observed <- setdiff(seq_len(t_max), observed)
    if (length(not_observed)) y[not_observed, idq] <- NA_real_
  }

  list(y = y, y_true = y_true, y_atrue = y_atrue,
       y0 = y0_sample, B0 = B0, A0 = A0, h = h,
       sigma = sigma, sigma2 = sigma^2, gamma = gamma, nu = nu,
       w = w_t[kept_draws, , drop = FALSE],
       q = q_t[kept_draws, , drop = FALSE],
       y_mean = y_mean[kept_draws, , drop = FALSE],
       volatility = volatility[kept_draws, , drop = FALSE],
       K = K, p = p, t_max = t_max, idq = idq, freq = freq,
       aggregation = aggregation, measurement_variance = measurement_variance,
       dist = "OST", SV = FALSE)
}
