#' @export
BMFVAR.OST.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OST", 
                       SV = FALSE, aggregation = "identity", idq = K)
  }
  
  if (!identical(prior$dist, "OST")) {
    stop("This sampler requires prior$dist == 'OST'.")
  }
  if (isTRUE(prior$SV)) {
    stop("This version implements OST without stochastic volatility.")
  }
  if (!all(dim(y) == c(nrow(y), K))) stop("K must equal ncol(y).")
  
  aggregation <- prior$aggregation
  idq <- prior$idq
  sdeviation <- prior$sdeviation
  ymiss_invprior <- 1/prior$ymiss_prior
  
  if (aggregation == "identity"){
    yt_vec <- vec(t(y))
    yobs_vec <- na.exclude(yt_vec)  
  }
  
  if (aggregation == "average" || aggregation == "triangular"){
    y[,idq] <- NA
    
    yt_vec <- vec(t(y))
    yobs_vec <- na.exclude(yt_vec)  
    z_obs_vec <- na.exclude(vec(t(y_raw[,idq])))  
    
    y_temp <- y
    y_temp[,idq] <- 0
    y_temp <- matrix(t(y_temp)[is.na(t(y))], ncol = 1)
    S_select_temp <- Make_SectMat(y_temp)
    # Account for the missing in other columns rather than in low freq columns
    M_aggregate <- kronecker(Make_bandSparse(t_max, aggregation)[!is.na(y_raw[,idq[1] ]),],
                             Diagonal(length(idq))) %*% Matrix::t(S_select_temp$S_obs)
  }
  
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  # prior for skewness parameters gamma
  gamma_prior <- prior$gamma_prior
  V_gamma_prior <- prior$V_gamma_prior
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  sigma2 <- sigma^2
  
  V_b_prior_inv <- solve(V_b_prior)
  
  # Multi degrees of freedom
  nu <- inits$nu
  gamma <- inits$gamma
  mu <- nu / (nu - 2)
  if (length(nu) != K || any(!is.finite(nu)) || any(nu <= 4)) {
    stop("inits$nu must contain K values strictly greater than 4.")
  }
  if (length(gamma) != K || any(!is.finite(gamma))) {
    stop("inits$gamma must contain K finite values.")
  }
  if (!isTRUE(all.equal(V_gamma_prior,
                        diag(diag(V_gamma_prior), nrow = K)))) {
    stop("Separate gamma updates require a diagonal prior covariance matrix.")
  }
  # Fixed random-walk proposal SD for nu (can be tuned via inits$nu_rw_sd).
  nu_rw_sd <- inits$nu_rw_sd
  if (is.null(nu_rw_sd)) nu_rw_sd <- rep(1, K)
  if (length(nu_rw_sd) == 1L) nu_rw_sd <- rep(nu_rw_sd, K)
  if (length(nu_rw_sd) != K || any(!is.finite(nu_rw_sd)) ||
      any(nu_rw_sd <= 0)) stop("nu_rw_sd must be positive for each equation.")
  
  # Initialize the scale mixture variables.
  w <- matrix(1, nrow = K, ncol = t_max)
  w_sqrt <- sqrt(w)
  
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior
  S_select <- Make_SectMat(y)
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  n_cond <- length(y_raw[,idq]) - sum(is.na(y_raw[,idq]))
  S_obs <- S_select$S_obs
  S_miss <- S_select$S_miss
  if (aggregation %in% c("average", "triangular")) {
    K_prior_y_miss <- Matrix::t(M_aggregate) %*%
      Diagonal(n_cond, 1/sdeviation) %*%
      M_aggregate + Diagonal(n_miss, ymiss_invprior)
    
    b_prior_y_miss <- Matrix::t(M_aggregate) %*%
      Diagonal(n_cond, 1/sdeviation) %*%
      z_obs_vec
  }
  
  # Output: B, A, sigma, gamma, nu, w (and missing values separately).
  if (samples <= inits$burnin || inits$thin < 1 ||
      inits$burnin < 0 || (samples - inits$burnin) %/% inits$thin < 1) {
    stop("Specify samples > burnin >= 0 and 1 <= thin <= samples - burnin.")
  }
  n_save <- (samples - inits$burnin) %/% inits$thin
  mcmc <- matrix(NA_real_, nrow = m*K + K*(K-1)/2 + 3*K + K*t_max,
                 ncol = n_save)
  My_miss <- y_miss_mu <- y_miss_sd <- matrix(NA, nrow = n_miss,
                                              ncol = n_save)
  id_miss <- seq_len(t_max*K)[is.na(t(y))] # Location of Y_miss
  id_miss_t <- rep(1:t_max, each = K)[id_miss] # Time of Y_miss
  id_miss_k <- rep(1:K, t_max)[id_miss]       # Variable k of Y_miss
  id_m_m <- which(!id_miss_k %in% idq)
  id_m_q <- which(id_miss_k %in% idq)
  if (aggregation %in% c("average", "triangular")) {
    M_forecast <- kronecker(
      Make_bandSparse(t_max, aggregation)[is.na(y[, idq[1]]), ],
      Diagonal(length(idq))
    ) %*% Matrix::t(S_select_temp$S_obs)
  }
  
  for (j in c(1:samples)){
    
    # Sample y missing
    H_B <- Make_HB(B, K, p ,t_max)
    c_B <- Make_cB(y0, B, K, p ,t_max)
    IA <- kronecker(Diagonal(t_max), as(A,Class = "TsparseMatrix") )
    mu <- nu / (nu - 2)
    q <- sweep(w, 1, mu, "-") * matrix(gamma, nrow = K, ncol = t_max)
    c_const <- IA %*% c_B + as.vector(q)
    IAH <- IA %*% H_B
    G_obs <-  IAH %*% S_obs
    G_miss <-  IAH %*% S_miss
    Wh <- 1/vec(w) * rep(1/sigma^2, t_max)
    if (n_miss == 0L) {
      ym_mu_post <- ym_sd_post <- ymiss_vec <- numeric(0)
    } else {
      K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% G_miss
      b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*%
        (c_const - G_obs %*% yobs_vec)
    }
    
    if (n_miss > 0L && aggregation == "identity"){
      # sample from N( K^-1 b, K^-1)
      U.chol <- Matrix::chol( K_y_miss )
      # ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_y_miss,
      #                                          upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      ym_mu_post <- backsolve(U.chol, backsolve(U.chol, b_y_miss, upper.tri = T, transpose = T )) 
      ym_Sigma_post <- backsolve(U.chol, backsolve(U.chol, diag(n_miss), upper.tri = T, transpose = T)) 
      ymiss_vec <- ym_mu_post + backsolve(U.chol, rnorm(n_miss) )
      ym_sd_post <- sqrt(diag(ym_Sigma_post))
    }
    if (n_miss > 0L && aggregation %in% c("average", "triangular")){
      K_bar_y_miss <- K_y_miss + K_prior_y_miss
      b_bar_y_miss <- b_y_miss + b_prior_y_miss
      U.chol <- Matrix::chol( K_bar_y_miss )
      ym_mu_post <- backsolve(U.chol, backsolve(U.chol, b_bar_y_miss, upper.tri = T, transpose = T )) 
      ym_Sigma_post <- backsolve(U.chol, backsolve(U.chol, diag(n_miss), upper.tri = T, transpose = T)) 
      ymiss_vec <- ym_mu_post + backsolve(U.chol, rnorm(n_miss) )
      #Forecast mean and sd
      ym_sd_post <- sqrt(diag(ym_Sigma_post))
      ym_sd_post[id_m_q] <- sqrt(Matrix::diag(M_forecast %*% ym_Sigma_post %*% Matrix::t(M_forecast)))
      ym_mu_post[id_m_q] <- M_forecast %*% ym_mu_post
    }
    
    yt_vec <- S_obs%*% yobs_vec + S_miss%*% ymiss_vec
    
    
    # Sample B
    yt <- matrix(yt_vec, nrow = K)
    y <- t(yt)
    xt <- makeRegressor(y, y0, t_max, K, p)
    
    # The structural response is A*y_t - q_t, q_it=(w_it-mu_i)*gamma_i.
    q <- sweep(w, 1, mu, "-") * matrix(gamma, nrow = K, ncol = t_max)
    wt <- as.vector(1 / (w_sqrt * matrix(sigma, K, t_max)))
    y.tilde <- as.vector(A %*% yt - q) * wt
    x.tilde <- kronecker(t(xt), A) * wt
    theta.prec.chol <- chol(V_b_prior_inv + crossprod(x.tilde))
    b_sample <- backsolve(
      theta.prec.chol,
      backsolve(theta.prec.chol,
                theta.prior.precmean + crossprod(x.tilde, y.tilde),
                upper.tri = TRUE, transpose = TRUE) + rnorm(K*m)
    )
    B <- matrix(b_sample, K, m)
    
    # Sample free elements of A.
    # For equation i, regress (u_it-q_it) on -u_jt, j<i.
    u_std <- yt - B %*% xt
    u_neg <- -u_std
    a_sample <- numeric(K*(K-1)/2)
    if (K > 1L) {
      for (i in 2:K) {
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        ids <- id_start:id_end
        a_sub <- a_prior[ids]
        V_a_sub <- V_a_prior[ids, ids, drop = FALSE]
        y_A <- (u_std[i,] - q[i,]) / (sigma[i]*w_sqrt[i,])
        x_A <- u_neg[1:(i-1), , drop = FALSE] /
          matrix(rep(sigma[i]*w_sqrt[i,], each = i-1), nrow = i-1)
        a_sample[ids] <- sample_A_ele(
          ysub = y_A, xsub = x_A, a_sub = a_sub, V_a_sub = V_a_sub
        )
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A_post)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1
    
    # Sample gamma
    # r_it = (w_it-mu_i)*gamma_i + sqrt(w_it)*epsilon_it.
    r <- A %*% (yt - B %*% xt)
    kappa <- sweep(w, 1, mu, "-")
    for (i in seq_len(K)) {
      prior_var <- V_gamma_prior[i, i]
      if (!is.finite(prior_var) || prior_var <= 0) {
        stop("Gamma prior variances must be positive.")
      }
      precision_gamma <- 1/prior_var +
        sum(kappa[i,]^2 / (w[i,]*sigma2[i]))
      mean_gamma <- (gamma_prior[i]/prior_var +
                       sum(kappa[i,]*r[i,] / (w[i,]*sigma2[i]))) /
        precision_gamma
      gamma[i] <- rnorm(1, mean = mean_gamma, sd = sqrt(1/precision_gamma))
    }
    q <- kappa * matrix(gamma, nrow = K, ncol = t_max)
    
    # Sample sigma_i^2
    # Structural residuals must be corrected by the latest gamma draw.
    s <- (r - q) / w_sqrt
    sigma_post_a <- sigma0_T0 + rep(t_max, K)
    sigma_post_b <- sigma0_S0 + rowSums(s^2)
    for (i in seq_len(K)) {
      sigma2[i] <- rinvgamma(1,
                             shape = sigma_post_a[i]/2,
                             rate = sigma_post_b[i]/2)
    }
    sigma <- sqrt(sigma2)
    
    # Sample w_it
    for (i in seq_len(K)) {
      lambda <- -(nu[i] + 1)/2
      psi <- gamma[i]^2 / sigma2[i]
      for (tt in seq_len(t_max)) {
        chi <- nu[i] + (r[i,tt] + mu[i]*gamma[i])^2 / sigma2[i]
        if (gamma[i] == 0) {
          w[i,tt] <- rinvgamma(1, shape = -lambda, rate = chi/2)
        } else {
          w[i,tt] <- GIGrvg::rgig(1, lambda = lambda,
                                  chi = chi, psi = psi)
        }
      }
    }
    w_sqrt <- sqrt(w)
    
    # Sample nu
    for (i in seq_len(K)) {
      proposed_nu <- nu[i] + nu_rw_sd[i]*rnorm(1)
      if (is.finite(proposed_nu) && proposed_nu > 4) {
        log_nu_target <- function(value) {
          mu_value <- value/(value-2)
          adjusted_r <- r[i,] - (w[i,] - mu_value)*gamma[i]
          dgamma(value, shape = nu_gam_a, rate = nu_gam_b,
                 log = TRUE) +
            sum(invgamma::dinvgamma(w[i,], shape = value/2,
                                    rate = value/2, log = TRUE)) -
            0.5*sum(adjusted_r^2 / (w[i,]*sigma2[i]))
        }
        log_alpha <- log_nu_target(proposed_nu) - log_nu_target(nu[i])
        if (!is.na(log_alpha) && log(runif(1)) < min(0, log_alpha)) {
          nu[i] <- proposed_nu
        }
      }
    }
    mu <- nu/(nu-2)
    
    if (j > inits$burnin && (j - inits$burnin) %% inits$thin == 0) {
      save_idx <- (j - inits$burnin) %/% inits$thin
      mcmc[, save_idx] <- c(b_sample, a_sample, sigma, gamma, nu, as.vector(w))
      My_miss[, save_idx] <- ymiss_vec
      y_miss_mu[, save_idx] <- ym_mu_post
      y_miss_sd[, save_idx] <- ym_sd_post
    }
    
    if (j %% 1000 == 0) {
      cat("Iteration", j, "nu:", round(nu, 2),
          "gamma:", round(gamma, 2), "\n")
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  
  if (aggregation == "average" || aggregation == "triangular"){
    y_raw[,idq] <- NA
  }
  nameY <- matrix(paste("y", repcol(c(1:K),t_max), reprow(c(1:t_max),K), sep = "_"), nrow = K)
  row.names(My_miss) <- row.names(y_miss_mu) <-  row.names(y_miss_sd) <- nameY[is.na(t(y_raw))]
  
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}