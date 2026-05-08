#' Bayesian inference of MF VAR model without SV
#'
#' Bayesian inference of MF VAR model without SV
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma eps_t}
#' @param y The input data as matrix T x K where T is the number of observations and K is the number of variables
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","OT").
#' @param y0 The number of observations
#' @param prior The prior specification of BMFVAR.
#' @param inits The initial values of BMFVAR.
#' @return A coda object of the posterior samples.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.MFVAR.novol(dist="Gaussian", idq = c(4), p = 2, K = 4, t_max = 500)
#' y <- datagen$y
#' prior <- get_prior(y, p = 2, dist="Gaussian", SV = F)
#' inits <- get_init(prior)
#' Chain1 <- BMFVAR.novol(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
#' plot(Chain1)
#' }
BMFVAR.novol <- function(y, K, p, dist, y0 = NULL, prior = NULL, inits = NULL){
  if (!(dist %in% c("Gaussian","Student",
                    "MT","OT") ))
    stop("dist is not implemented.")

  if (prior$SV == FALSE){
    Start = Sys.time()
    if (dist == "Gaussian") Chain <- BMFVAR.Gaussian.novol(y, K, p, y0, prior, inits)
    if (dist == "Student") Chain <- BMFVAR.Student.novol(y, K, p, y0, prior, inits)
    if (dist == "MT") Chain <- BMFVAR.MT.novol(y, K, p, y0, prior, inits)
    if (dist == "OT") Chain <- BMFVAR.OT.novol(y, K, p, y0, prior, inits)
    elapsedTime = Sys.time() - Start
    print(elapsedTime)
    out <- list(mcmc = Chain,
                y = y,
                y0 = y0,
                K = K,
                p = p,
                dist = dist,
                prior = prior,
                inits = inits,
                esttime = elapsedTime)
    class(out) <- c("MFVART")

  } else {
    warning("prior$SV is TRUE")
  }
  return(out)
}

###########################################################################
#' @export
BMFVAR.Gaussian.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Gaussian", 
                       SV = FALSE, aggregation = "identity", idq = K)
  }
  
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


  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(sigma, nrow = length(sigma))
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior
  S_select <- Make_SectMat(y)
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  n_cond <- length(y_raw[,idq]) - sum(is.na(y_raw[,idq]))
  S_obs <- S_select$S_obs
  S_miss <- S_select$S_miss
  K_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% M_aggregate + Diagonal(n_miss, ymiss_invprior)
  b_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% ( z_obs_vec)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  My_miss <- y_miss_mu <- y_miss_sd <- matrix(NA, nrow = n_miss,
                                              ncol = (samples - inits$burnin)%/% inits$thin)
  id_miss <- (1:(t_max*K))[is.na(t(y))] # Location of Y_miss
  id_miss_t <- rep(1:t_max, each = K)[id_miss] # Time of Y_miss
  id_miss_k <- rep(1:K, t_max)[id_miss]       # Variable k of Y_miss
  id_m_m <- (1:length(id_miss))[! id_miss_k  %in% idq ] # ID of ragged monthly Y_miss = S_select_temp$S_miss@i + 1
  id_m_q <- (1:length(id_miss))[id_miss_k  %in% idq ] # ID of ragged quarterly Y_miss = S_select_temp$S_obs@i + 1
  M_forecast <- kronecker(Make_bandSparse(t_max, aggregation)[is.na(y[,idq[1] ]),],
                          Diagonal(length(idq))) %*%  Matrix::t(S_select_temp$S_obs)
  
  for (j in c(1:samples)){
    
    # # Sample y missing
    H_B <- Make_HB(B, K, p ,t_max)
    c_B <- Make_cB(y0, B, K, p ,t_max)
    IA <- kronecker(Diagonal(t_max), as(A,Class = "TsparseMatrix") )
    c_const <- IA %*% c_B
    IAH <- IA %*% H_B
    G_obs <-  IAH %*% S_obs
    G_miss <-  IAH %*% S_miss
    K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, 1/sigma^2) %*% G_miss
    b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, 1/sigma^2) %*% ( c_const - G_obs %*% yobs_vec)
    
    if (aggregation == "identity"){
      # sample from N( K^-1 b, K^-1)
      U.chol <- Matrix::chol( K_y_miss )
      # ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_y_miss,
      #                                          upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      ym_mu_post <- backsolve(U.chol, backsolve(U.chol, b_y_miss, upper.tri = T, transpose = T )) 
      ym_Sigma_post <- backsolve(U.chol, backsolve(U.chol, diag(n_miss), upper.tri = T, transpose = T)) 
      ymiss_vec <- ym_mu_post + backsolve(U.chol, rnorm(n_miss) )
      ym_sd_post <- sqrt(diag(ym_Sigma_post))
    }
    if (aggregation == "average" || aggregation == "triangular"){
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
    
    # sample B
    yt = matrix(yt_vec, nrow = K)
    y = t(yt)
    xt <- makeRegressor( y, y0, t_max, K, p)

    A.tmp <- diag(1/sigma, K) %*% A
    y.tilde <- as.vector( A.tmp %*% yt )
    x.tilde <- kronecker( t(xt), A.tmp )
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                         + rnorm(K*m) )
    B <- matrix(b_sample,K,m)

    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% (yt - B %*%xt)
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
    u_std <- (yt - B %*%xt)
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i],
                                                     xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1

    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma)
      My_miss[, (j - inits$burnin) %/% inits$thin] <- ymiss_vec
      y_miss_mu[, (j - inits$burnin) %/% inits$thin] <- ym_mu_post
      y_miss_sd[, (j - inits$burnin) %/% inits$thin] <- ym_sd_post
    }
      
    if (j %% 1000 == 0) { cat(" Iteration ", j, " \n")}
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
                        paste("sigma",c(1:K), sep = ""))
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
         ))
}

###########################################################################
#' @export
BMFVAR.Student.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  yt_vec <- vec(t(y))
  yobs_vec <- na.exclude(yt_vec)
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Student", 
                       SV = FALSE, aggregation = "identity", idq = K)
  }
  
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
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  
  samples <- inits$samples
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(inits$sigma, nrow = K)
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)
  
  nu <- inits$nu
  logsigma_nu <- 0
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior
  S_select <- Make_SectMat(y)
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  n_cond <- length(y_raw[,idq]) - sum(is.na(y_raw[,idq]))
  S_obs <- S_select$S_obs
  S_miss <- S_select$S_miss
  K_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% M_aggregate + Diagonal(n_miss, ymiss_invprior)
  b_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% ( z_obs_vec)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + 1 + t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  My_miss <- y_miss_mu <- y_miss_sd <- matrix(NA, nrow = n_miss,
                                              ncol = (samples - inits$burnin)%/% inits$thin)
  id_miss <- (1:(t_max*K))[is.na(t(y))] # Location of Y_miss
  id_miss_t <- rep(1:t_max, each = K)[id_miss] # Time of Y_miss
  id_miss_k <- rep(1:K, t_max)[id_miss]       # Variable k of Y_miss
  id_m_m <- (1:length(id_miss))[! id_miss_k  %in% idq ] # ID of ragged monthly Y_miss = S_select_temp$S_miss@i + 1
  id_m_q <- (1:length(id_miss))[id_miss_k  %in% idq ] # ID of ragged quarterly Y_miss = S_select_temp$S_obs@i + 1
  M_forecast <- kronecker(Make_bandSparse(t_max, aggregation)[is.na(y[,idq[1] ]),],
                          Diagonal(length(idq))) %*%  Matrix::t(S_select_temp$S_obs)
  
  for (j in c(1:samples)){
    
    # # Sample y missing
    H_B <- Make_HB(B, K, p ,t_max)
    c_B <- Make_cB(y0, B, K, p ,t_max)
    IA <- kronecker(Diagonal(t_max), as(A,Class = "TsparseMatrix") )
    c_const <- IA %*% c_B
    IAH <- IA %*% H_B
    G_obs <-  IAH %*% S_obs
    G_miss <-  IAH %*% S_miss
    Wh <- rep(1/w_sample, each = K) * rep(1/sigma^2, t_max)
    K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% G_miss
    b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% ( c_const - G_obs %*% yobs_vec)
    
    if (aggregation == "identity"){
      # sample from N( K^-1 b, K^-1)
      U.chol <- Matrix::chol( K_y_miss )
      # ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_y_miss,
      #                                          upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      ym_mu_post <- backsolve(U.chol, backsolve(U.chol, b_y_miss, upper.tri = T, transpose = T )) 
      ym_Sigma_post <- backsolve(U.chol, backsolve(U.chol, diag(n_miss), upper.tri = T, transpose = T)) 
      ymiss_vec <- ym_mu_post + backsolve(U.chol, rnorm(n_miss) )
      ym_sd_post <- sqrt(diag(ym_Sigma_post))
    }
    if (aggregation == "average" || aggregation == "triangular"){
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
    
    # sample B
    yt = matrix(yt_vec, nrow = K)
    y = t(yt)
    xt <- makeRegressor( y, y0, t_max, K, p)
    
    A.tmp <- diag(1/sigma, K) %*% A
    wt <- as.vector(1/w_sqrt)
    y.tilde <- as.vector( A.tmp %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A.tmp ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)
    
    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% ((yt - B %*%xt)/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)
    
    # Sample A0
    u_std <- (yt - B %*%xt) / w_sqrt # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i] ,
                                                     xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1
    
    
    
    # Sample w
    u <- (yt - B %*%xt)
    shape_w <- nu*0.5 + K*0.5
    rate_w <- as.numeric(nu*0.5 + 0.5 * colSums((A %*% u / sigma)^2))
    w_sample <- rinvgamma( n=t_max, shape = shape_w, rate = rate_w )
    w <- reprow(w_sample, K)
    w_sqrt <- sqrt(w)
    
    
    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(1)
    if (nu_temp > 4 && nu_temp < 100){
      num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
        sum(dinvgamma(w_sample, shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))
      denum_mh = dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
        sum(dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T))
      alpha = num_mh - denum_mh;
      temp = log(runif(1));
      if (alpha > temp){
        nu = nu_temp
        acount_nu = acount_nu + 1
      }
      
    }
    
    if(j %% batchlength == 0 ){
      if (acount_nu > batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu + adaptamount(j %/% batchlength);
      }
      if (acount_nu < batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu - adaptamount(j %/% batchlength);
      }
      acount_nu = 0
    }
    
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w_sample))
      My_miss[, (j - inits$burnin) %/% inits$thin] <- ymiss_vec
      y_miss_mu[, (j - inits$burnin) %/% inits$thin] <- ym_mu_post
      y_miss_sd[, (j - inits$burnin) %/% inits$thin] <- ym_sd_post
    }
    
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", round(nu,2)," \n")
      acount_w <- rep(0,t_max)
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
                        paste("nu"),
                        paste("w",c(1:t_max), sep = ""))
  
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}


###########################################################################
#' @export
BMFVAR.MT.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "MT", 
                       SV = FALSE, aggregation = "identity", idq = K)
  }
  
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
    na_id <- is.na(yt_vec)
    obs_id <- !na_id
    
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
  
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  V_b_prior_inv <- solve(V_b_prior)
  
  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior
  S_select <- Make_SectMat(y)
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  n_cond <- length(y_raw[,idq]) - sum(is.na(y_raw[,idq]))
  S_obs <- S_select$S_obs
  S_miss <- S_select$S_miss
  K_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% M_aggregate + Diagonal(n_miss, ymiss_invprior)
  b_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% ( z_obs_vec)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  My_miss <- y_miss_mu <- y_miss_sd <- matrix(NA, nrow = n_miss,
                                              ncol = (samples - inits$burnin)%/% inits$thin)
  id_miss <- (1:(t_max*K))[is.na(t(y))] # Location of Y_miss
  id_miss_t <- rep(1:t_max, each = K)[id_miss] # Time of Y_miss
  id_miss_k <- rep(1:K, t_max)[id_miss]       # Variable k of Y_miss
  id_m_m <- (1:length(id_miss))[! id_miss_k  %in% idq ] # ID of ragged monthly Y_miss = S_select_temp$S_miss@i + 1
  id_m_q <- (1:length(id_miss))[id_miss_k  %in% idq ] # ID of ragged quarterly Y_miss = S_select_temp$S_obs@i + 1
  M_forecast <- kronecker(Make_bandSparse(t_max, aggregation)[is.na(y[,idq[1] ]),],
                          Diagonal(length(idq))) %*%  Matrix::t(S_select_temp$S_obs)
  for (j in c(1:samples)){
    # # Sample y missing
    H_B <- Make_HB(B, K, p ,t_max)
    
    # Another calculation
    c_B <- Make_cB(y0, B, K, p ,t_max)  # Change Make_cB
    W_mat <- bandSparse(K*t_max, k = -c(0:(K-1)), 
                        diag = lapply(X = c(1:K), FUN = function(i) 1/as.numeric(w_sqrt)[i:(K*t_max)] ),
                        symm=FALSE)
    IWA <- W_mat * kronecker(Diagonal(t_max), as(A,Class = "TsparseMatrix") )
    c_const <- IWA %*% c_B
    IAH <- IWA %*% H_B
    G_obs <-  IAH %*% S_obs
    G_miss <-  IAH %*% S_miss
    Wh <- rep(1/sigma^2, t_max)
    K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% G_miss
    b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% ( c_const - G_obs %*% yobs_vec)

    if (aggregation == "identity"){
      # sample from N( K^-1 b, K^-1)
      U.chol <- Matrix::chol( K_y_miss )
      # ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_y_miss,
      #                                          upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      ym_mu_post <- backsolve(U.chol, backsolve(U.chol, b_y_miss, upper.tri = T, transpose = T )) 
      ym_Sigma_post <- backsolve(U.chol, backsolve(U.chol, diag(n_miss), upper.tri = T, transpose = T)) 
      ymiss_vec <- ym_mu_post + backsolve(U.chol, rnorm(n_miss) )
      ym_sd_post <- sqrt(diag(ym_Sigma_post))
    }
    if (aggregation == "average" || aggregation == "triangular"){
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
    
    # sample B
    yt = matrix(yt_vec, nrow = K)
    y = t(yt)
    xt <- makeRegressor( y, y0, t_max, K, p)
    
    wt <- rep(1/sigma, t_max)
    y.tilde <- as.vector( A %*% ( w^(-0.5) * yt ) ) * wt
    # make and stack X matrices
    x.tilde <- kronecker( t(xt), A )
    tmp <- kronecker( t(w^(-1/2)), matrix(1,K,1) ) * wt # W^(-1/2)
    ii <- 0
    for (i in 1:m) {
      x.tilde[,(ii+1):(ii+K)] <- x.tilde[,(ii+1):(ii+K)] * tmp
      ii <- ii+K
    }
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)
    
    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% ((yt - B %*%xt)/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    # sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)
    
    # Sample A0
    u_std <- (yt - B %*%xt) / w_sqrt # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i],
                                                     xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1
    
    Sigma <- solve(A) %*% diag(sigma, nrow = length(sigma))
    Sigma2 <- Sigma %*% t(Sigma)
    Sigma2_inv <- solve(Sigma2)
    Sigma2_inv <- (Sigma2_inv + t(Sigma2_inv))*0.5
    
    # Sample w
    u <- (yt - B %*%xt)
    u_proposal <- A %*% (yt - B %*%xt)
    
    a_target <- (nu*0.5 + 1*0.5) * 0.75 # adjust by 0.75
    b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sigma^2) * 0.75 # adjust by 0.75
    w_temp <- matrix(rinvgamma( n = K*t_max, shape = a_target, rate = b_target), nrow = K)
    w_temp_sqrt <- sqrt(w_temp)
    prior_num <- colSums(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T))
    prior_denum <- colSums(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T))
    proposal_num <- colSums(dinvgamma(w_temp, shape = a_target, rate = b_target, log = T))
    proposal_denum <- colSums(dinvgamma(w, shape = a_target, rate = b_target, log = T))
    
    A.w.u.s.prop <- ( A %*% ( u/w_temp_sqrt ) )*( 1/ sigma )
    num_mh <- prior_num - proposal_num - 0.5 * colSums(log(w_temp)) - 0.5 * colSums(A.w.u.s.prop^2)
    A.w.u.s.curr <- ( A %*% ( u/w_sqrt ) )*( 1/ sigma )
    denum_mh <- prior_denum - proposal_denum - 0.5 * colSums(log(w)) - 0.5 * colSums(A.w.u.s.curr^2)
    
    alpha = num_mh - denum_mh
    temp = log(runif(t_max))
    acre = alpha > temp
    w[,acre] <- w_temp[,acre]
    w_sqrt[,acre] <- w_temp_sqrt[,acre]
    #    w_sqrt_inv <- 1/w_sqrt
    acount_w[acre] <- acount_w[acre] + 1
    
    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(K)
    for (k in c(1:K)){
      if (nu_temp[k] > 4 && nu_temp[k] < 100){
        num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
        denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
        alpha = num_mh - denum_mh;
        temp = log(runif(1));
        if (alpha > temp){
          nu[k] = nu_temp[k]
          acount_nu[k] = acount_nu[k] + 1
        }
        
      }
    }
    
    if(j %% batchlength == 0 ){
      for (jj in c(1:K)) {
        if (acount_nu[jj] > batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
        }
        if (acount_nu[jj] < batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
        }
        acount_nu[jj] = 0
      }
    }
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w))
      My_miss[, (j - inits$burnin) %/% inits$thin] <- ymiss_vec
      y_miss_mu[, (j - inits$burnin) %/% inits$thin] <- ym_mu_post
      y_miss_sd[, (j - inits$burnin) %/% inits$thin] <- ym_sd_post
    }
    
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " \n")
      acount_w <- rep(0,t_max)
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
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )
  
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}


###########################################################################
#' @export
BMFVAR.OT.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OT", 
                       SV = FALSE, aggregation = "identity", idq = K)
  }
  
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
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  w_sqrt_inv <- 1/w_sqrt
  
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior
  S_select <- Make_SectMat(y)
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  n_cond <- length(y_raw[,idq]) - sum(is.na(y_raw[,idq]))
  S_obs <- S_select$S_obs
  S_miss <- S_select$S_miss
  K_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% M_aggregate + Diagonal(n_miss, ymiss_invprior)
  b_prior_y_miss <- Matrix::t(M_aggregate) %*% Diagonal(n_cond, 1/sdeviation) %*% ( z_obs_vec)
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  My_miss <- y_miss_mu <- y_miss_sd <- matrix(NA, nrow = n_miss,
                                              ncol = (samples - inits$burnin)%/% inits$thin)
  id_miss <- (1:(t_max*K))[is.na(t(y))] # Location of Y_miss
  id_miss_t <- rep(1:t_max, each = K)[id_miss] # Time of Y_miss
  id_miss_k <- rep(1:K, t_max)[id_miss]       # Variable k of Y_miss
  id_m_m <- (1:length(id_miss))[! id_miss_k  %in% idq ] # ID of ragged monthly Y_miss = S_select_temp$S_miss@i + 1
  id_m_q <- (1:length(id_miss))[id_miss_k  %in% idq ] # ID of ragged quarterly Y_miss = S_select_temp$S_obs@i + 1
  M_forecast <- kronecker(Make_bandSparse(t_max, aggregation)[is.na(y[,idq[1] ]),],
                          Diagonal(length(idq))) %*%  Matrix::t(S_select_temp$S_obs)
  
  for (j in c(1:samples)){
    
    # # Sample y missing
    H_B <- Make_HB(B, K, p ,t_max)
    c_B <- Make_cB(y0, B, K, p ,t_max)
    IA <- kronecker(Diagonal(t_max), as(A,Class = "TsparseMatrix") )
    c_const <- IA %*% c_B
    IAH <- IA %*% H_B
    G_obs <-  IAH %*% S_obs
    G_miss <-  IAH %*% S_miss
    Wh <- 1/vec(w) * rep(1/sigma^2, t_max)
    K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% G_miss
    b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, Wh) %*% ( c_const - G_obs %*% yobs_vec)
    
    if (aggregation == "identity"){
      # sample from N( K^-1 b, K^-1)
      U.chol <- Matrix::chol( K_y_miss )
      # ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_y_miss,
      #                                          upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
      ym_mu_post <- backsolve(U.chol, backsolve(U.chol, b_y_miss, upper.tri = T, transpose = T )) 
      ym_Sigma_post <- backsolve(U.chol, backsolve(U.chol, diag(n_miss), upper.tri = T, transpose = T)) 
      ymiss_vec <- ym_mu_post + backsolve(U.chol, rnorm(n_miss) )
      ym_sd_post <- sqrt(diag(ym_Sigma_post))
    }
    if (aggregation == "average" || aggregation == "triangular"){
      K_bar_y_miss <- K_y_miss + K_prior_y_miss
      b_bar_y_miss <- b_y_miss + b_prior_y_miss
      U.chol <- Matrix::chol( K_bar_y_miss )
      ymiss_vec <- backsolve(U.chol, backsolve(U.chol, b_bar_y_miss,
                                               upper.tri = T, transpose = T ) + rnorm(n_miss) ) 
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
    yt = matrix(yt_vec, nrow = K)
    y = t(yt)
    xt <- makeRegressor( y, y0, t_max, K, p)
    
    wt <- rep( 1/sigma, t_max ) / as.numeric(w_sqrt)
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)
    
    # Sample sigma
    u_ort <- (A %*% (yt - B %*%xt))/ w_sqrt    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)
    
    # Sample A0
    u_std <- (yt - B %*%xt) # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i]/ w_sqrt[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / sigma[i] / reprow(w_sqrt[i,], i-1), nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1
    
    # Sample w
    u <-  A %*% (yt - B %*%xt)
    w <- matrix( rinvgamma( n = K*t_max, shape = (nu+1)*0.5, rate = 0.5*( nu + (u^2)/sigma2 ) ), K, t_max )
    w_sqrt <- sqrt(w)
    
    
    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(K)
    for (k in c(1:K)){
      if (nu_temp[k] > 4 && nu_temp[k] < 100){
        num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
        denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
        alpha = num_mh - denum_mh;
        temp = log(runif(1));
        if (alpha > temp){
          nu[k] = nu_temp[k]
          acount_nu[k] = acount_nu[k] + 1
        }
        
      }
    }
    
    if(j %% batchlength == 0 ){
      for (jj in c(1:K)) {
        if (acount_nu[jj] > batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
        }
        if (acount_nu[jj] < batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
        }
        acount_nu[jj] = 0
      }
    }
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w))
      My_miss[, (j - inits$burnin) %/% inits$thin] <- ymiss_vec
      y_miss_mu[, (j - inits$burnin) %/% inits$thin] <- ym_mu_post
      y_miss_sd[, (j - inits$burnin) %/% inits$thin] <- ym_sd_post
    }
    
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " \n")
      acount_w <- rep(0,t_max)
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
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}


