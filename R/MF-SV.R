#' Bayesian inference of MF VAR model with RW-SV
#'
#' Bayesian inference of MF VAR model with RW-SV
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) H_t eps_t}
#' @param y The input data as matrix T x K where T is the number of observations and K is the number of variables
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","OT").
#' @param y0 The number of observations
#' @param prior The prior specification of BMFVAR.
#' @param freq The frequency of observed quarter variables.
#' @return A coda object of the posterior samples.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.VAR.SV(dist="Gaussian", idq = c(4), p = 2, K = 4, t_max = 500)
#' y <- datagen$y
#' prior <- get_prior(y, p = 2, dist="Gaussian", SV = T)
#' inits <- get_init(prior)
#' Chain1 <- BMFVAR.SV(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
#' plot(Chain1)
#' }
BMFVAR.SV <- function(y, K, p, dist, y0 = NULL, prior = NULL, inits = NULL){
  if (!(dist %in% c("Gaussian","Student",
                    "MT","OT","cG","cT") ))
    stop("dist is not implemented.")
  if (prior$SV == TRUE){
    Start = Sys.time()
    if (dist == "Gaussian") Chain <- BMFVAR.Gaussian.SV(y, K, p, y0, prior, inits)
    if (dist == "Student") Chain <- BMFVAR.Student.SV(y, K, p, y0, prior, inits)
    if (dist == "MT") Chain <- BMFVAR.MT.SV(y, K, p, y0, prior, inits)
    if (dist == "OT") Chain <- BMFVAR.OT.SV(y, K, p, y0, prior, inits)
    if (dist == "cG") Chain <- BMFVAR.cG.SV(y, K, p, y0, prior, inits)
    if (dist == "cT") Chain <- BMFVAR.cT.SV(y, K, p, y0, prior, inits)
    
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
    return(out)
  } else {
    warning("prior$SV is TRUE")
  }
}

###########################################################################
#' @export
BMFVAR.Gaussian.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Gaussian", 
                       SV = TRUE, aggregation = "identity", idq = K)
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
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  h <- inits$h
  sigma_h <- inits$sigma_h
  
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
    K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, as.vector(exp(-h)) ) %*% G_miss
    b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, as.vector(exp(-h)) ) %*% ( c_const - G_obs %*% yobs_vec)
    
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
    # Sys.sleep(0.1)
    # plot(as.numeric(Make_bandSparse(t_max, aggregation="triangular") %*% ymiss_vec)[600:t_max], type ='l')

    # sample B
    yt = matrix(yt_vec, nrow = K)
    y = t(yt)
    xt <- makeRegressor( y, y0, t_max, K, p)
    wt <- as.vector(exp(-h/2)) # move up
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)
    
    # Sample vol
    ytilde <- A%*% (yt - B %*% xt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h
    
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
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1
    
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, diag(sigma_h), h0, as.numeric(h))
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
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}

###########################################################################
#' @export
BMFVAR.Student.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  yt_vec <- vec(t(y))
  yobs_vec <- na.exclude(yt_vec)
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Student", 
                       SV = TRUE, aggregation = "identity", idq = K)
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
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h
  
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
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + 1 + K + K + K*t_max + t_max,
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
    Wh <- rep(1/w_sample, each = K) * as.vector(exp(-h))
    
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
    
    wt <- as.vector(1/(exp(h/2)*w_sqrt))
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)
    
    # Sample vol
    ytilde <- A%*% ((yt - B %*% xt)/w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h
    
    
    # Sample A0
    u_std <- (yt - B %*%xt)/ w_sqrt # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
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
    w_sample <- rinvgamma( n = t_max, shape = (nu+K)*0.5, rate = 0.5*( nu + colSums((u^2)/exp(h)) ) )
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w_sample))
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
                        paste("nu"),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))
  
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}

###########################################################################
#' @export
BMFVAR.MT.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "MT", 
                       SV = TRUE, aggregation = "identity", idq = K)
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
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h
  
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
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max + K*t_max,
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
    c_B <- Make_cB(y0, B, K, p ,t_max)  # Change Make_cB
    W_mat <- bandSparse(K*t_max, k = -c(0:(K-1)), 
                        diag = lapply(X = c(1:K), FUN = function(i) 1/as.numeric(w_sqrt)[i:(K*t_max)] ),
                        symm=FALSE)
    IWA <- W_mat * kronecker(Diagonal(t_max), as(A,Class = "TsparseMatrix") )
    c_const <- IWA %*% c_B
    IAH <- IWA %*% H_B
    G_obs <-  IAH %*% S_obs
    G_miss <-  IAH %*% S_miss
    Wh <- as.numeric(exp(-h))
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
    
    wt <- as.vector(exp(-h/2))
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
    
    # Sample vol
    ytilde <- A%*% ((yt - B %*%xt)/ w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h
    
    
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
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
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
    u_proposal <- A %*% (yt - B %*%xt)
    a_target <- (nu*0.5 + 1*0.5) * 0.75
    b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sqrtvol^2)*0.75  # adjust by 0.75
    w_temp <- matrix( rinvgamma( n = K*t_max, shape = a_target, rate = b_target ), K, t_max ) # rinvgamma recycles arguments
    w_temp_sqrt <- sqrt(w_temp)
    
    A.w.u.s.prop <- ( A %*% ( u/w_temp_sqrt ) )/sqrtvol
    num_mh <- colSums(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T)) - #prior
      0.5*( colSums(log(w_temp)) + colSums(A.w.u.s.prop^2) ) - # posterior
      colSums( dinvgamma(w_temp, shape = a_target, rate = b_target, log = T)) # proposal
    A.w.u.s.curr <- ( A %*% ( u/w_sqrt ) )/sqrtvol
    denum_mh <- colSums(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) - #prior
      0.5*( colSums(log(w)) + colSums(A.w.u.s.curr^2) ) - # posterior
      colSums( dinvgamma(w, shape = a_target, rate = b_target, log = T)) # proposal
    acc <- (num_mh-denum_mh) > log(runif(t_max))
    w[,acc] <- w_temp[,acc]
    w_sqrt[,acc] <- w_temp_sqrt[,acc]
    acount_w[acc] <- acount_w[acc] + 1
    
    
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))
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
                        paste("nu",c(1:K), sep = ""),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))
  
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}

###########################################################################
#' @export
BMFVAR.OT.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OT", 
                       SV = TRUE, aggregation = "identity", idq = K)
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
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h
  
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
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max + K*t_max,
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
    Wh <- 1/vec(w) * vec(exp(-h))
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
    
    
    # Sample B
    yt = matrix(yt_vec, nrow = K)
    y = t(yt)
    xt <- makeRegressor( y, y0, t_max, K, p)
    
    # Sample B
    wt <- as.vector(1/(exp(h/2)*w_sqrt))
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)
    
    # Sample vol
    ytilde <- (A %*% (yt - B %*% xt)) / w_sqrt
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h
    
    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.numeric(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)
    
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
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,] / w_sqrt[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / reprow(sqrtvol[i,],i-1) / reprow(w_sqrt[i,], i-1), nrow = i-1),
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
    w <- matrix( rinvgamma( n = K*t_max, shape = nu*0.5 + 1*0.5, rate = nu*0.5 + 0.5 * (u^2) / exp(h) ),
                 K, t_max )
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))
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
                        paste("nu",c(1:K), sep = ""),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))
  
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}

###########################################################################
#' @export
BMFVAR.cG.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Gaussian", 
                       SV = TRUE, aggregation = "identity", idq = K)
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
  # prior h0_mean
  h0_mean <- mean(2 * log(prior$sigma))
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
  h_row <- colMeans(inits$h)*10
  h <- reprow(h_row, K)
  
  sigma_h <- inits$sigma_h[1,1]
  
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
  Hrho <-  Diagonal(t_max) - Matrix::bandSparse(t_max,t_max, -1,list(rep(1,t_max-1))) # T matrix
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + 1 + 1 + 1*t_max,
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
    K_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, as.vector(exp(-h)) ) %*% G_miss
    b_y_miss <- Matrix::t(G_miss) %*% Diagonal(t_max*K, as.vector(exp(-h)) ) %*% ( c_const - G_obs %*% yobs_vec)
    
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
    # Sys.sleep(0.1)
    # plot(as.numeric(Make_bandSparse(t_max, aggregation="triangular") %*% ymiss_vec)[600:t_max], type ='l')
    
    # sample B
    yt = matrix(yt_vec, nrow = K)
    y = t(yt)
    xt <- makeRegressor( y, y0, t_max, K, p)
    wt <- as.vector(exp(-h/2)) # move up
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)
    
    # Sample vol
    ytilde <- A%*% (yt - B %*% xt)
    s2 <- colSums(ytilde^2)
    # is_ForcedAccept <- ifelse(j < 20,TRUE,FALSE)
    
    {
      invSh <- Diagonal(t_max,as.numeric(1/sigma_h)) # Sigma_h
      HiSH <- Matrix::t(Hrho) %*% invSh %*% Hrho
      errh <- 1
      ht <- h_row
      
      while (errh > 1e-3) {
        eht <- exp(ht)
        sieht <- s2 / eht
        fh <- -K/2 + 0.5 * sieht
        Gh <- 0.5 * sieht
        Kh <- HiSH + Matrix::Diagonal(t_max, Gh)
        newht <- Matrix::solve(Kh, fh + Gh * ht)
        errh <- max(abs(newht - ht))
        ht <- newht
      }
      
      CKh <- Matrix::chol(Kh)
      
      # AR-step
      hstar <- ht
      logc <- -0.5 * t(hstar) %*% HiSH %*% hstar - K/2 * sum(hstar) - 0.5 * sum(exp(-hstar) * s2) + log(3)
      flag <- 0
      
      while (flag == 0) {
        hc <- ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max))
        alpARc <- as.numeric(-0.5 * t(hc) %*% HiSH %*% hc - K/2 * sum(hc) - 0.5 * sum(exp(-hc) * s2) +
          0.5 * t(hc - ht) %*% Kh %*% (hc - ht) - logc)
        if ( alpARc > log(runif(1))) {
          flag <- 1
        }
      }
      
      # MH-step
      alpAR <- as.numeric( -0.5 * t(h_row) %*% HiSH %*% h_row - K/2 * sum(h_row) - 0.5 * sum(exp(-h_row) * s2) + 
        0.5 * t(h_row - ht) %*% Kh %*% (h_row - ht) - logc)
      
      if (alpAR < 0) {
        alpMH <- 1
      } else if (alpARc < 0) {
        alpMH <- -alpAR
      } else {
        alpMH <- alpARc - alpAR
      }
      
      if (alpMH > log(runif(1)) ) {
      # if (alpMH > log(runif(1)) || is_ForcedAccept) {
        h_row <- hc
        h <- reprow(h_row, K)
        sqrtvol <-  exp(0.5*h)
        h0 <- h_row[1]
      }
      
    }
    sse_2 <- sum( (h_row[2:t_max] - h_row[1:(t_max-1)] )^2)
    sigma_h <- GIGrvg::rgig(n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                            psi = 1/prior$sigma_S0 )
    
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
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1
    
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma_h, h0, as.numeric(h_row))
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
                        paste("sigma_h",c(1), sep = ""),
                        paste("lh0",c(1), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1),t_max), reprow(c(1:t_max),1))
  )
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}

###########################################################################
#' @export
BMFVAR.cT.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  y_raw <- y
  yt_vec <- vec(t(y))
  yobs_vec <- na.exclude(yt_vec)
  
  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Student", 
                       SV = TRUE, aggregation = "identity", idq = K)
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
  # prior h0_mean
  h0_mean <- mean(2 * log(prior$sigma))
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
  B <- inits$B0
  h_row <- colMeans(inits$h)*10
  h <- reprow(h_row, K)
  
  sigma_h <- inits$sigma_h[1,1]
  
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
  Hrho <-  Diagonal(t_max) - Matrix::bandSparse(t_max,t_max, -1,list(rep(1,t_max-1))) # T matrix
  
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + 1 + 1 + 1 + 1*t_max + t_max,
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
    Wh <- rep(1/w_sample, each = K) * as.vector(exp(-h))
    
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
    # Sys.sleep(0.1)
    # plot(as.numeric(Make_bandSparse(t_max, aggregation="triangular") %*% ymiss_vec)[600:t_max], type ='l')
    
    # sample B
    yt = matrix(yt_vec, nrow = K)
    y = t(yt)
    xt <- makeRegressor( y, y0, t_max, K, p)
    wt <- as.vector(1/(exp(h/2)*w_sqrt))
    
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)
    
    # Sample vol
    ytilde <- A%*% ((yt - B %*% xt)/w_sqrt)
    s2 <- colSums(ytilde^2)
    # is_ForcedAccept <- ifelse(j < 20,TRUE,FALSE)
    
    {
      invSh <- Diagonal(t_max,as.numeric(1/sigma_h)) # Sigma_h
      HiSH <- Matrix::t(Hrho) %*% invSh %*% Hrho
      errh <- 1
      ht <- h_row
      
      while (errh > 1e-3) {
        eht <- exp(ht)
        sieht <- s2 / eht
        fh <- -K/2 + 0.5 * sieht
        Gh <- 0.5 * sieht
        Kh <- HiSH + Matrix::Diagonal(t_max, Gh)
        newht <- Matrix::solve(Kh, fh + Gh * ht)
        errh <- max(abs(newht - ht))
        ht <- newht
      }
      
      CKh <- Matrix::chol(Kh)
      
      # AR-step
      hstar <- ht
      logc <- -0.5 * t(hstar) %*% HiSH %*% hstar - K/2 * sum(hstar) - 0.5 * sum(exp(-hstar) * s2) + log(3)
      flag <- 0
      
      while (flag == 0) {
        hc <- ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max))
        alpARc <- as.numeric(-0.5 * t(hc) %*% HiSH %*% hc - K/2 * sum(hc) - 0.5 * sum(exp(-hc) * s2) +
                               0.5 * t(hc - ht) %*% Kh %*% (hc - ht) - logc)
        if ( alpARc > log(runif(1))) {
          flag <- 1
        }
      }
      
      # MH-step
      alpAR <- as.numeric( -0.5 * t(h_row) %*% HiSH %*% h_row - K/2 * sum(h_row) - 0.5 * sum(exp(-h_row) * s2) + 
                             0.5 * t(h_row - ht) %*% Kh %*% (h_row - ht) - logc)
      
      if (alpAR < 0) {
        alpMH <- 1
      } else if (alpARc < 0) {
        alpMH <- -alpAR
      } else {
        alpMH <- alpARc - alpAR
      }
      
      if (alpMH > log(runif(1)) ) {
        # if (alpMH > log(runif(1)) || is_ForcedAccept) {
        h_row <- hc
        h <- reprow(h_row, K)
        sqrtvol <-  exp(0.5*h)
        h0 <- h_row[1]
      }
      
    }
    sse_2 <- sum( (h_row[2:t_max] - h_row[1:(t_max-1)] )^2)
    sigma_h <- GIGrvg::rgig(n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                            psi = 1/prior$sigma_S0 )
    
    # Sample A0
    u_std <- (yt - B %*%xt)/ w_sqrt # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
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
    w_sample <- rinvgamma( n = t_max, shape = (nu+K)*0.5, rate = 0.5*( nu + colSums((u^2)/exp(h)) ) )
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, sigma_h, h0, as.numeric(h_row), as.numeric(w_sample))
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
                        paste("nu"),
                        paste("sigma_h",c(1), sep = ""),
                        paste("lh0",c(1), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1),t_max), reprow(c(1:t_max),1)),
                        paste("w",c(1:t_max), sep = "")
  )
  return(list( param = as.mcmc(t(mcmc)),
               y_miss = as.mcmc(t(My_miss)),
               y_miss_mu =  as.mcmc(t(y_miss_mu)),
               y_miss_sd =  as.mcmc(t(y_miss_sd))
  ))
}
