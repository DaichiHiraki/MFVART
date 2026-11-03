#' Simulate data from MF VAR model without SV
#'
#' This function simulate data from a MF VAR model without SV.
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma eps_t}
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","Skew.Student","Skew.Student", "MT","Skew.MT","MST").
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param t_max The number of observations
#' @param b0 The coefficients of B matrix.
#' @param a0 The coefficients of A matrix.
#' @param h The log diag(Sigma) matrix.
#' @param nu The degree of freedom.
#' @param idq The indicator positions of the quarter variables
#' @param freq The frequency of observed quarter variables.
#' @param aggregation The method of aggregation: "identity", "average", "triangular"
#' @param seednum The default seed.
#' @param burn_in The discarded observations.
#' @return A list of simulated data with its specification.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.MFVAR.novol(dist="Gaussian")
#' y <- datagen$y
#' }
sim.MFVAR.novol <- function(dist, K = 5, p = 2, t_max = 1000,
                            b0 = 0.5, a0 = 0.5, h = 0, nu = 6, gamma = 0.5,
                            y0 = matrix(0, ncol = K, nrow = p), sigma_G = NULL,
                            idq = rep(0,K), freq = 3, aggregation = "identity", 
                            seednum = 0, burn_in = 0){
  if (!(dist %in% c("Gaussian","Student",
                    "MT","OT","OST") ))
    stop("dist is not implemented.")
  
  if (dist == "Gaussian") datagen <- sim.MFVAR.Gaussian.novol(K, p, t_max, b0, a0, h, y0, idq, freq, aggregation, seednum, burn_in)
  if (dist == "Student") datagen <- sim.MFVAR.Student.novol(K, p, t_max, b0, a0, h, y0, nu, idq, freq, aggregation, seednum, burn_in)
  
  if (dist == "MT") datagen <- sim.MFVAR.MT.novol(K, p, t_max, b0, a0, h, y0, nu, idq, freq, aggregation, seednum, burn_in)
  if (dist == "OT") datagen <- sim.MFVAR.OT.novol(K, p, t_max, b0, a0, h, y0, nu, idq, freq, aggregation, seednum, burn_in)
  if (dist == "OST") datagen <- sim.MFVAR.OST.novol(K, p, t_max, b0, a0, h, y0, nu, gamma, idq, freq, aggregation, seednum, burn_in)
  
  return(datagen)
}

#' @export
sim.MFVAR.Gaussian.novol <- function(K = 5, p = 2, t_max = 1000,
                                     b0 = 0.5, a0 = 0.5, h = 0,
                                     y0 = matrix(0, ncol = K, nrow = p),
                                     idq = rep(0,K), freq = 3, aggregation = "identity", 
                                     seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }
  
  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # No tail
  # No volatility
  
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  
  eps <- matrix(rnorm(t_max*K), ncol = K)
  
  Sigma <- solve(A0) %*% diag(as.numeric(exp(0.5*h)), nrow = K)
  Sigma2 <- Sigma %*% t(Sigma)
  
  y_var <- reprow(Sigma[lower.tri(Sigma, diag = T)], t_max) # Sigma[lower.tri(Sigma, diag = T)] <- y_var[1,]
  volatility <- reprow(diag(Sigma2), t_max)
  
  for (i in c(1:t_max)){
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }
  if (aggregation == "identity"){
    t_max = t_max - burn_in
    y = y_true = y_atrue = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  if (aggregation == "average" || aggregation == "triangular"){
    y_aggre <- as.matrix(Make_bandSparse(nrow(ystar), aggregation) %*% ystar[,idq])
    t_max = t_max - burn_in
    y_true = y_atrue = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y_atrue[,idq] <- y_aggre[(p+burn_in+1):(p+burn_in+t_max),]
    y = y_atrue
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  
  list(y = y, y_true = y_true, y_atrue = y_atrue,
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       aggregation = aggregation,
       dist = "Gaussian", SV = FALSE)
}
#' @export
sim.MFVAR.Student.novol <- function(K = 5, p = 2, t_max = 1000,
                                    b0 = 0.5, a0 = 0.5, h = 0,
                                    y0 = matrix(0, ncol = K, nrow = p),
                                    nu = 6, idq = rep(0,K), freq = 3, aggregation = "identity", 
                                    seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }
  
  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)
  
  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  
  eps <- matrix(rnorm(t_max*K), ncol = K)
  
  Sigma <- solve(A0) %*% diag(as.numeric(exp(0.5*h)), nrow = K)
  Sigma2 <- Sigma %*% t(Sigma)
  
  for (i in c(1:t_max)){
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    Sigma_t <- diag(w_sqrt_t[i], nrow = K) %*% Sigma
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  B0 %*% xt + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }
  
  if (aggregation == "identity"){
    t_max = t_max - burn_in
    y = y_true = y_atrue = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  if (aggregation == "average" || aggregation == "triangular"){
    y_aggre <- as.matrix(Make_bandSparse(nrow(ystar), aggregation) %*% ystar[,idq]) 
    t_max = t_max - burn_in
    y_true = y_atrue = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y_atrue[,idq] <- y_aggre[(p+burn_in+1):(p+burn_in+t_max),]
    y = y_atrue
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  list(y = y, y_true = y_true, y_atrue = y_atrue,
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max)],
       aggregation = aggregation,
       dist = "Student", SV = FALSE)
}

#' @export
sim.MFVAR.MT.novol <- function(K = 5, p = 2, t_max = 1000,
                               b0 = 0.5, a0 = 0.5, h = 0,
                               y0 = matrix(0, ncol = K, nrow = p),
                               nu = 6, idq = rep(0,K), freq = 3, aggregation = "identity", 
                               seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }
  
  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)
  
  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  
  eps <- matrix(rnorm(t_max*K), ncol = K)
  
  Sigma <- solve(A0) %*% diag(as.numeric(exp(0.5*h)), nrow = K)
  Sigma2 <- Sigma %*% t(Sigma)
  
  for (i in c(1:t_max)){
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    Sigma_t <- diag(w_sqrt_t[i,],nrow = K) %*% Sigma
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  B0 %*% xt + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }
  
  if (aggregation == "identity"){
    t_max = t_max - burn_in
    y = y_true = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  if (aggregation == "identity"){
    t_max = t_max - burn_in
    y = y_true = y_atrue = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  if (aggregation == "average" || aggregation == "triangular"){
    y_aggre <- as.matrix(Make_bandSparse(nrow(ystar), aggregation) %*% ystar[,idq]) 
    t_max = t_max - burn_in
    y_true = y_atrue = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y_atrue[,idq] <- y_aggre[(p+burn_in+1):(p+burn_in+t_max),]
    y = y_atrue
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  list(y = y, y_true = y_true, y_atrue = y_atrue,
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       aggregation = aggregation,
       dist = "MT", SV = FALSE)
}

#' @export
sim.MFVAR.OT.novol <- function(K = 5, p = 2, t_max = 1000,
                               b0 = 0.5, a0 = 0.5, h = 0,
                               y0 = matrix(0, ncol = K, nrow = p),
                               nu = 6, idq = rep(0,K), freq = 3, aggregation = "identity", 
                               seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }
  
  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)
  
  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  
  eps <- matrix(rnorm(t_max*K), ncol = K)
  
  inv_A0 <- solve(A0)
  
  for (i in c(1:t_max)){
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*h), nrow = K)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  B0 %*% xt + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }
  
  if (aggregation == "identity"){
    t_max = t_max - burn_in
    y = y_true = y_atrue = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  if (aggregation == "average" || aggregation == "triangular"){
    y_aggre <- as.matrix(Make_bandSparse(nrow(ystar), aggregation) %*% ystar[,idq]) 
    t_max = t_max - burn_in
    y_true = y_atrue = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max)
    y_atrue[,idq] <- y_aggre[(p+burn_in+1):(p+burn_in+t_max),]
    y = y_atrue
    y[-seq(1, t_max, by = freq) ,idq] = NA  
  }
  
  list(y = y, y_true = y_true, y_atrue = y_atrue,
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       aggregation = aggregation,
       dist = "OT", SV = FALSE)
}


