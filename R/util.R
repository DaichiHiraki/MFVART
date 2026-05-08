##########################################################################
# Adatptive metropolis hasting #
##########################################################################
TARGACCEPT = 0.3
batchlength = 10

#' @export
adaptamount <- function(iteration){
  return( min( 0.01, 1.0 / sqrt(iteration) ) );
}

# This function is borrowed from https://github.com/FK83/bvarsv
#' @export
getmix <- function(){
  # 7 components
  q <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)      # probabilities
  m <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) # means
  u2 <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)    #variances

  # # 10 components
  # q <- c(0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575, 0.00115)      # probabilities
  # m <- c(3.13417, 2.55484, 1.94244, 1.23006, 0.35567, -0.76538, -2.26048, -4.34506, -7.47644, -13.44260)
  # u2 <- c(0.11265, 0.17788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342)    #variances
  # # m <- c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000) # means in the Omori paper
  return(list(q=q,m=m,u2=u2))
}
##########################################################################

##########################################################################
# Utility functions  #
##########################################################################
#' @export
makeRegressor <- function(y, y0, t_max, K, p){
  # xt <- matrix(NA, nrow = K * p, ncol = t_max)
  yt = t(y)

  if (is.null(y0)){
    y0 <- matrix(0, ncol = K, nrow = p)
  } else {
    y0 <- tail(y0, p)
  }

  if (p > 0){
    # y0 <- rbind(y0, y[1:p,])
    # for (i in 1:p) {
    #   xt[, i] <- as.vector(t(y0)[, (p+i-1):i])
    # }
    # for (i in (p+1):t_max) {
    #   xt[, i] <- as.vector(yt[, (i-1):(i-p)])
    # }
    
    y0 <- t(rbind(y0, y))
    # for (i in 1:t_max) {
    #   xt[, i] <- as.vector(y0[, (p+i-1):i])
    # }
    
    id_col <- t(sapply(p:1, FUN = function(i) i:(i+t_max-1 )))
    xt <- matrix(y0[, id_col], ncol = t_max)
    xt <- rbind( rep(1,t_max), xt)
  } else {
    xt <- matrix(1, nrow = 1, ncol = t_max)
  }

  return(xt)
}

#' @export
Vec_to_Mat <- function(b, K, p){
  return(matrix(b, nrow = K))
}

#' @export
reprow<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#' @export
repcol<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#' @export
sample_A_ele <- function(ysub, xsub, a_sub, V_a_sub){
  t_max = length(ysub)
  n = nrow(xsub)

  a_post = rep(0, n)
  V_a_sub_inv = solve(V_a_sub)

  V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
  V_a_post <- solve(V_a_post_inv)
  a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)
  return(a_post  + t(chol(V_a_post)) %*% rnorm(n))
}

#' @export
sample_h_ele <- function(ytilde, sigma_h = 0.0001*diag(K), h0_mean = rep(0,K),
                         h = matrix(0, nrow = t_max, ncol = K), K, t_max, prior){
  tmp <- getmix()
  q <- tmp$q
  m_mean <- tmp$m
  u2 <- tmp$u2
  # {
  #   Zs <- matrix(1,t_max,1) %x% diag(K)
  #   sigma_prmean <- h0_mean # mean h_0
  #   sigma_prvar <- 4*diag(K)   # variance h_0
  #
  #
  #   aux <- sigmahelper4(t(ytilde^2), q, m_mean, u2, h, Zs, sigma_h, sigma_prmean, sigma_prvar)
  #   h <- aux$Sigtdraw
  #   h0 <- as.numeric(aux$h0)
  # }

  {
    h0 <- rep(0,K)
    cond_var_sigma <- rep(0,K)
    cond_mean_sigma <- rep(0,K)
    Zs <- matrix(1,t_max,1) %x% diag(1)
    for (i in c(1:K)){
      sigma_prmean <- h0_mean[i] # mean h_0
      sigma_prvar <- matrix(4)   # variance h_0
      aux <- sigmahelper4(t(ytilde[ i,, drop =FALSE]^2), q, m_mean, u2, h[ i,, drop =FALSE], Zs, matrix(sigma_h[i,i]), sigma_prmean, sigma_prvar)
      h[i,] <- aux$Sigtdraw
      h0[i] <- as.numeric(aux$h0)
      # vart <- aux$vart # Variance of components
      # yss1 <- aux$yss1 # demean of components
      # h_tilde <- (h[i,1:t_max] - h0[i])/sqrt(sigma_h[i,i])
      # cond_var_sigma[i] <- 1/( 1 + sum(h_tilde^2/vart))
      # cond_mean_sigma[i] <- cond_var_sigma[i] * ( 0 + sum(h_tilde * yss1 / vart))
    }
  }

  # sqrtvol <- aux$sigt
  # [TODO] fix this
  # sse_2 <- apply( (h[,2:t_max] - h[,1:(t_max-1)])^2, MARGIN = 1, FUN = sum)
  if (K>1) {
    sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
  } else {
    sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
  }

  # # Normal prior
  # # Equation 9 in https://doi.org/10.1016/j.csda.2013.01.002
  # sigma_post_a <- rep(t_max,K) # prior of sigma_h Gamma(1,0.0001)
  # sigma_post_b <- sse_2 # prior of sigma_h
  #
  # for (i in c(1:K)){
  #   sigma_new <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
  #   alpha = (sigma_h[i,i] - sigma_new) / 2 / prior$sigma_S0 + 0.5 * (log(sigma_new) - log(sigma_h[i,i])) # B_sigma = 1
  #   temp = log(runif(1))
  #   if (alpha > temp){
  #     sigma_h[i,i] <- sigma_new
  #   }
  #   #log_sigma_den[]
  # }

  sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                       psi = 1/prior$sigma_S0 ) , nrow = K)



  # # Invgamma conjugate prior
  # sigma_post_a <- 1 + rep(t_max,K) # prior of sigma_h Gamma(1,0.0001)
  # sigma_post_b <- 0.01 + sse_2 # prior of sigma_h
  #
  #   for (i in c(1:K)){
  #     sigma_h[i,i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
  #   }

  aux <- list(sigma_h = sigma_h,
              h0 = h0,
              Sigtdraw = h,
              sigt = exp(0.5*h))
  return(aux)
}


#' @export
sample_h_Chan <- function(ytilde, sigma_h = 0.0001*diag(K), h0_mean = rep(0,K),
                         h = matrix(0, nrow = K, ncol = t_max), K, t_max, prior){
  pi <- c(0.0073, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.2575)
  mi <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) - 1.2704  # already adjusted
  sigi <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)
  sqrtsigi <- sqrt(sigi)
  
  Hphi <- sparseMatrix(i = c(1:t_max, 2:t_max), 
                       j = c(1:t_max, 1:(t_max-1)), 
                       x = c(rep(1, t_max), rep(-1, t_max-1)), 
                       dims = c(t_max, t_max))
  
  # Zs <- matrix(1,t_max,1) %x% diag(1)
  
  for (ii in 1:K) {
    sigma_prmean <- h0_mean[ii] # mean h_0
    sigma_prvar <- matrix(4)   # variance h_0
    
    ystar <- log(ytilde[ii,]^2 + 0.0001)
    
    # Sample S from a 7-point discrete distribution
    temprand <- runif(t_max)
    q <- matrix(rep(pi, each = t_max), nrow = t_max) * dnorm(matrix(rep(ystar, times = 7), nrow = t_max),
                                                     matrix(rep(h[ii,], times = 7), nrow = t_max) + 
                                                       matrix(rep(mi, each = t_max), nrow = t_max),
                                                     matrix(rep(sqrtsigi, each = t_max), nrow = t_max))
    q <- q / rowSums(q)
    S <- 7 - rowSums(matrix(rep(temprand, times = 7), nrow = t_max) < t(apply(q, 1, cumsum))) + 1
    
    ## sample h using the precision-based algorithm
    # y* = h + d + epsilon*, epsilon* ~ N(0,Sig_y*),
    # Hphi h = alpha-tilde + zeta, zeta ~ N(0,Sig_h),
    # where 
    # d_t = E(epsilon*_t), Sig_y* = diag(var(epsilon*_1), ..., epsilon*_T),
    # Sig_h = diag(sig2/(1-phi^2), sig2, ..., sig2)

    sig2 <- sigma_h[ii,ii]
    invSigh <- sparseMatrix(i = 1:t_max, 
                            j = 1:t_max, 
                            x = rep(1/sig2, t_max))
    
    d <- mi[S]
    invSigystar <- Diagonal(t_max, 1/sigi[S])
    
    alpha <- c(sigma_prmean, rep(0, t_max-1))
    
    Kh <- Matrix::t(Hphi) %*% invSigh %*% Hphi + invSigystar
    dh <- Matrix::t(Hphi) %*% invSigh %*% alpha + invSigystar %*% (ystar - d)
    Ch <- chol(Kh)
    

    h[ii,] <- backsolve( Ch,
                         backsolve( Ch, dh,
                                    upper.tri = T, transpose = T )
                         + rnorm(t_max) )
  }
  h0 <- rep(0, K)
  

  
  if (K>1) {
    sse_2 <- apply( (h[,2:t_max] - h[,1:(t_max-1)] )^2, MARGIN = 1, FUN = sum)
  } else {
    sse_2 <- sum( (h[,2:t_max] - h[,1:(t_max-1)] )^2)
  }

  sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max -1 - 1)*0.5, chi = sse_2,
                          psi = 1/prior$sigma_S0 ) , nrow = K)

  
  aux <- list(sigma_h = sigma_h,
              h0 = h0, # here do not sample h0
              Sigtdraw = h,
              sigt = exp(0.5*h))
  return(aux)
}
#' @export
sample_h_mod <- function(ytilde, sigma_h = 0.0001*diag(K), h0_mean = rep(0,K),
                         h = matrix(0, nrow = t_max, ncol = K), K, t_max, prior){
  tmp <- getmix()
  q <- tmp$q
  m_mean <- tmp$m
  u2 <- tmp$u2

    h0 <- rep(0,K)
    h0mean <- rep(0,K)
    h0var <- rep(0,K)

    cond_var_sigma <- rep(0,K)
    cond_mean_sigma <- rep(0,K)
    Zs <- matrix(1,t_max,1) %x% diag(1)
    for (i in c(1:K)){
      sigma_prmean <- h0_mean[i] # mean h_0
      sigma_prvar <- matrix(4)   # variance h_0
      aux <- sigmahelper4(t(ytilde[ i,, drop =FALSE]^2), q, m_mean, u2, h[ i,, drop =FALSE], Zs, matrix(sigma_h[i,i]), sigma_prmean, sigma_prvar)
      h[i,] <- aux$Sigtdraw
      h0[i] <- as.numeric(aux$h0)
      h0mean[i] <- as.numeric(aux$h0mean)
      h0var[i] <- as.numeric(aux$h0var)

    }

  aux <- list(h0 = h0,
              Sigtdraw = h,
              sigt = exp(0.5*h),
              h0mean = h0mean,
              h0var = h0var)
  return(aux)
}
##########################################################################
# Plot functions  #
##########################################################################

#' @export
plot.MFVART <- function(MFVARTobj, element = NULL, namex = "", namey = "", namez = ""){
  if (is.null(element)) {
    y_miss <- MFVARTobj$mcmc$y_miss
    num_sample = nrow(y_miss)
    num_latent = ncol(y_miss)
    aggregation = MFVARTobj$prior$aggregation
    mcmc_tmp <- matrix(NA, nrow = num_sample, ncol = num_latent)
    for (i in c(1:num_sample)){
      mcmc_tmp[i,] <- as.numeric(Make_bandSparse(num_latent, aggregation=aggregation) %*% y_miss[i,])
    }
      pdat <- t(apply(mcmc_tmp, 2, function(z) quantile(z, c(0.025, 0.25, 0.5, 0.75, 0.975))))
      xax <- tail(Data$Time, num_latent)
      # matplot(x = xax, y = pdat*1.25, type = "n", ylab = "", xlab = "Time")
      # polygon(c(xax, rev(xax)), c(pdat[,5], rev(pdat[,4])), col = "grey60", border = NA)
      # polygon(c(xax, rev(xax)), c(pdat[,4], rev(pdat[,3])), col = "grey30", border = NA)
      # polygon(c(xax, rev(xax)), c(pdat[,3], rev(pdat[,2])), col = "grey30", border = NA)
      # polygon(c(xax, rev(xax)), c(pdat[,2], rev(pdat[,1])), col = "grey60", border = NA)
      # lines(x = xax, y = pdat[,3], type = "l", col = 1, lwd = 1)
      # abline(h = 0, lty = 2)
      
      years <- tail(Data$Time, num_latent)
      ObsValue <- tail(Data$GDP_growth, num_latent)
      df <- data.frame(Date = years, Value = pdat[,3], ObsValue = ObsValue)
      p1 <- ggplot(df, aes(x = Date, y = Value)) +
        geom_line(color = "black") +
        geom_point(aes(y = ObsValue), color = "red", size = 1) + # Scatter plot
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = namex, y = namey, title = namez) +
        ylim(c(-10,10)) + 
        theme_minimal() 
      p1
      return(p1)
      

  } else {
    
    y_miss <- get_post(MFVARTobj$mcmc$y_miss, element = element)
    
    num_sample = nrow(y_miss)
    num_latent = ncol(y_miss)
    
    mcmc_tmp <- matrix(NA, nrow = num_sample, ncol = num_latent)
    for (i in c(1:num_sample)){
      mcmc_tmp[i,] <- as.numeric(Make_bandSparse(num_latent, aggregation="triangular") %*% y_miss[i,])
    }
    pdat <- t(apply(mcmc_tmp, 2, function(z) quantile(z, c(0.025, 0.25, 0.5, 0.75, 0.975))))
    xax <- tail(Data$Time, num_latent)
    matplot(x = xax, y = pdat*1.25, type = "n", ylab = "", xlab = "Time")
    polygon(c(xax, rev(xax)), c(pdat[,5], rev(pdat[,4])), col = "grey60", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,4], rev(pdat[,3])), col = "grey30", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,3], rev(pdat[,2])), col = "grey30", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,2], rev(pdat[,1])), col = "grey60", border = NA)
    lines(x = xax, y = pdat[,3], type = "l", col = 1, lwd = 1)
    abline(h = 0, lty = 2)
  }
}


##########################################################################
# Get posterior functions  #
##########################################################################
#' Get Posterior samples of BVAR model
#'
#' This function returns Posterior samples of BVAR-SV-fatTail model.
#' @param obj The Chain/mcmc obtained by MFVARTobj$mcmc
#' @param element The name of parameters.
#' @return The Posterior samples of BVAR-SV-fatTail model in mcmc object
#' @export
#' @examples
#' \dontrun{
#' B_samples <- get_post(Chain, element = "B")
#' }
get_post <- function(obj, element = NULL, ...) {
  UseMethod("get_post", obj)
}

#' @export
get_post.MFVART <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj$mcmc)
  } else {
    mcmc_name <- substr(colnames(obj$mcmc$param), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj$mcmc$param[,mcmc_id, drop = FALSE])
  }
}

#' @export
get_post.numeric <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(names(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj[mcmc_id])
  }
}

#' @export
get_post.matrix <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(colnames(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(as.matrix(obj[,mcmc_id], nrow = nrow(obj)))
  }
}

#' @export
get_post.mcmc <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(colnames(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj[,mcmc_id, drop = FALSE])
  }
}

#' @export
Make_SectMat <- function(y){
  n = ncol (y)
  t_max = nrow(y)
  
  yt <- t(y)
  yt_vec <- vec(yt)
  
  n_obs <- sum(!is.na(y))
  n_miss <- sum(is.na(y))
  yobs_vec <- na.exclude(yt_vec)
  ymiss_vec <- matrix(NA, nrow = n_miss, ncol = 1)
  
  nT <- n_obs+n_miss
  i_obs <- j_obs <- rep(0,n_obs)
  i_miss <- j_miss <- rep(0,n_miss)
  c_obs <- c_miss <- 0
  
  for (t in c(1:t_max))
    for (i in c(1:n)){
      if (is.na(yt[i,t])){
        c_miss <- c_miss + 1 
        i_miss[c_miss] <- i + n*(t-1)
        j_miss[c_miss] <- c_miss
      } else{
        c_obs <- c_obs + 1
        i_obs[c_obs] <- i + n*(t-1)
        j_obs[c_obs] <- c_obs
      }
    }
  S_obs <- sparseMatrix(i = i_obs, j = j_obs, x = 1, dims =c(nT,n_obs), repr = "T")
  S_miss <- sparseMatrix(i_miss, j_miss, x = 1, dims =c(nT,n_miss), repr = "T")
  # cbind(yt_vec, S_obs%*% yobs_vec + S_miss%*% ymiss_vec)
  return(list(S_obs = S_obs,
              S_miss = S_miss))
}

#' @export
Make_HB <- function(B,K,p,t_max){
  b_intercept = B[,1]
  m =  B[,2:ncol(B)]
  B_t = - t(matrix(t(m), nrow=nrow(m))[, c(matrix(1:ncol(m), nrow(m), byrow=T)) ])
  
  nT = K*t_max
  
  #########################
  # H_B1 <- sparseMatrix(i = c(1:nT), j = c(1:nT), x = 1, dims =c(nT,nT), repr = "T")
  # for (t in c(1:(t_max-1))){
  #   p_max <- min(p, t_max-t)
  #   row_range <- t*K+c(1:(K*p_max))
  #   H_B1[row_range,(t-1)*K+c(1:K)] <- head(B_t, length(row_range))
  # }
  #########################
  #H_B <- sparseMatrix(i = c(1:nT), j = c(1:nT), x = 1, dims =c(nT,nT), repr = "T")
  i_stack <- c(1:nT)
  j_stack <- c(1:nT)
  x_stack <- rep(1, nT)
  for (i in c(1:p)){
    # tmp <- as(bdiag(replicate(t_max-i,B_t[(i-1)*K + c(1:K),1:K],simplify=FALSE)),Class = "TsparseMatrix")
    # # tmp2 <- sparseMatrix(i = tmp@i+K*i+1, j = tmp@j+1, x = tmp@x, dims = c(nT,nT), repr = "T")
    # i_stack <- c(i_stack, tmp@i+K*i+1)
    # j_stack <- c(j_stack, tmp@j+1)
    # x_stack <- c(x_stack, tmp@x)
    
    i_post <- matrix(rep( (i*K+1):(t_max*K), each = K), nrow = K)
    i_post1 <- i_post[,t(matrix(1:ncol(i_post), nrow = K))]
    i_post2 <- vec(matrix(i_post1, nrow = K, byrow = T))
    
    i_stack <- c(i_stack, i_post2)
    j_stack <- c(j_stack, rep(1:((t_max-i)*K), each = K))
    x_stack <- c(x_stack, rep(vec(B_t[(i-1)*K + c(1:K),1:K]), t_max-i))
  }
  H_B <- sparseMatrix(i = i_stack, j = j_stack, x = x_stack, dims = c(nT,nT), repr = "T")
  
  return(H_B)
}

#' @export
Make_cB <- function(y0,B,K,p,t_max, LF = rep(0,K*t_max)){
  # w_t is used only in MT model to account for the W_t^{1/2} * A
  b_intercept = B[,1]
  B_x =  B[,2:ncol(B)]
  #makeRegressor
  y0 <- tail(y0, p)
  xt <- NULL
  if (p > 0){
    y0_add <- rbind(y0, matrix(0, nrow = p, ncol = K ))
    for (i in c(1:p)){
      xt <- cbind(xt, rbind(vec( as.matrix( t(y0_add)[,(p+i-1):i]))) )
    }
  } 
  c_B <- rep(b_intercept, t_max) + LF
  for (t in c(1:p)){
    c_B[(t-1)*K + c(1:K) ] <- c_B[(t-1)*K + c(1:K) ] + (B_x %*% xt[,t])
  }
  return(c_B)
}

#' @export
Make_bandSparse <- function(t_row, aggregation){
  if (aggregation == "average"){
    M_aggregate <-   bandSparse(t_row, k = -c(0:2), 
                                diag = list(rep(1/3, t_row), 
                                            rep(1/3, t_row),
                                            rep(1/3, t_row)),
                                symm=FALSE)
  }
  if (aggregation == "triangular"){
    M_aggregate <-   bandSparse(t_row, k = -c(0:4), 
                                diag = list(rep(1/3, t_row), 
                                            rep(2/3, t_row),
                                            rep(1, t_row),
                                            rep(2/3, t_row),
                                            rep(1/3, t_row)),
                                symm=FALSE)
  }
  return(M_aggregate)
}
