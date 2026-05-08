#' @export
metrics <- function(Y_obs_i_j, Y_sim_i_j, Y_mu_i_j, Y_sd_i_j){
  MSE <- (mean(Y_sim_i_j) - Y_obs_i_j)^2
  MAE <- abs(mean(Y_sim_i_j) - Y_obs_i_j)
  Quantile <- quantile(Y_sim_i_j, probs = c(0.05,0.10,0.5,0.90, 0.95))
  emp_CDF <- ecdf(x = Y_sim_i_j)(Y_obs_i_j)
  # predict_den <- stats::density(Y_sim_i_j, adjust = 1)
  # 
  # id <- which(predict_den$y == 0);
  #     if (length(id)>0) {
  #       for (fix in id){
  #         predict_den$y[fix] <- mean(predict_den$y[(fix-1):(fix+1)])
  #       }
  #     }
  # appximate_density <- smooth.spline(x = predict_den$x, y = log(predict_den$y), df = 10)
  # log_ps1 <- - predict(appximate_density, Y_obs_i_j)$y # Smaller is better
  
  log_ps1 <- - log(mean(dnorm(Y_obs_i_j, mean = Y_mu_i_j, sd = Y_sd_i_j)))
  
  # fitdistrplus::fitdist(Y_sim_i_j, distr = "sstd", 
  #                       method="mle", start=list(mu = 0, sigma = 1, skew = 1, shape = 5))
  predict_den_skewt <- rugarch::fitdist(distribution = "sstd", Y_sim_i_j, control=list())
  log_ps2 <- - log(rugarch::ddist(distribution = "sstd", y = Y_obs_i_j,
                mu = predict_den_skewt$pars[1],
                sigma = predict_den_skewt$pars[2],
                skew = predict_den_skewt$pars[3],
                shape = predict_den_skewt$pars[4]))
  # log_ps2 <- - dnorm(x = Y_obs_i_j, log = TRUE, 
  #                    mean = mean(Y_sim_i_j), sd = sd(Y_sim_i_j)) # Smaller is better
  
  samples1 <- Y_sim_i_j
  samples2 <- sample(samples1, size = length(samples1), replace = T)
  CRPS <- mean(abs(samples1 - Y_obs_i_j)) - 0.5*mean(abs(samples1 - samples2)) # Smaller is better
  
  alpha <- seq(0.05, 0.95, by = 0.05)
  Q_tau <- quantile(Y_sim_i_j, probs = alpha)
  QS <- (Q_tau - Y_obs_i_j) * ( (Y_obs_i_j < Q_tau) - alpha ) # Smaller is better
  qwCRPS_2t <- 2/length(alpha) * sum(  (2 * alpha - 1)^2 * QS )
  qwCRPS_rt <- 2/length(alpha) * sum(  alpha^2 * QS )
  qwCRPS_lt <- 2/length(alpha) * sum(  (1 - alpha)^2 * QS )
  
  return(list(MSE = MSE, MAE = MAE,
              Quantile = Quantile, emp_CDF= emp_CDF,
              log_ps1 = log_ps1, log_ps2 = log_ps2,
              CRPS = CRPS, qwCRPS_2t = qwCRPS_2t,
              qwCRPS_rt = qwCRPS_rt, qwCRPS_lt = qwCRPS_lt))
}

#' #' @export
#' recursive_seperate <- function(y_current, y_next, t_pred, 
#'                                K, p, dist = "OT", SV = T, 
#'                                aggregation = "triangular", idq, outname = NULL, 
#'                                samples = 60000, burnin = 10000, thin = 5){
#'   t_max = nrow(y)
#'   if (is.null(outname)) outname = paste("Recursive_", dist, "_", SV, "_", t_max, ".RData", sep = "")
#'   time_current <- t_max
#'   
#'   y = as.matrix(y_current[(p+1):t_max,])
#'   y0 = as.matrix(y_current[1:p,, drop = FALSE ])
#'   for (i in c(1:ncol(y))){
#'     y0[is.na(y0[,i]) ,i] <- mean(y[,i], na.rm = TRUE)
#'   }
#'   prior <- get_prior(y, p = p, dist=dist, SV = SV, aggregation = aggregation, idq = idq)
#'   inits <- get_init(prior, samples = samples, burnin = burnin, thin = thin)
#'   
#'   if (SV) {
#'     Chain <- BMFVAR.SV(y = y, K = K, p = p, dist = dist,
#'                      y0 = y0, prior = prior, inits = inits)
#'   } else {
#'     Chain <- BMFVAR.novol(y = y, K = K, p = p, dist = dist, 
#'                           y0 = y0, prior = prior, inits = inits)
#'   }
#' 
#'   forecast_err <- forecast_density(Chain = Chain, y_next = y_next, t_current = t_max, t_pred = t_pred)
#'   out_recursive <- list(time_id = time_current,
#'                         forecast_err = forecast_err,
#'                         dist = dist,
#'                         SV = SV,
#'                         y_current = y_current,
#'                         y_next = y_next)
#'   save(out_recursive, file = outname)
#'   return( out_recursive)
#' 
#' }
#' 
