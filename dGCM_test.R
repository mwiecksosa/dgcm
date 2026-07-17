library(ggplot2)

dgcm_test <- function(x_on_z_resid, y_on_z_resid, lw_candidates, mv_se_param,
                      fp, alpha = 0.05, nsim = 4999L,
                      plot.residuals = TRUE) {

  x_on_z_resid <- data.frame(x_on_z_resid)
  y_on_z_resid <- data.frame(y_on_z_resid)
  n_x <- dim(x_on_z_resid)[1]
  n_y <- dim(y_on_z_resid)[1]
  d_x <- dim(x_on_z_resid)[2]
  d_y <- dim(y_on_z_resid)[2]
  if (n_x != n_y) {
    stop("Residuals must be the same length")
  }

  # number of residuals
  n <- n_x

  # check lag windows
  if (length(lw_candidates) == 0) {
    stop("Lag-window candidates lw_candidates must contain at least one value")
  }

  if (any(is.na(lw_candidates))) {
    stop("Lag-window candidates lw_candidates cannot contain NA values")
  }

  if (any(lw_candidates < 1) || max(lw_candidates) > n) {
    stop("Lag-window candidates lw_candidatesmust be between 1 and n")
  }

  lw_candidates <- sort(unique(lw_candidates)) # enforce increasing order

  # check number of simulations is positive
  if (nsim < 1) {
    stop("Number of simulations nsim must be positive")
  }

  # check significance level is between 0 and 1
  if (alpha <= 0 || alpha >= 1) {
    stop("Significance level alpha must be between 0 and 1")
  }

  # check plot.residuals is a logical value
  if (!is.logical(plot.residuals)) {
    stop("Plot residuals decision variable plot.residuals must be a logical value")
  }

  # check mv_se_param is positive 
  if (mv_se_param < 1) {
    stop("Minimum volatility parameter mv_se_param must be positive")
  }

  # check fp is a string
  if (!is.character(fp)) {
    stop("File path fp must be a string")
  }


  resid_prod <- x_on_z_resid * y_on_z_resid
  lw <- mv_method(resid_prod, n, lw_candidates, mv_se_param)
  cumul_var_hat <- estimate_cumul_var(r = resid_prod, lw = lw, n = n)
  s_hat_vec <- rep(NA, nsim)
  TnL <- n - lw + 1
  for (sim in 1:nsim) {    
    # gaussian partial sum process
    g_psum_proc <- rep(NA, n)
    g_psum_t <- 0
    # no estimate for 1:(lw-1)
    for (t in lw:n) {
      if (t == lw) {
        sigma2_hat_t <- cumul_var_hat[t]
      } else {
        sigma2_hat_t <- cumul_var_hat[t] - cumul_var_hat[t - 1]
      }
      g_psum_t <- g_psum_t + rnorm(mean = 0, sd = sqrt(sigma2_hat_t), n = 1)
      g_psum_proc[t] <- g_psum_t
    }
    
    s_hat <- max(abs(g_psum_proc[lw:n])) / sqrt(TnL)
    s_hat_vec[sim] <- s_hat
  }

  boot_quantile <- quantile(s_hat_vec, probs = (1 - alpha))
  max_abs_partial_sum <- max(abs(cumsum(resid_prod[lw:n,1])))
  test_stat <- max_abs_partial_sum / sqrt(TnL)

  p_value <- (sum(s_hat_vec >= test_stat) + 1) / (nsim + 1)
  reject <- test_stat > boot_quantile
  if (plot.residuals) {
    png(paste0(fp, "dgcm_absNormCumsumResid.png"), width = 800, height = 600)
    plot(lw:n, abs((cumsum(resid_prod[lw:n,1]))) / sqrt(TnL),
         xlab = "Index",
         ylab = "Normalized Absolute Value of Cumulative Sum of Residual Products",
         main = "Normalized Absolute Value of Cumulative Sum of Residual Products")
    abline(h = boot_quantile, col = "red")
    dev.off()
  }
  return(list(p_value = p_value,
              test_stat = test_stat,
              reject = reject,
              max_abs_partial_sum = max_abs_partial_sum,
              n = n,
              boot_quantile = boot_quantile,
              lw = lw,
              lw_candidates = lw_candidates,
              mv_se_param = mv_se_param,
              nsim = nsim))
}

mv_method <- function(resid_prod, n, lw_candidates, mv_se_param) {
  # minimum volatility method for choosing lag window size
  num_lws_h <- length(lw_candidates) # number of lag window candidates H
  sigma2_hat_m <- compute_sigma2_hat(resid_prod, lw_candidates, n, num_lws_h)
  vol_vec <- compute_vol_vec(sigma2_hat_m, lw_candidates,
                             num_lws_h, n, mv_se_param)
  lw_star <- lw_candidates[which.min(vol_vec)]
  return(lw_star)
}

compute_sigma2_hat <- function(resid_prod, lw_candidates, n, num_lws_h) {
  # utility function for estimating covariance at each time
  # using each of the different candidate lag windows
  sigma2_hat_m <- matrix(NA, nrow = n, ncol = num_lws_h)
  for (h in 1:num_lws_h) {
    lw_h <- lw_candidates[h]
    cumul_var_hat <- estimate_cumul_var(r = resid_prod, lw = lw_h, n = n)
    for (t in lw_h:n) {
      if (t == lw_h) {
        sigma2_hat_m[t, h] <- cumul_var_hat[t]
      } else {
        sigma2_hat_m[t, h] <- cumul_var_hat[t] - cumul_var_hat[t - 1]
      }
    }
  }
  return(sigma2_hat_m)
}

compute_vol_vec <- function(sigma2_hat_m, lw_candidates,
                            num_lws_h, n, mv_se_param) {
  # utility function for calculating the vol for each candidate lag window
  vol_vec <- rep(NA, num_lws_h)
  lw_max <- max(lw_candidates)
  for (j in 1:num_lws_h) {
    # number of lag windows for covariances
    h2 <- min(num_lws_h, j + mv_se_param)
    h1 <- max(1, j - mv_se_param)
    nc <- 1 / (h2 - h1 + 1)
    se_vec <- rep(NA, n)
    for (t in lw_max:n) {
      # avg covariance
      avg_sigma2_hat_t <- nc * sum(sigma2_hat_m[t, h1:h2])
      se_vec[t] <- sqrt(nc * sum((sigma2_hat_m[t, h1:h2] - avg_sigma2_hat_t)^2))
    }
    vol_vec[j] <- max(se_vec, na.rm = TRUE)
  }
  return(vol_vec)
}

estimate_cumul_var <- function(r, lw, n) {
  cumul_var_hat <- rep(NA, n)
  psum_sqr_sum <- 0
  # no estimate for 1:(lw-1)
  for(k in lw:n) {
    start_idx <- k - lw + 1
    end_idx <- k
    sqr_sum <- (1 / lw) * (sum(r[start_idx:end_idx, 1]))^2
    psum_sqr_sum <- psum_sqr_sum + sqr_sum
    cumul_var_hat[k] <- psum_sqr_sum
  }
  return(cumul_var_hat)
}
