library(dplyr)
library(rlist)
library(parallel)

save_path <- # set path
load_path <- # set path
dev_path <- # set path
source(paste0(dev_path, "dGCM_test.R"))
source(paste0(dev_path, "sieve.R"))

run_sims <- function(params) {
  print(params)

  ci_test <- params$ci_tests
  n_sample <- params$n_samples
  alpha <- params$alphas
  k_complexity <- params$k_complexities
  correlation <- params$correlations
  n_realization <- params$n_realizations

  # for minimum volatility method for selecting lag window for cov estimation
  # lag window candidates to consider
  mv_seq_by <- 1 # consider contiguous lag window sizes, increments of 1
  mv_num_neighbors <- 12 # for calculation of volatility criterion
  power_lw <- 3 / 4 # power for largest candidate lag window size
  lw_candidates <- seq(1, max(1, floor(n_sample^power_lw)), mv_seq_by)

  # number of simulations for the Gaussian process for estimating quantile
  n_boot <- 4999

  # grid of numbers of basis functions for the sieve estimator
  c_vec <- seq(2, 10, 2)
  d_vec <- seq(2, 10, 2)

  # for sieve estimator grid (z,u)
  n_esti <- 250

  # Get the current datetime
  start_time <- Sys.time()
  if (ci_test == "sievedgcm") {
    res <- run_sieve_dgcm_sims(n_samples = n_sample,
                               n_esti = n_esti, n_realization = n_realization,
                               n_boot = n_boot, alpha = alpha,
                               k_complexity = k_complexity,
                               correlation = correlation,
                               c_vec = c_vec, d_vec = d_vec,
                               cv_buffer = 1, lw_candidates = lw_candidates,
                               fp = dev_path, mv_se_param = mv_num_neighbors,
                               plot_residuals = FALSE,
                               x_name = c("x"), y_name = c("y"),
                               z_name = c("z"))
  } else if (ci_test == "oracledgcm") {
    res <- run_oracle_dgcm_sims(n_samples = n_sample,
                                n_realization = n_realization,
                                n_boot = n_boot, alpha = alpha,
                                k_complexity = k_complexity,
                                correlation = correlation,
                                lw_candidates = lw_candidates, fp = dev_path,
                                mv_se_param = mv_num_neighbors,
                                plot_residuals = FALSE,
                                x_name = c("x"), y_name = c("y"),
                                z_name = c("z"))

  } else if (ci_test == "gcm" || ci_test == "rpt" || ci_test == "kci") {
    res <- run_iid_ci_tests_sims(n_samples = n_sample,
                                 n_realization = n_realization,
                                 alpha = alpha, k_complexity = k_complexity,
                                 correlation = correlation, fp = dev_path,
                                 plot_residuals = FALSE,
                                 x_name = c("x"), y_name = c("y"),
                                 z_name = c("z"), ci_test = ci_test)
  } else {
    print("Only ci_tests: sievedgcm, oracledgcm, gcm, rpt, kci are supported")
    break
  }

  end_time <- Sys.time()

  # Format datetimes as YYYYMMDDHHMMSS to avoid invalid characters in filenames
  start_time_formatted <- format(start_time, "%Y%m%d%H%M%S")
  end_time_formatted <- format(end_time, "%Y%m%d%H%M%S")


  if (alpha == 0.025) {
    alpha_name <- "025"
  }else if (alpha == 0.05) {
    alpha_name <- "05"
  }else {
    print("Only alphas: 0.025, 0.05 are supported")
    break
  }




  if (correlation == 0) {
    correlation_name <- "00"
  }else if (correlation == 0.3) {
    correlation_name <- "03"
  }else if (correlation == 0.6) {
    correlation_name <- "06"
  }else if (correlation == 0.9) {
    correlation_name <- "09"
  }else {
    print("Only parameters: 0, 0.3, 0.6, 0.9 are supported")
    break
  }

  k_str <- gsub("\\.", "", as.character(k_complexity))


  # Create a descriptive filename including parameter values
  file_name <- paste0("ci_test_", ci_test,
                      "_dgp_", "corrshock",
                      "_n_sample_", n_sample,
                      "_alpha_", alpha_name,
                      "_k_", k_str,
                      "_corr_", correlation_name,
                      "_st_", start_time_formatted,
                      "_et_", end_time_formatted)
  print(file_name)

  # Save the list as an RData file
  list.save(res, file = paste0(save_path, file_name, ".RData"))
  print(paste("Saved", file_name))
}




run_sieve_dgcm_sims <- function(n_samples, n_esti, n_realization, n_boot,
                                alpha, k_complexity, correlation,
                                c_vec, d_vec, cv_buffer,
                                lw_candidates, fp, mv_se_param,
                                plot_residuals, x_name, y_name, z_name) {


  rej_counter <- 0
  res_list <- list()
  for (realization in 1:n_realization){
    print(paste("n_realization:", realization, "/", n_realization))
    df <- load_xyz(n = n_samples, k = k_complexity,
                   rho = correlation, r = realization)

    cv_res_x <- cv_skip_sampling(regressors = data.frame(df[, z_name]), # nolint
                                 regressand = df[, x_name],
                                 c_vec = c_vec, d_vec = d_vec, b_time = "Legen",
                                 b_timese = "Legen", mp_type = "algeb",
                                 cv_buffer = cv_buffer,
                                 s = 1, n_esti = n_esti)
    c_sieve_x <- cv_res_x$c_min_mse
    d_sieve_x <- cv_res_x$d_min_mse

    x_fitted <- esti_ts(regressors = data.frame(df[, z_name]), # nolint
                        regressand = df[, x_name],
                        c = c_sieve_x, d = d_sieve_x,
                        b_time = "Legen", b_timese = "Legen",
                        mp_type = "algeb", s = 1, n_esti = n_esti)

    cv_res_y <- cv_skip_sampling(regressors = data.frame(df[, z_name]), # nolint
                                 regressand = df[, y_name],
                                 c_vec = c_vec, d_vec = d_vec,
                                 b_time = "Legen", b_timese = "Legen",
                                 mp_type = "algeb", cv_buffer = cv_buffer,
                                 s = 1, n_esti = n_esti)
    c_sieve_y <- cv_res_y$c_min_mse
    d_sieve_y <- cv_res_y$d_min_mse
    y_fitted <- esti_ts(regressors = data.frame(df[, z_name]), # nolint
                        regressand = df[, y_name],
                        c = c_sieve_y, d = d_sieve_y,
                        b_time = "Legen", b_timese = "Legen",
                        mp_type = "algeb", s = 1, n_esti = n_esti)
    x_resid <- df[, x_name] - x_fitted
    y_resid <- df[, y_name] - y_fitted

    res <- dgcm_test(x_on_z_resid = x_resid, # nolint
                     y_on_z_resid = y_resid,
                     lw_candidates = lw_candidates, mv_se_param, fp = dev_path,
                     alpha = alpha, nsim = n_boot,
                     plot.residuals = plot_residuals)

    res["c_sieve_x"] <- c_sieve_x
    res["d_sieve_x"] <- d_sieve_x
    res["c_sieve_y"] <- c_sieve_y
    res["d_sieve_y"] <- d_sieve_y

    if (res$reject[[1]]) {
      rej_counter <- rej_counter + 1
    }
    res_list[[paste0("realization", realization)]] <- res
  }
  return(list(rej_counter = rej_counter, res_list = res_list))
}

run_iid_ci_tests_sims <- function(n_samples, n_realization, alpha,
                                  k_complexity, correlation, fp,
                                  plot_residuals,
                                  x_name, y_name, z_name, ci_test) {

  if (ci_test == "gcm") {
    library(GeneralisedCovarianceMeasure)
  }
  if (ci_test == "rpt") {
    library(CondIndTests)
  }
  if (ci_test == "kci") {
    library(CondIndTests)
  }

  rej_counter <- 0
  res_list <- list()

  for (realization in 1:n_realization){
    print(paste("n_realization:", realization, "/", n_realization))
    df <- load_xyz(n = n_samples, k = k_complexity,
                   rho = correlation, r = realization)

    if (ci_test == "gcm") {
      res <- gcm.test(df[, x_name], df[, y_name], df[, z_name],
                      regr.method = "gam", alpha = alpha)
      if (res$reject) {
        rej_counter <- rej_counter + 1
      }
      res_list[[paste0("realization", realization)]] <- res
    }
    if (ci_test == "rpt") {
      res <- CondIndTest(df[, x_name], df[, y_name], df[, z_name],
                         method = "ResidualPredictionTest", alpha = alpha)
      if (res$pvalue < alpha) {
        rej_counter <- rej_counter + 1
      }
      res_list[[paste0("realization", realization)]] <- res
    }
    if (ci_test == "kci") {
      res <- CondIndTest(df[, x_name], df[, y_name], df[, z_name],
                         method = "KCI", alpha = alpha)
      if (res$pvalue < alpha) {
        rej_counter <- rej_counter + 1
      }
      res_list[[paste0("realization", realization)]] <- res
    }
  }
  return(list(rej_counter = rej_counter, res_list = res_list))
}




run_oracle_dgcm_sims <- function(n_samples, n_realization, n_boot, alpha,
                                 k_complexity, correlation,
                                 lw_candidates, fp, mv_se_param,
                                 plot_residuals, x_name, y_name, z_name) {
  rej_counter <- 0
  res_list <- list()
  for (realization in 1:n_realization){
    print(paste("n_realization:", realization, "/", n_realization))

    df <- load_xyz(n = n_samples, k = k_complexity,
                   rho = correlation, r = realization)

    # perfectly estimated time-varying regression functions
    f_z_perf_est <- mapply(f, df[, z_name], seq(1, n_samples) / n_samples,
                           MoreArgs = list(k = k_complexity))
    g_z_perf_est <- mapply(g, df[, z_name], seq(1, n_samples) / n_samples,
                           MoreArgs = list(k = k_complexity))                           
    x_fitted <- f_z_perf_est
    y_fitted <- g_z_perf_est
    x_resid <- df[, x_name] - x_fitted
    y_resid <- df[, y_name] - y_fitted

    res <- dgcm_test(x_on_z_resid = x_resid, # nolint
                     y_on_z_resid = y_resid,
                     lw_candidates = lw_candidates,
                     mv_se_param, fp = dev_path,
                     alpha = alpha, nsim = n_boot,
                     plot.residuals = plot_residuals)

    if (res$reject[[1]]) {
      rej_counter <- rej_counter + 1
    }
    res_list[[paste0("realization", realization)]] <- res
  }
  return(list(rej_counter = rej_counter, res_list = res_list))
}


# up to date as of april 23, 2025
f <- function(z, u, k) {
  (0.5 + 0.25 * cos(2 * pi * u)) * exp(-z^2) * sin(k * z)
}

g <- function(z, u, k) {
  (0.3 + 0.15 * sin(pi * u)) * exp(-z^2) * cos(k * z)
}

load_xyz <- function(n, k, rho, r) {
  # Load the saved data frame

  # not great naming conventions but ok
  # 0 becomes "0", 0.3 becomes "3", 0.6 becomes "6", 0.9 becomes "9"
  corr_name <- as.character(rho * 10)
  # 1 becomes "1", 1.5 becomes "15", 2 becomes "2"
  k_str <- gsub("\\.", "", as.character(k))
  
  filename <- paste("dgp", "corrshock", "n", n,
                    "k", k_str, "corr", corr_name,
                    "r", r,
                    sep = "_")
  obj_name <- paste0(load_path, filename, ".rds")
  print(paste0("Loading: ", obj_name))

  df <- readRDS(obj_name)

  return(df)
}

######################## run simulations ########################

# parameters
# conditional independence test to use
ci_tests <- c("gcm", "rpt", "oracledgcm", "sievedgcm") # "kci" is far too slow
alphas <- c(0.05) # significance levels
k_complexities <- c(1, 2, 3, 4) # regression complexity parameters
n_samples <- c(250, 500, 750, 1000) # number of samples
n_realizations <- c(100) # number of realizations for empirical rejection rates

### parallelized simulations
correlations <- c(0.9) # 0, 0.3, 0.6, 0.9 # correlations
set.seed(8)
### seeds for reproducibility
# 5 for correlation 0
# 6 for correlation 0.3
# 7 for correlation 0.6
# 8 for correlation 0.9

param_grid <- expand.grid(alphas = alphas,
                          k_complexities = k_complexities,
                          n_samples = n_samples,
                          correlations = correlations,
                          n_realizations = n_realizations,
                          ci_tests = ci_tests)

# Convert to list for easier parallel processing
param_list <- split(param_grid, seq_len(nrow(param_grid)))

# Detect available cores
num_cores <- detectCores() - 2  # Use all but 2 cores

# Run simulations in parallel
mclapply(param_list, run_sims, mc.cores = num_cores)
