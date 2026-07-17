library(MASS)

save_path <- # set path
load_path <- # set path

f <- function(u, k) {
  (0.5 + 0.25 * cos(k * pi * u)) 
}

g <- function(u, k) {
  (0.3 + 0.15 * sin(k * pi * u))
}

theta_ep <- function(u) {
  0.4 + 0.2 * sin(pi * u)
}

theta_xi <- function(u) {
  0.5 + 0.25 * sin(2 * pi * u)
}

sigma_ep <- function(u) {
  0.2 + (0.5 + 0.25 * sin(2 * pi * u))
}

sigma_xi <- function(u) {
  0.5 + (0.4 + 0.2 * cos(2 * pi * u)) 
}

get_err_shocks <- function(n, rho) {
  mean_vec <- c(0, 0)
  corr_mat <- matrix(c(1, rho, rho, 1), nrow = 2)
  shocks <- mvrnorm(n = n, mu = mean_vec, Sigma = corr_mat)
  return(shocks)
}


get_err_processes <- function(n, rho) {
  err_shocks <- get_err_shocks(n, rho)
  eta_ep <- err_shocks[, 1]
  eta_xi <- err_shocks[, 2]
  ep <- rep(NA, n)
  xi <- rep(NA, n)
  ep_init <- 0.0001
  xi_init <- 0.0001
  ep[1] <- theta_ep(1 / n) * ep_init + eta_ep[1]
  xi[1] <- theta_xi(1 / n) * xi_init + eta_xi[1]
  for (t in 2:n){
    ep[t] <- theta_ep(t / n) * ep[t - 1] + eta_ep[t]
    xi[t] <- theta_xi(t / n) * xi[t - 1] + eta_xi[t]
  }
  return(list(ep = ep, xi = xi))
}



save_xyz <- function(n, rho, k, n_realization) {
  for (r in 1:n_realization) {
    err_processes <- get_err_processes(n, rho)
    ep <- err_processes$ep
    xi <- err_processes$xi
    x <- rep(NA, n)
    y <- rep(NA, n)
    z <- rep(1, n) # constant
    for (t in 1:n){
      x[t] <- f(t / n, k) + sigma_ep(t / n) * ep[t]
      y[t] <- g(t / n, k) + sigma_xi(t / n) * xi[t]
    }
    df <- data.frame(x = x, y = y, z = z, ep = ep, xi = xi)

    # Save the data frame to an RData file

    # not great naming conventions but ok
    # 0 becomes "0", 0.3 becomes "3", 0.6 becomes "6", 0.9 becomes "9"
    corr_name <- as.character(rho * 10)
    # 1 becomes "1", 1.5 becomes "15", 2 becomes "2"
    k_str <- gsub("\\.", "", as.character(k))

    filename <- paste("dgp", "corrshock", "n", n,
                      "k", k_str, "corr", corr_name,
                      "r", r,
                      sep = "_")
    obj_name <- paste0(save_path, filename, ".rds")
    print(paste0("Saving: ", obj_name))
    saveRDS(df, obj_name)
  }
}




load_xyz <- function(n, rho, r, k) {
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

# parameters
# conditional independence test to use
k_complexities <- c(1, 2, 3, 4) # regression complexity parameters
n_samples <- c(250, 500, 750, 1000) # number of samples
n_realizations <- c(100) # number of realizations for empirical rejection rates
correlations <- c(0, 0.3, 0.6, 0.9) # correlations
set.seed(2029) # seed for reproducibility and so same errors across complexities


param_grid <- expand.grid(n_samples = n_samples,
                          correlations = correlations,
                          k_complexities = k_complexities,
                          n_realizations = n_realizations)

for (i in seq_len(nrow(param_grid))) {
  df <- save_xyz(n = param_grid[i, ]$n_samples,
                 rho = param_grid[i, ]$correlations,
                 k = param_grid[i, ]$k_complexities,
                 n_realization = param_grid[i, ]$n_realizations)

}

# load
# df <- load_xyz(n=500, rho=0.9, r=3)