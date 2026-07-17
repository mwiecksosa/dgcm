save_path <- # set path
load_path <- # set path

f <- function(z, u, k) {
  (0.4 + 0.2 * sin(2 * pi * u)) * exp(-z^2) * sin(k * z)
}

theta_e <- function(u) {
  0.45 + 0.3 * sin(2 * pi * u)
}

theta_z <- function(u) {
  0.5 + 0.25 * cos(pi * u)
}

get_err_processes <- function(n, z) {
  eta_ep <- rnorm(n)
  eta_xi <- rnorm(n)
  ep <- rep(NA, n)
  xi <- rep(NA, n)
  ep_init <- 0.0001
  xi_init <- 0.0001
  ep[1] <- theta_e(1 / n) * ep_init + eta_ep[1]
  xi[1] <- theta_e(1 / n) * xi_init + eta_xi[1]
  for (t in 2:n){
    ep[t] <- theta_e(t / n) * ep[t - 1] + eta_ep[t]
    xi[t] <- theta_e(t / n) * xi[t - 1] + eta_xi[t]
  }
  return(list(ep = ep, xi = xi))
}

get_z <- function(n) {
  eta_z <- rnorm(n)
  z <- rep(NA, n)
  z_init <- 0.0001
  z[1] <- theta_z(1 / n) * z_init + eta_z[1]
  for (t in 2:n){
    z[t] <- theta_z(t / n) * z[t - 1] + eta_z[t]
  }
  return(z)
}

sigma_e <- function() {
  0.3
}

save_xyz <- function(n, k, beta, n_realization) {
  for (r in 1:n_realization) {
    z <- get_z(n)
    err_processes <- get_err_processes(n, z)
    ep <- err_processes$ep
    xi <- err_processes$xi
    x <- rep(NA, n)
    y <- rep(NA, n)
    for (t in 1:n){
      x[t] <- f(z[t], t / n, k) + sigma_e() * ep[t]
      y[t] <- f(z[t], t / n, k) + beta * x[t] + sigma_e() * xi[t]
    }
    df <- data.frame(x = x, y = y, z = z, ep = ep, xi = xi)


    # Save the data frame to an RData file

    # not great naming conventions but ok
    # 0 becomes "0", 0.3 becomes "3", 0.6 becomes "6", 0.9 becomes "9"
    effect_size_name <- as.character(beta * 10) 
    # 1 becomes "1", 1.5 becomes "15", 2 becomes "2"
    k_str <- gsub("\\.", "", as.character(k))

    filename <- paste("dgp", "addeffect", "n", n,
                      "k", k_str, "eff", effect_size_name,
                      "r", r,
                      sep = "_")
    obj_name <- paste0(save_path, filename, ".rds")
    print(paste0("Saving: ", obj_name))
    saveRDS(df, obj_name)
  }
}


load_xyz <- function(n, k, beta, r) {
  # Load the saved data frame

  # not great naming conventions but ok
  # 0 becomes "0", 0.3 becomes "3", 0.6 becomes "6", 0.9 becomes "9"
  effect_size_name <- as.character(beta * 10) 
  # 1 becomes "1", 1.5 becomes "15", 2 becomes "2"
  k_str <- gsub("\\.", "", as.character(k))

  filename <- paste("dgp", "addeffect", "n", n,
                    "k", k_str, "eff", effect_size_name,
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
effect_sizes <- c(0, 0.3, 0.6, 0.9) # additive effect sizes
set.seed(2025) # seed for reproducibility and so same errors across complexities


param_grid <- expand.grid(k_complexities = k_complexities,
                          n_samples = n_samples,
                          effect_sizes = effect_sizes,
                          n_realizations = n_realizations)

for (i in seq_len(nrow(param_grid))) {
  df <- save_xyz(n = param_grid[i, ]$n_samples,
                 k = param_grid[i, ]$k_complexities,
                 beta = param_grid[i, ]$effect_sizes,
                 n_realization = param_grid[i, ]$n_realizations)

}

# load
# df <- load_xyz(n=750, k=1, beta=0.9, r=3)