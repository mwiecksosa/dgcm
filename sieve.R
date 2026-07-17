#' "Conditional Independence Testing with a Single Realization of a Multivariate Nonstationary Nonlinear Time Series"
#' Submitted to the Journal of Business and Economic Statistics, 2025
#' 
#' Original code for the sieve model is from the SIMle R package from:
#' "Simultaneous Sieve Inference for Time-Inhomogeneous Nonlinear Time Series Regression."
#' https://arxiv.org/pdf/2112.08545.pdf
#' 
#' functions esti_ts(), general_esti(), esti_beta() are adapted from the SIMle R package.
#' Purpose: original code for sieve estimator is for the autoregressive case only.
#'          We changed the original code to allow for general covariates.
#' functions skip_sample_time_series(), cv_skip_sampling(), predict_with_last_time() are added.
#' Purpose: We added a utility function for forecasting using the 
#'          Also, we introduced a new cross-validation approach for selecting the
#'          parameters of the sieve estimator in this setting.
#'          Lastly, we added a utility function for forecasting using the 
#'          regression function at the last time.
#' User-specified estimation of nonlinear time series regression
#' @description This function estimates nonlinear time series regression
#'
#' @param x (dataframe) regressors with relevant lags already included
#' @param y (vector) regressand
#' @param c number of basis for time input
#' @param d number of basis for variate input 
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param type select type of estimation."nfix" refers to no fix estimation. "fixt" indicates fix time t estimation. 
#'             "fixx" represents fix variate estimation
#' @param fix_num fix_num indicates the use of fixed-value nonlinear time series regression. The default value is 0, which is employed for non-fixed estimation. 
#'        If "fixt" is chosen, it represents a fixed time value. Otherwise, if not selected, it pertains to a fixed variate value
#' @param s s is a positive scaling factor, the default is 1
#' @param n_esti number of points for estimation, the default is 2000
#'
#' @return If "nfix" is selected, the function returns a list where each element is a matrix representing the estimation function in two dimensions. Otherwise, 
#'         if "nfix" is not selected, the function returns a list where each element is a vector representing the estimation function.
#' @export
#'
#' @examples
#' res = fix.fit(regressors=x, regressand=y,c=2, d=2, b_time="tri", b_timese="tri", 
#'               mp_type="algeb", type="nfix", s=1, n_esti=2000)
#' 
#' 

library(stringr)
library(RCurl)
library(splines)

fix.fit <- function(regressors, regressand, c, d, b_time, b_timese, mp_type, type="nfix", fix_num = 0, s = 1, n_esti = 2000){
  basis_candi = c("Legen","Cheby","tri", "cos", "sin", "Cspli", "db1", "db2", "db3", "db4", "db5",
                  "db6", "db7", "db8", "db9", "db10",
                  "db11", "db12", "db13", "db14", "db15",
                  "db16", "db17", "db18", "db19", "db20",
                  "cf1", "cf2", "cf3", "cf4", "cf5"
  )
  
  
  
  if((b_time %in% basis_candi) && (b_timese %in% basis_candi)){
    if(type == "nfix"){
      res = general_esti(regressors=regressors, regressand=regressand, c, d, b_time, b_timese, mp_type, s, n_esti = n_esti)
      return(res)
    }
  }
}



# Estimation: overall (3D plot), fix_t(2D plot), fix_x(2D plot) (auto version, automatically choosing c and d)
general_esti <- function(regressors, regressand, c, d, b_time, b_timese, mp_type, s = 1, n_esti = 2000){ 
  colMax <- function(data){sapply(data, max, na.rm = TRUE)}
  colMin <- function(data){sapply(data, min, na.rm = TRUE)}

  #print("general_esti")
  #print("regressors")
  #print(regressors)
  #print("regressand")
  #print(regressand)

  x <- regressors
  y <- regressand
  
  n <- dim(x)[1]
  r <- dim(x)[2]
  
  uppers <- colMax(x)
  lowers <- colMin(x)
  basis_ti <- list()
  x_i <- list()
  for(k in 1:r){
    # typically, scalar function is s. 
    pos_map = c( "algebp", "logarip") # positive mapping functions
    if(mp_type %in% pos_map){ # positive only
      x_i[[k]] = seq(0, uppers[[k]], length.out = n_esti)
    }else{ # default to this (mapping function with real line as domain)
      x_i[[k]] = seq(lowers[[k]], uppers[[k]], length.out = n_esti) 
    }
  }
  
  basis_ti = matrix(nrow = n_esti) # n_estix1 matrix of NAs
  
  esti_df = esti_beta(x, y, c, d, b_time, b_timese, mp_type, s)
  beta_hat = esti_df[[1]]
  W =  esti_df[[2]]

  df_basis = 0
  if(b_timese == "Cspli"){
    df_basis = Cspline_table(d) # true is d - 2 + or 
    d = dim(df_basis)[2]
  }
  
  # change for wavelet
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  
  basis_xi <- list()
  for(k in 1:r){
    basis_xi[[k]] = matrix(ncol = d_tri) # 1 x d_tri vector of NAs
    
    for(i in 1:n_esti){
      basis_xi[[k]] = rbind(basis_xi[[k]], sqrt(mp_selection(mp_type, x_i[[k]][i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[[k]][i], s)[1], b_timese, df_basis = df_basis))
    }
    basis_xi[[k]] = basis_xi[[k]][-1, ] # all but first row of NAs
    basis_xi[[k]] = as.matrix(basis_xi[[k]], ncol = d_tri) # convert back to matrix 
  }
  
  
  
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      aux_bti = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      aux_bti = as.matrix(df_basis_1[, l1])
    }else{
      aux_bti = bs.gene(b_time, l1, n_esti)
    }
    basis_ti = cbind(basis_ti, aux_bti)
  }
  
  basis_ti = as.matrix(basis_ti[,-1], ncol = c) # all but the first column of NAs
  colnames(basis_ti) = NULL
  
  # 2 dimensinal integrate
  
 
  # n_esti x n_esti matrix of zeros 
  m_hat_ij = matrix(rep(0, (n_esti)^2), nrow = n_esti, ncol = n_esti) 
  # c*n_esti x c*n_esti matrix of zeros 
  B_j = matrix(rep(0, (c*d_tri)^2), nrow = c*d_tri, ncol = c*d_tri)
  
  res_m_hat_ij = list()

  ### fit each 2 dimensional function of time and r-th covariate 
  r <- dim(x)[2] # dimension
  for(k in 1:r){ # for each dimension
    for(i in 1:n_esti){ # for each time
      for(j in 1:n_esti){ # for each value
        ti = matrix(as.numeric(basis_ti[i,]), ncol = 1)
        xj = matrix(as.numeric(basis_xi[[k]][j,]), ncol = 1)
        B_j = B_j + kronecker(ti, xj)%*%t(kronecker(ti, xj))
        aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
        aux_m_hat_ij = t(aux_beta_hat)%*%kronecker(ti, xj) # (1 x rcd) x (rcd x 1)
        m_hat_ij[i, j] = aux_m_hat_ij # scalar
      }
    }
    res_m_hat_ij[[k]] = m_hat_ij # dim k, time i, value j
  }
  
  
  return(res_m_hat_ij)
  
}


esti_beta <- function(regressors, regressand, c, d, b_time, b_timese, mp_type, s = 1){  # b_time indicates type of basis for time t 
  #print("esti_beta")
  x <- regressors
  y <- regressand
  #print("regressand")
  #print(regressand)
           
  # b_timese indicates name of basis for time series
  r = dim(x)[2]
  n = length(y)
  #print("y")
  #print(y)
  Y = matrix(y, ncol = 1)
  W = matrix(nrow = n, ncol = 1) # nx1 vector of NAs
  

  
  # change for Cspline
  df_basis = 0
  if(b_timese == "Cspli"){
    df_basis = Cspline_table(d) # true is d - 2 + or 
    d = dim(df_basis)[2]
  }
  
  # change for wavelet
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  #}
  
  phi_x = matrix(nrow = n, ncol = 1) # nx1 vector of NAs
  aux_phi_x = matrix(ncol = d_tri) # 1xn vector of NAs
  
  
  
  for(j in 1:r){ # j for dimension of x 
    for(i in 1:n){
      aux_phi_x = rbind(aux_phi_x, sqrt(mp_selection(mp_type, x[i,j], s)[2])*select_basis_timese(d, mp_selection(mp_type, x[i,j], s)[1], b_timese, df_basis))
    }
    aux_phi_x = aux_phi_x[-1, ]
    aux_phi_x = as.matrix(aux_phi_x, ncol = d_tri)
    #print(paste("nrow(aux_phi_x)",nrow(aux_phi_x)))
    #print(paste("nrow(phi_x)",nrow(phi_x)))
    phi_x = cbind(phi_x, aux_phi_x)
    aux_phi_x = matrix(ncol = d_tri)
  }
  
  # all but first column of NAs
  phi_x = matrix(phi_x[,-1], ncol = r*d_tri)
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(j in 1:r){ # j for dimension
    for(l1 in 1:c){
      if(b_time == "Cspli"){
        aux_basis_t = as.matrix(df_basis[1:n, l1])
      }else if(b_time %in% wavelet_basis){
        aux_basis_t = as.matrix(df_basis_1[1:n, l1])
      }else{
        aux_basis_t = as.matrix(bs.gene(b_time, l1, n)[1:n,])
      }
      
      
      for(l2 in 1:d_tri){
        W = cbind(W, aux_basis_t*matrix(phi_x[, (j-1)*d_tri + l2], ncol = 1))
      }
    }
  }
  
  
  
  W = W[,-1] # all but first col of NAs
  W = as.matrix(W, ncol = r*c*d_tri)
  colnames(W) = NULL
  #print("W")
  #print(W)
  #print("Y")
  #print(Y)
  beta_hat = solve(t(W)%*%W, tol = 1e-40)%*%t(W)%*%Y # (W^T W)^{-1} W^T Y 
  return(list(fit_beta = beta_hat, d_m = W)) # estimation need to change
}



# Function to skip-sample a time series by skipping samples
skip_sample_time_series <- function(x_ts,y_ts,n_ts,step_size) {
  #print("skip_sample_time_series")
  #print(paste("step_size",step_size))
  #print(paste("n_ts",n_ts))

  #print("x_ts")
  #print(x_ts)
  #print("y_ts")
  #print(y_ts)
  
 
  # List to store skip-sampled series
  x_skip_samples <- list()
  y_skip_samples <- list()
  # Loop over possible starting points (from 1 to step_size)
  for (start in 1:step_size) {
    # skip-sample the series by skipping samples
    x_skip_sample <- x_ts[seq(start, n_ts, by = step_size),]
    y_skip_sample <- y_ts[seq(start, n_ts, by = step_size)]
    
    # Store skip-sampled series into a list
    x_skip_samples[[paste0("skip_sample_", start)]] <- x_skip_sample
    y_skip_samples[[paste0("skip_sample_", start)]] <- y_skip_sample
    
  }
  
  return(list(x_skip_samples=x_skip_samples,y_skip_samples=y_skip_samples))
}


cv_skip_sampling <- function(regressors, regressand, 
                             c_vec, d_vec, b_time, 
                             b_timese, mp_type, cv_buffer,
                             s=1, n_esti = 2000){  

  #print("cv_skip_sampling")
  
  # cross-validation via skip-sampling
  colMax <- function(data){sapply(data, max, na.rm = TRUE)}
  colMin <- function(data){sapply(data, min, na.rm = TRUE)}
  

    
  
  x <- regressors
  y <- regressand

  #print("x")
  #print(x)
  #print("y")
  #print(y)
  
  uppers <- colMax(x)
  lowers <- colMin(x)
  
  n = dim(x)[1] # sample size
  r = dim(x)[2] # dimension

  #print(paste("r",r))
  #print(paste("n",n))

  x_i <- list()
  for(k in 1:r){
    pos_map = c( "algebp", "logarip") 
    if(mp_type %in% pos_map){
      x_i[[k]] = seq(0, uppers[[k]], length.out = n_esti)
    }else{
      x_i[[k]] = seq(lowers[[k]], uppers[[k]], length.out = n_esti) 
    }
  }
  
  
  # if you want a gap of cv_buffer between the time series 
  # then this is the number you should use for subsampling the original time series
  skip_by <- 2*(cv_buffer+1) 
  #print(paste("skip_by",skip_by))
  #print(paste("r",r))
  
  skip_sample_res <- skip_sample_time_series(x_ts=x,y_ts=y,n_ts=n, 
                                             step_size=skip_by) 
  x_skip_samples <- skip_sample_res$x_skip_samples
  y_skip_samples <- skip_sample_res$y_skip_samples

  #print("x_skip_samples")
  #print(x_skip_samples)
  #print("y_skip_samples")
  #print(y_skip_samples)
  #print("after skip_sample_time_series")
  #for(c in c_vec){
  #  for(d in d_vec){
  config_mse_vec<-c()
  param_grid <- expand.grid(c_sieves=c_vec,d_sieves=d_vec)
  for(config_num in 1:nrow(param_grid)){
    
    #print("loop1")
    c <- param_grid$c_sieves[config_num]
    d <- param_grid$d_sieves[config_num]  
    print(paste("c_d_config (",config_num,"/",nrow(param_grid), "): c =", c, "d =", d))

    folds_mse_vec<-c()
    for(main_start in 1:skip_by){
      #print("loop2")
      
      main_x_ss <- data.frame(x_skip_samples[[paste0("skip_sample_", main_start)]])
      main_y_ss <- y_skip_samples[[paste0("skip_sample_", main_start)]]

      #print("main_x_ss")
      #print(main_x_ss)
      #print("main_y_ss")
      #print(main_y_ss)

      # estimated function 
      esti_func = general_esti(regressors=main_x_ss, regressand=main_y_ss, 
                               c, d, b_time, 
                               b_timese, mp_type, s, n_esti=n_esti) 
      
      
      #twin_start <- main_start+cv_buffer+1 # start time of the twin ts
      #twin_x_ss <- x_skip_samples[[paste0("skip_sample_", twin_start)]]
      #twin_y_ss <- y_skip_samples[[paste0("skip_sample_", twin_start)]]
      
      fited_res = c()
      sqr_err = c()
      
      for(i in seq(main_start, n, by = skip_by)){
        #print("loop3")
        
        rsum = 0
        #print(paste("r",r))
        for(j in 1:r){ 
          #print("loop4")
          
          # j for dimension
          # i for time
          # last for value of this covariate (dimension)
          # evaluate estimated function f at closest point x_i to x i.e. f(x_i)
          

          # "twin series" from the middle times             
          closest_time <- unique(which(abs((1:n_esti)/n_esti - (i+cv_buffer+1)/n )
                                       == min(abs((1:n_esti)/n_esti - (i+cv_buffer+1)/n ) )))[1]
          # covariate value
          closest_val <- unique(which(abs(x_i[[j]]-x[i+cv_buffer+1,j]) 
                                      == min(abs(x_i[[j]]-x[i+cv_buffer+1,j]))))[1]
          rsum = rsum + esti_func[[j]][closest_time, closest_val]
        }
        
        fited_res[i] = rsum
        sqr_err[i] = (y[i+cv_buffer+1] - fited_res[i])^2
      }
      folds_mse_vec[main_start]<-mean(sqr_err,na.rm=TRUE)
      
    }
    config_mse_vec[config_num]<-mean(folds_mse_vec)
  }
  param_grid['mse']<-config_mse_vec
  c_min_mse <- param_grid[which.min(param_grid[['mse']]),'c_sieves']
  d_min_mse <- param_grid[which.min(param_grid[['mse']]),'d_sieves']
  
  
  


  return(list(c_min_mse=c_min_mse,d_min_mse=d_min_mse))
}



# make prediction on a new sample at rescaled time 1 e.g. for nowcasting or forecasting
predict_with_last_time <- function(new_data, esti_func, mp_type, n_esti = 2000){

  new_data <- as.matrix(new_data) # data to predict

  # esti_func is the estimated function from the training data from general_esti
    # esti_func = general_esti(regressors=main_x_ss, regressand=main_y_ss, 
    #                          c, d, b_time, 
    #                          b_timese, mp_type, s, n_esti=n_esti) 

  # this code is last part of esti_ts
  # but changed so used on new_data not the data used to fit the model
  # and only uses the time-varying regression function at rescaled time 1 (i.e. last time) 

  x <- new_data
  # get maximum and minimum of each column of x
  # Get max and min of each column
  uppers <- apply(x, 2, max)
  lowers <- apply(x, 2, min)


  n = dim(x)[1] 
  r = dim(x)[2]

  x_i <- list()
  for(k in 1:r){
    pos_map = c( "algebp", "logarip") 
    if(mp_type %in% pos_map){
      x_i[[k]] = seq(0, uppers[[k]], length.out = n_esti)
    }else{
      x_i[[k]] = seq(lowers[[k]], uppers[[k]], length.out = n_esti) 
    }
  }
  
  preds = c() 
  for(i in 1:n){
    rsum = 0
    for(j in 1:r){ 
      # j for dimension
      # i for time
      # last for value of this covariate (dimension)
        # originally it was:
        # unique(which(abs((1:n_esti)/n_esti - i/n) == min(abs((1:n_esti)/n_esti - i/n) )))[1]
        # however, - 1 instead of - i/n for closest_time to get index of rescaled time 1 
        # this is equal to n_esti for all times i, so just set closest_time to n_esti      
        # next, evaluate estimated function f at closest point x_i to x i.e. f(x_i)
        # [1] at end to choose the first of the minima if there are ties
      closest_time <- n_esti
      closest_val <- unique(which(abs(x_i[[j]]-x[i,j]) == min(abs(x_i[[j]]-x[i,j]))))[1]
      #print(paste("x[i,j]: ", x[i,j], " x_i[[j]]: ", x_i[[j]][closest_val]))
      #print(paste("j: ", j, " i: ", i, " closest_time: ", closest_time, " closest_val: ", closest_val))
      rsum = rsum + esti_func[[j]][closest_time, closest_val]
    }
    preds[i] = rsum
    
  }
  return(preds)
}

# Get estimated time series (for prediction)
esti_ts <- function(regressors, regressand, c, d, b_time, b_timese, mp_type, s=1, n_esti = 2000){  
  colMax <- function(data){sapply(data, max, na.rm = TRUE)}
  colMin <- function(data){sapply(data, min, na.rm = TRUE)}
  
  x <- regressors
  y <- regressand
           
  uppers <- colMax(x)
  lowers <- colMin(x)

  n = dim(x)[1] # sample size
  r = dim(x)[2]
  fited_res = c()
  x_i <- list()
  for(k in 1:r){
    pos_map = c( "algebp", "logarip") 
    if(mp_type %in% pos_map){
      x_i[[k]] = seq(0, uppers[[k]], length.out = n_esti)
    }else{
      x_i[[k]] = seq(lowers[[k]], uppers[[k]], length.out = n_esti) 
    }
  }
  
  esti_func = general_esti(regressors=x, regressand=y, c, d, b_time, b_timese, mp_type, s, n_esti) # estimated function 
  for(i in 1:n){
    rsum = 0
    for(j in 1:r){ 
      # j for dimension
      # i for time
      # last for value of this covariate (dimension)
      # evaluate estimated function f at closest point x_i to x i.e. f(x_i)
      # [1] at end to choose the first of the minima if there are ties

      
      # do we need to map it to the closest time too???
      closest_time <- unique(which(abs((1:n_esti)/n_esti - i/n) == min(abs((1:n_esti)/n_esti - i/n) )))[1]
      closest_val <- unique(which(abs(x_i[[j]]-x[i,j]) == min(abs(x_i[[j]]-x[i,j]))))[1]
      rsum = rsum + esti_func[[j]][closest_time, closest_val]
    }
    fited_res[i] = rsum
    
  }
  return(fited_res)
}





# Mapping selection 

mp_selection <- function(mp_type, x, s = 1){ # 
  if(mp_type == "algeb"){
    return(c(algeb_u(x, s), algeb_up(x, s))) # first is mapping function, second is derivative of mapping function 
    
  }else if(mp_type == "logari"){
    return(c(logari_u(x, s), logari_up(x, s)))
    
    #else if(mp_type == "algebp"){
    #return(c(algeb_u_posi(x, s), algeb_up_posi(x, s)))
    #}else if(mp_type == "logarip"){
    #  return(c(logari_u_posi(x, s), logari_up_posi(x, s)))
    
  }else{
    return(stop("Invalid option!"))
  }
  
}



# algebraic mapping R to [0,1]
algeb_u <- function(x, s){
  return(x/(2*sqrt(x^2+s^2)) + 1/2) 
}
algeb_up <- function(x, s){
  return((s^2)/(2*((x^2+s^2)^(3/2))))
}



# algebraic mapping R+ to [0,1]
algeb_u_posi <- function(x, s){
  return(x/(x+s)) 
}
algeb_up_posi <- function(x, s){
  return(s/((x+s)^2))
}


# The logarithmic mapping from R to [0, 1]:
logari_u <- function(x, s){
  return((tanh(x/s) + 1)/2) 
}
logari_up <- function(x, s){
  return((1-(tanh(x/s)^2))/(2*s))
}



# The logarithmic mapping from R+ to [0, 1]:
logari_u_posi <- function(x, s){
  return(tanh(x/s)) 
}
logari_up_posi <- function(x, s){
  return((1-(tanh(x/s)^2))/(s))
}



# Adding Cspline and wavelets
select_basis_timese <- function(d, x, b_timese, df_basis){ # d: number-th of basis
  # x: value corresponding 
  # b_timese: type of basis used
  # we should confirm that b_timese is a valid option in this function 
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  
  if(b_timese == "Legen"){ # https://github.com/xcding1212/SIMle/blob/main/SIMle/R/SIMle.Legen.v1.R
    return(Legendre_basis(d, x))
    
  }else if(b_timese == "Cheby"){ # https://github.com/xcding1212/SIMle/blob/main/SIMle/R/SIMle.Cheby.v1.R
    return(Chebyshev_basis(d, x))
    
  }else if(b_timese %in% c("tri", "cos", "sin")){ # https://github.com/xcding1212/SIMle/blob/main/SIMle/R/SIMle.Four.v1.R
    return(select.basis(d, x, ops = b_timese))
    
  }else if(b_timese == "Cspli"){ # https://github.com/xcding1212/SIMle/blob/main/SIMle/R/SIMle.Csp.v1.R
    return(Spline_basis(d, x, df_basis)) # this d is c for spline and the true number of basis is c-2+or
    
  }else{ # for wavelet # https://github.com/xcding1212/SIMle/blob/main/SIMle/R/SIMle.db1-20.v1.R
    return(Wavelet_basis(d, x, df_basis))
  }
}






db_table = function(w, ops){ # c indicate the order of db and w indicate the number of basis(the true basis is 2^w).
  if(ops == "db1"){
    point = 10000
    dbt = valdb(w,psi.f = 0, 1, point)
    n = dim(dbt)[1]
    x = seq(0, 1, length=n)
    
    df = matrix(nrow = n) 
    for (i in 2:(2^w + 1)){
      res = dbt[, i] #unlist(poly_val(coeffi, i))
      df = cbind(df, res)
    }
    df = df[,-1]   
    return(df)
  } else{
    #library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }
      
    } else{
      c = as.numeric(aux_str[3])*3
    }
    
    # library(RCurl)
    ### e.g. https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/db2_fa_table
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    
    dbt = dbt[,2]
    n = length(dbt)
    df = matrix(nrow = n) 
    t=seq(0, 2*c - 1, length.out = length(dbt))
    n = length(t)
    x =seq(0,1, length.out = n)
    for(h in 0:(2^w-1)){
      p=rep(0, n)
      for (k in 1:n){
        for (l in -70:70){
          p[k]=p[k]+ w_find(dbt, t, 2*c - 1, (2^(w)*(x[k]+l)-h))
        }
        p[k]=2^(w/2)*p[k]
      }
      
      df= cbind(df, p)
      p=rep(0, n)
    }
    df = df[,-1]
    return(df)
  }
}


Cspline_table = function(c, n=1000, or = 4){
  #library(splines)
  x = seq(0, 1, length=n)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ## need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)                ## default order: ord =4, corresponds to cubic splines
  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
  }
  
  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
  }
  B = csptable[, -1]
  return(B)
}



# c is number of basis function, x is the inputs value. 
Spline_basis <- function(k, x_val, df){ # true k is c-2+or.    k = c-2+or   
  n = dim(df)[1]
  x = seq(0,1, length.out = n)
  return(df[unique(which(abs(x - x_val) == min(abs(x -x_val)))),])
}



# Chebyshev basis function (kind 1)
# input is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order. The out put is this polynomial multiple by x.
move_order = function(poly){ # get the coefficient when polynomial multiple 1 times of x
  return(c(0,poly))
}

# we can get the plot from the coefficients
# input is n (the number of basis function)
# output is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.


# No normalized, kind only can be taken 1 and 2.
# The Chebyshev polynomials of the first kind are a special case of the Jacobi polynomials where alpha=beta=-1/2


Chebyshev_coeff = function(n, kind = 1){
  
  p_c = list()
  p_c[1] = c(1)
  if(n == 1){
    return(p_c[[1]])
  } else if(n==2){
    if(kind == 1){
      p_c[[2]] = c(0,1)
      return(p_c)
    } else{
      p_c[[2]] = c(0,2)
      return(p_c)
    }
    
  } else{
    if (kind == 1){
      for(i in 3:n){
        aux.i = i - 1
        p_c[[2]] = c(0,1)
        p_c[[i]] = 2*move_order(p_c[[i-1]]) - c(p_c[[i-2]], 0, 0)
      }
    } else{
      for(i in 3:n){
        aux.i = i - 1
        p_c[[2]] = c(0,2)
        p_c[[i]] = 2*move_order(p_c[[i-1]]) - c(p_c[[i-2]], 0, 0)
      }
    }
  }
  return(p_c)
  
}


# input is n (the number of basis function). coeffi is the list and contain the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
poly_val = function(coeffi, x){
  for(i in 1:length(coeffi)){
    aux = 0
    for(j in 1:length(coeffi[[i]])){
      aux = aux + coeffi[[i]][j]*(2*x-1)^(j - 1)
    }
    coeffi[[i]] = aux
  }
  return(coeffi)
}

# c is number of basis function, x is the inputs value.
Chebyshev_basis = function(c, x, kind = 1){
  
  aux_li = Chebyshev_coeff(c, kind)
  res = c()
  res.aux = unlist(poly_val(aux_li, x))
  normcos = c()
  
  for(i in 1:c){
    fx = function(x){
      aux = 0
      len = length(Chebyshev_coeff(c)[[i]])
      for(j in 1:len){
        aux = aux + Chebyshev_coeff(c)[[i]][j]*((2*x-1)^(j-1))
      }
      return(aux^2)
    }
    normcos[i] = sqrt(1/integrate(fx, 0, 1)[[1]])
  }
  
  for(n in 1:c){
    res[n] = res.aux[n]*normcos[n]
  }
  return(res)
}



# Cubic-spline Demo
# Write function for C-spline, # c is the number of the cubic spline.In detail, splines library is called.

## ord options default 4   number of basis c-2+r
# ord = 5 we need add 4 additional points

# c is at least 2, which is the number of the basis function - 2, the true number of basis is c+2

phi_csp <- function(cspline, beta, b){
  c = length(cspline)
  b_res = list()
  for(i in 0:b){
    B.aux = matrix(c(rep(0, c*i), cspline, rep(0, c*(b-i))), ncol = 1)
    b_res[[i+1]] = as.numeric(t(beta)%*%B.aux)
  }
  
  return(b_res)
}

csp_basis.f <- function(b.table, x){
  return(as.numeric(b.table[which(abs(b.table[,1]-x) == min(abs(b.table[,1]-x)))[1], -1]))
}

Cspline_table = function(c, n=1000, or = 4){
  #library(splines)
  x = seq(0, 1, length=n)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ## need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)                ## default order: ord =4, corresponds to cubic splines
  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
  }
  
  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
  }
  B = csptable[, -1]
  return(B)
}



# c is number of basis function, x is the inputs value. 
Spline_basis <- function(k, x_val, df){ # true k is c-2+or.    k = c-2+or   
  n = dim(df)[1]
  x = seq(0,1, length.out = n)
  return(df[unique(which(abs(x - x_val) == min(abs(x -x_val)))),])
}


# c indicate the number of tri basis function. The true number of basis functions are 2*c - 1. For the rest basis, c expresses the number of basis function.
# Compared with version1, F_general() is removed

F_tri = function(c, x){
  aux = c()
  if(c == 1){
    aux[c] = 1
    return(aux)
  }
  ind = 1
  aux[ind] = 1
  ind = ind + 1
  for(i in 1:(c-1)){
    aux[ind] = sqrt(2)*sin(2*i*pi*x)
    aux[ind+1] = sqrt(2)*cos(2*i*pi*x)
    ind = ind + 2
  }
  return(aux)
  
}

F_cosPol = function(c, x){
  aux = c()
  if(c == 1){
    aux[c] = 1
    return(aux)
  }
  for(i in 1:c){
    if(i == 1){
      aux[i] = 1
    } else{
      aux[i] = sqrt(2)*cos((i-1)*pi*x)
    }
  }
  return(aux)
  
}

F_sinPol = function(c, x){
  aux = c()
  for(i in 1:c){
    aux[i] = sqrt(2)*sin(i*pi*x)
    
  }
  return(aux)
}

# F_general = function(c, x){
#   aux = c()
#   for(i in 1:c){
#     if (i == 1){
#       aux[i] = 1
#     } else if (i %% 2 == 0){
#       aux[i] = sqrt(2)*sin(i*pi*x)
#     } else{
#       aux[i] = sqrt(2)*cos((i-1)*pi*x)
#     }
#   }
#   return(aux)
# }

# If the option parameter equals tri, it means we choose trigometric basis, cos means cospol, sin means sinpol.
select.basis = function(c,x, ops = "tri"){
  if(ops == "tri"){
    return(F_tri(c, x))
  } else if (ops == "cos"){
    return(F_cosPol(c,x))
  } else{
    return(F_sinPol(c,x))
  }
}


# tri --> it gives you (2*d - 1) basis 




# This is the first function in the SIMle package and here are several tasks remain. 

# 
# 2. Created guideline 
# 3. Finish function to SIMle

# Remaining: 
# 1. Addiing Cspline and wavelet for basis 
# 2. Change matrix structure, now only allow r =1 but r can be any number. 
# 4. Not including real data.
# 1. There are 4 mapped functions, x should close to 0 to apply this approach.




# Legendre
move_order <- function(poly){ # get the coefficient when polynomial multiple 1 times of x
  return(c(0,poly))
}
# we can get the plot from the coefficients
# input is n (the number of basis function)
# output is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
# always normalized
legendre_coeff <- function(n){
  p_c = list()
  p_c[1] = c(1)
  if(n == 1){
    return(p_c[[1]])
  } else if(n==2){
    p_c[[2]] = c(0,sqrt(3))
    return(p_c)
    
  } else{
    for(i in 3:n){
      aux.i = i - 1
      p_c[[2]] = c(0, 1)
      p_c[[i]] = ((2*aux.i-1)/aux.i)*move_order(p_c[[i-1]]) - ((aux.i-1)/aux.i)*c(p_c[[i-2]], 0, 0)
    }
    p_n = list(p_c[[1]])
    for (i in 2:length(p_c)){
      p_n[[i]] = sqrt((2*i-1))*p_c[[i]]
    }
    # p_n[[1]] = 1/sqrt(2)
    return(p_n)
  }
}
# x \in [0,1]
# input is n (the number of basis function). coeffi is the list and contain the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
poly_val <- function(coeffi, x){
  for(i in 1:length(coeffi)){
    aux = 0
    for(j in 1:length(coeffi[[i]])){
      aux = aux + coeffi[[i]][j]*(2*x-1)^(j - 1)
    }
    coeffi[[i]] = aux
  }
  return(coeffi)
}
# c is number of basis function, x is the inputs value.
Legendre_basis <- function(c, x){
  aux_li = legendre_coeff(c)
  return(unlist(poly_val(aux_li, x)))
}



# Daubechies Wavelet 1 and some general functions to be used for Daubechies, the basis of wavelet is formulated by the method of Meyer (S.2 in the paper)
# However, we also try the S.1 for estimating coefficients and choose the J_0 = 0.
#library(ggplot2)


# find the minimal value correspond, inter is the interval created, valtable is the db2 value table which is generated in advance.upper is the upper bound of basis. val express the input x.
w_find = function(valtable, inter, upper, val){
  if(val < 0 | val >= upper){
    return(0)
  } else{
    return(valtable[which(abs(inter-val) == min(abs(inter-val)))])
  }
}

db1_f = function(t){
  return(ifelse(t<1 & t>=0, 1,0))
}


wavelet_kth_b = function(k, point, ops){ # the k-th basis in 2^k total basis, n refers number of points
  w = k
  if(ops == "db1"){
    dbt = valdb(w,psi.f = 0, 1, point)
    n = dim(dbt)[1]
    x = seq(0, 1, length = n)
    
    df = data.frame(x)
    for (i in 2:(2^w + 1)){
      res = as.data.frame(dbt[, i]) #unlist(poly_val(coeffi, i))
      df = cbind(df, res)
    }
    return(data.matrix(df[,-1]))
  } else{
    #library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }
      
    } else{
      c = as.numeric(aux_str[3])*3
    }
    
    #library(RCurl)
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    dbt = dbt[,2]
    t=seq(0, 2*c - 1, length.out = length(dbt))
    n = length(t)
    x =seq(0,1, length.out = n)
    df = data.frame(x)
    for(h in 0:(2^w-1)){
      p=rep(0, n)
      for (k in 1:n){
        for (l in -70:70){
          p[k]=p[k]+ w_find(dbt, t, 2*c - 1, (2^(w)*(x[k]+l)-h))
        }
        p[k]=2^(w/2)*p[k]
      }
      
      df=cbind(df, data.frame(value = p))
      p=rep(0, n)
    }
    db.leng = seq(0, 1, length.out = point)
    new_x = c()
    for (i in 1:point){
      new_x[i] = unique(which(abs(df[,1] - db.leng[i]) == min(abs(df[,1] - db.leng[i]))))
    }
    df = df[new_x, -1]
    return(data.matrix(df))
  }
}


# res = wavelet_kth_b(1, 2000, "db3")



# This basis is formulated by Meyer for Daubechies1
valdb1 = function(w, n){
  x =seq(0,1, length.out = n)
  df.db1 = data.frame(x = x)
  
  for(h in 0:(2^w-1)){
    p=rep(0, n);
    for (k in 1:n){
      for (l in 0:70){
        p[k]=p[k]+ db1_f((2^(w)*(x[k]+l)-h))
      }
      p[k]=2^(w/2)*p[k]
    }
    df.db1[,h+2] = p
  }
  return(df.db1)
}


# get the res with different w.

# db_number represent which order of Daubechies to be used. 1-20 right now could be chosen. w indicates the number of basis functions.

valdb = function(w, psi.f=0, db_number, len.n){ # len.n indicates the point chosen
  if(db_number == 1){
    return(valdb1(w, len.n))
  } else{
    n = length(psi.f)
    x =seq(0,1, length.out = n)
    df.db = data.frame(x)
    t=seq(0, 2*db_number - 1, length.out = n)
    
    for(k in 0:(2^w-1)){
      p=rep(0, n);
      for (ind in 1:n){
        for (l in 0:70){
          p[ind]=p[ind]+ w_find(psi.f, t, 2*db_number - 1, (2^(w)*(x[ind]+l)-k))
        }
        p[ind]=(2^(w/2))*p[ind]
      }
      df.db[,k+2] = p
    }
    return(df.db)
  }
  
}


library(stringr)
library(RCurl)

db_table = function(w, ops){ # c indicate the order of db and w indicate the number of basis(the true basis is 2^w).
  if(ops == "db1"){
    point = 10000
    dbt = valdb(w,psi.f = 0, 1, point)
    n = dim(dbt)[1]
    x = seq(0, 1, length=n)
    
    df = matrix(nrow = n) 
    for (i in 2:(2^w + 1)){
      res = dbt[, i] #unlist(poly_val(coeffi, i))
      df = cbind(df, res)
    }
    df = df[,-1]
    return(df)
  } else{
    #library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }
      
    } else{
      c = as.numeric(aux_str[3])*3
    }
    
    # library(RCurl)
    
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    
    dbt = dbt[,2]
    n = length(dbt)
    df = matrix(nrow = n) 
    t=seq(0, 2*c - 1, length.out = length(dbt))
    n = length(t)
    x =seq(0,1, length.out = n)
    for(h in 0:(2^w-1)){
      p=rep(0, n)
      for (k in 1:n){
        for (l in -70:70){
          p[k]=p[k]+ w_find(dbt, t, 2*c - 1, (2^(w)*(x[k]+l)-h))
        }
        p[k]=2^(w/2)*p[k]
      }
      
      df= cbind(df, p)
      p=rep(0, n)
    }
    df = df[,-1]
    return(df)
  }
}



# c is number of basis function, x is the inputs value. 
Wavelet_basis <- function(k, x_val, df){ # true k is c-2+or.    k = c-2+or   
  n = dim(df)[1]
  x = seq(0,1, length.out = n)
  return(df[unique(which(abs(x - x_val) == min(abs(x -x_val)))),])
}



bs.gene = function(type, k, point = 200, c=10, or = 4, ops = "auto"){
  # library(ggplot2)
  # Bspline and Cspline indicate the k-th basis under the total k+1 basis. Point number is fixed for wavelet.
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )
  
  if(ops == "non-auto"){
    return(cat("general table first"))
    
  } else{
    if(type == "Legen"){
      return(Legendre_kth_b(k, point))
    } else if (type == "Cheby"){
      return(Chebyshev_kth_b(k, point))
    } else if (type %in% c("tri", "cos", "sin")){
      return(Fourier_kth_b(k, point, ops = type))
    } else if (type == "Cspli"){
      return(Cspline_kth_b(k, point, c, or=or))
    } else if (type %in% wavelet_basis){
      return(wavelet_kth_b(k, ops = type))
    } else{
      return(stop("Invalid option!"))
    }
  }
  
}



# c indicate the number of tri basis function. The true number of basis functions are 2*c - 1. For the rest basis, c expresses the number of basis function.
# Compared with version1, F_general() is removed

F_tri = function(c, x){
  aux = c()
  if(c == 1){
    aux[c] = 1
    return(aux)
  }
  ind = 1
  aux[ind] = 1
  ind = ind + 1
  for(i in 1:(c-1)){
    aux[ind] = sqrt(2)*sin(2*i*pi*x)
    aux[ind+1] = sqrt(2)*cos(2*i*pi*x)
    ind = ind + 2
  }
  return(aux)
  
}

F_cosPol = function(c, x){
  aux = c()
  if(c == 1){
    aux[c] = 1
    return(aux)
  }
  for(i in 1:c){
    if(i == 1){
      aux[i] = 1
    } else{
      aux[i] = sqrt(2)*cos((i-1)*pi*x)
    }
  }
  return(aux)
  
}

F_sinPol = function(c, x){
  aux = c()
  for(i in 1:c){
    aux[i] = sqrt(2)*sin(i*pi*x)
    
  }
  return(aux)
}

# F_general = function(c, x){
#   aux = c()
#   for(i in 1:c){
#     if (i == 1){
#       aux[i] = 1
#     } else if (i %% 2 == 0){
#       aux[i] = sqrt(2)*sin(i*pi*x)
#     } else{
#       aux[i] = sqrt(2)*cos((i-1)*pi*x)
#     }
#   }
#   return(aux)
# }

# If the option parameter equals tri, it means we choose trigometric basis, cos means cospol, sin means sinpol.
select.basis = function(c,x, ops = "tri"){
  if(ops == "tri"){
    return(F_tri(c, x))
  } else if (ops == "cos"){
    return(F_cosPol(c,x))
  } else{
    return(F_sinPol(c,x))
  }
}


Fourier_kth_b = function(k, n, ops){
  df = data.frame()
  aux_x = seq(0,1, length.out = n)
  # coeffi = legendre_coeff(k)[k]
  for (i in aux_x){
    res = select.basis(k, i, ops)
    df = rbind(df, res)
  }
  df = data.frame(basis_value = df[,k])
  return(df)
}


fourier_plot = function(c, ops = "tri", title){
  
  df = data.frame()
  aux_x = seq(0,1,0.005)
  
  for (i in aux_x){
    res = select.basis(c, i, ops)
    df = rbind(df, res)
  }
  
  new = c()
  for(i in 1:dim(df)[2]){
    new = c(new, df[, i])
  }
  f.df = as.data.frame(new)
  f.df$x = rep(aux_x, dim(df)[2])
  f.df$order = as.factor(rep(0:(dim(df)[2]-1), each = 201))
  theme_update(plot.title = element_text(hjust = 0.5))
  p1 <- ggplot(f.df, aes(x=x, y=new, group=order, colour = order))+ geom_line()  + ggtitle(title) +
    xlab("") + ylab("") + scale_colour_discrete(name  ="order")+theme(plot.title = element_text(size=18, face="bold"),
                                                                      legend.text=element_text(size=24, face = "bold"),
                                                                      axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                      axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                      axis.title.x=element_text(size=22,face='bold'),
                                                                      axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                      legend.title = element_text(face = "bold"))
  return(p1)
}

# fourier_plot(2)

# In the function beta_f(), alpha.fou(), basis depends on the options, we have 4 choice for fourier basis.
# beta_f returns the list which contains beta and design matrix Y

beta_f = function(ts, c, b, ops = "tri"){
  n = length(ts)
  X = matrix(ts[(b+1):n], ncol = 1)
  aux = c(select.basis(c, (b+1)/n, ops))
  for(j in 1:b){
    aux = c(aux, select.basis(c, ((b+1)/n), ops)*ts[b+1-j])
  }
  Y = matrix(aux, nrow = 1) # i =2 j >= 1
  
  for(i in (b+2):n){
    aux = c(select.basis(c, i/n, ops))
    for(j in 1:b){
      aux = c(aux, select.basis(c, (i/n),ops)*ts[i-j])
    }
    
    aux_Y = matrix(aux, nrow = 1)
    Y = rbind(Y, aux_Y)
  }
  beta = solve(t(Y)%*%Y, tol = 1e-40)%*%t(Y)%*%X
  return(list(beta, Y))
}

phi_f = function(fourier, beta, b){
  c = length(fourier)
  b_res = list()
  for(i in 0:b){
    B.aux = matrix(c(rep(0, c*i), fourier, rep(0, c*(b-i))), ncol = 1)
    b_res[[i+1]] = as.numeric(t(beta)%*%B.aux)
  }
  return(b_res)
  
}


# we want to generate m points of coefficients of time series and the default number is 500.
alpha.fou = function(ts, c, b, m=500, ops){
  
  l.alpha = list()
  aux.alpha = c()
  beta.es = beta_f(ts, c, b, ops)
  
  for(j in 1:(b+1)){
    for(i in 1:m){
      aux.alpha[i] = phi_f(select.basis(c, i/m, ops), beta.es[[1]], b)[[j]]
    }
    l.alpha[[j]] = aux.alpha
    aux.alpha = c()
  }
  return(l.alpha)
}

alpha.loocv.f = function(ts, c, b, ops){
  n = length(ts)
  aux.true = ts[(b+1):n]
  aux.esti = c()
  leve.i = c()
  
  beta.es = beta_f(ts, c, b, ops)
  hat = beta.es[[2]]%*%solve(t(beta.es[[2]])%*%beta.es[[2]])%*%t(beta.es[[2]])
  
  for(i in (b+1):n){
    aux.esti[i-b]  = matrix(beta.es[[2]][i-b,], nrow = 1)%*%beta.es[[1]]
    leve.i[i-b] = as.numeric(hat[i-b,i-b]) # hii is the diagonal of the hat matrix
  }
  error = (aux.true - aux.esti)^2
  lever = sum(error/((1-leve.i)^2))/(n-b)
  return(c(c,b,lever))
}

#  CV in the paper
alpha.cv.f = function(ts, c, b, ops){
  n = length(ts)
  l = floor(3*log2(n))
  aux.train = ts[1:(n-l)]
  aux.vali = ts[(n-l+1):n]
  tt = fix.fit.four(aux.train, c, b, length(aux.train), ops)
  pre = predict.Four(aux.train, tt, length(aux.vali))
  error = sum((aux.vali - pre)^2)/l
  return(c(c,b,error))
}


# prediction
predict.Four = function(ts, esti.li, k){ # k indicates the number of predictions
  ts.pre = c()
  phi.h = esti.li[[2]]
  n = length(ts)
  b = length(phi.h)
  
  for(h in 1:k){
    aux.pre = phi.h[[1]][n]
    for(j in 2:b){
      aux.pre = aux.pre + phi.h[[j]][n]*ts[n-h-j]
    }
    ts.pre[h] = aux.pre
  }
  return(ts.pre)
}


# The return of fit.ts.f() is the list contains 4 parts named Estimate, cv, coefficients and bc. Rather, Estimate is the estimate of coefficient for the time series
# cv is the cross validation matrix, Coefficients is the estimate of coefficient for each basis function, BC contains the best number of b and c basis on LOOCV method.
# fit.ts.f() automatically choose the best b and c for the time series and get the estimate basis on that.


fix.fit.four = function(ts, c, b, m, ops){
  
  error.s = c()
  n = length(ts)
  es.alpha = alpha.fou(ts, c, b, n, ops)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*ts[i-j+1]
    }
    error.s[i-b] = ts[i] - val.aux
  }
  
  return(list(ols.coef = beta_f(ts, c, b, ops)[[1]], ts.coef = alpha.fou(ts, c, b, m, ops), Residuals = error.s))
  
}

auto.fit.four = function(ts, c = 10, b = 3, m = 500, ops, method = "LOOCV", threshold = 0){
  res.bc = matrix(ncol = 3, nrow = c*b)
  ind = 1
  for(i in 1:c){
    for(j in 1:b){
      if(method == "CV"){
        res.bc[ind, ] = alpha.cv.f(ts, i, j, ops)
      } else{
        res.bc[ind, ] = alpha.loocv.f(ts, i, j, ops)
      }
      ind = ind + 1
    }
  }
  
  colnames(res.bc) = c("c", "b", "cv")
  
  if(method == "Elbow"){
    
    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    res.bc = res.bc[which(res.bc[,2] == b.s), ]
    
    if(threshold == 0){
      c.s = 1 + which(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1) == max(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1)))
    } else{
      c.s = max(which(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1) >= threshold)) + 1
    }
    estimate = alpha.fou(ts, c.s, b.s, m, ops)
    
  }else{
    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    c.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),1]
    estimate = alpha.fou(ts, c.s, b.s, m, ops)
  }
  
  
  return(list(Estimate = estimate, CV = res.bc, Coefficients = beta_f(ts, c.s, b.s)[[1]], BC = c(c.s, b.s)))
}

# Testing
mv_method.four = function(timese, c, b, ops){
  h.0 = 3
  m.li = c(1:25)
  #library(Matrix)
  # Design matrix
  Y = beta_l(timese, c, b)[[2]]
  n = length(timese)
  # li.res = list()
  # m = 6
  
  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.fou(timese, c, b, n, ops)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }
  
  Phi.li = list()
  for (m in m.li){
    aux_Phi=0
    Phi = 0
    for(i in (b+1):(n-m)){
      h = 0
      for(j in i:(i+m)){
        aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
        h = h + aux.h
      }
      B = matrix(select.basis(c, i/n, ops), ncol = 1)
      Phi = Phi + kronecker(h, B)
      aux_Phi = aux_Phi + Phi%*%t(Phi)
    }
    Phi.li[[m]] = 1/((n-m-b+1)*m)*aux_Phi
  }
  
  se.li = list()
  for(mj in (min(m.li)+h.0):(max(m.li)-h.0)){
    av.Phi = 0
    se = 0
    for (k in -3:3){
      av.Phi = av.Phi + Phi.li[[mj + k]]
    }
    av.Phi = av.Phi/7
    
    for(k in -3:3){
      se = se + norm(av.Phi - Phi.li[[mj + k]], "2")^2
    }
    se.li[[mj-3]] = sqrt(se/6)
  }
  return(m.op = which(unlist(se.li) == min(unlist(se.li))) + 3)
  # return(unlist(se.li))
}


fix.test.four = function(timese, c, b, ops, B.s, m){
  #library(Matrix)
  # Design matrix
  Y = beta_f(timese, c, b, ops)[[2]]
  n = length(timese)
  # li.res = list()
  # m = 6
  if(m == 0){
    m = mv_method.four(timese, c, b, ops) #mv_method(timese, c, b)  #floor(n^(1/3))
  }
  
  esti = alpha.fou(timese, c, b, 10000, ops) # the estimate of coefficients
  
  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.fou(timese, c, b, n, ops)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }
  length(error.s)
  
  
  # B
  inte = select.basis(c, 1/10000, ops)*(1/10000)
  for(i in 2:10000){
    inte = inte + select.basis(c, i/10000, ops)*(1/10000)
  }
  if(ops == "tri"){
    r.c = 2*(c-1)+1
    
  }else{
    r.c = c # 2*(c-1)+1 only for tri basis.
  }
  
  # I.bc
  I = matrix(rep(0, ((b+1)*r.c)^2), ncol = (b+1)*r.c)
  for(ind in 0:(b*r.c-1)){
    I[dim(I)[1] - ind, dim(I)[2] -ind] = 1
  }
  
  nT = 0
  for (i in 2:length(esti)){
    nT = nT + sum(((esti[[i]] - sum(esti[[i]]/10000))^2)/10000)
  }
  nT = n*nT
  
  Sigma = n*solve(t(Y)%*%Y, tol = 1e-40)
  inte = matrix(inte, ncol =1)
  W = diag(r.c) - inte%*%t(inte)
  W = matrix(bdiag(replicate(b+1,W,simplify=FALSE)), ncol = (b+1)*r.c)
  
  Tao = Sigma%*%I%*%W%*%Sigma
  
  
  # hist(unlist(Sta))
  # print(ite)
  
  Sta = list()
  Phi.li = list()
  
  for(k in 1:B.s){
    R = rnorm(n-m-b, 0, 1)
    Phi = 0
    for(i in (b+1):(n-m)){
      h = 0
      for(j in i:(i+m)){
        aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
        h = h + aux.h
      }
      B = matrix(select.basis(c, i/n, ops), ncol = 1)
      Phi = Phi + kronecker(h, B)*R[i-b]
    }
    
    Phi = (1/sqrt((n-m-b+1)*m))*Phi
    Phi.li[[k]] = Phi
  }
  # image(W)
  # W[(c+1):dim(W)[1], (c+1):dim(W)[2]] = 0
  for(k in 1:B.s){
    Sta[[k]] = t(Phi.li[[k]])%*%Tao%*%Phi.li[[k]]
  }
  
  # nT > sort(unlist(Sta))[950] if TRUE reject the null
  return(1 - sum(unlist(Sta) <= nT)/B.s) # P value
  
}



# testing b


fit.testing.b.four = function(timese, c, b.0 = 3, ops, b = 8, B.s, m){
  if(b.0 >= b){return(FALSE)}
  #library(Matrix)
  # Design matrix
  Y = beta_f(timese, c, b, ops)[[2]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.four(timese, c, b, ops) # does m influenced by b ??? mv_method(timese, c, b)  #floor(n^(1/3))
  }
  
  esti = alpha.fou(timese, c, b, 10000, ops) # the estimate of coefficients
  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.fou(timese, c, b, n, ops)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }
  if(ops == "tri"){
    r.c = 2*(c-1)+1
    
  }else{
    r.c = c # 2*(c-1)+1 only for tri basis.
  }
  
  aux.pval = list()
  
  
  # B
  for(k.aux in 0:(b.0-1)){ # 0 ---> 1-15
    # 1 ---> 2-15
    nT = 0
    for (i in (2+k.aux):length(esti)){
      nT = nT + sum((esti[[i]]^2)/10000)
    }
    nT = n*nT
    
    # I.bc
    I = matrix(rep(0, ((b+1)*r.c)^2), ncol = (b+1)*r.c)
    for(ind in 0:((b-k.aux)*r.c-1)){
      I[dim(I)[1] - ind, dim(I)[2] -ind] = 1
    }
    
    Sigma = n*solve(t(Y)%*%Y, tol = 1e-40)
    Tao = Sigma%*%I%*%Sigma
    
    Sta = list()
    Phi.li = list()
    
    for(k in 1:B.s){
      R = rnorm(n-m-b, 0, 1)
      Phi = 0
      for(i in (b+1):(n-m)){
        h = 0
        for(j in i:(i+m)){
          aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
          h = h + aux.h
        }
        B = matrix(select.basis(c, i/n, ops), ncol = 1)
        Phi = Phi + kronecker(h, B)*R[i-b]
      }
      
      Phi = (1/sqrt((n-m-b+1)*m))*Phi
      Phi.li[[k]] = Phi
    }
    for(k.aux2 in 1:B.s){
      Sta[[k.aux2]] = t(Phi.li[[k.aux2]])%*%Tao%*%Phi.li[[k.aux2]]
      
    }
    aux.pval[[k.aux+1]] = 1 - sum(unlist(Sta) <= nT)/B.s
  }
  # nT > sort(unlist(Sta))[950] if TRUE reject the null
  return(aux.pval)
}


# Legendre basis function
# input is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order. The out put is this polynomial multiple by x.
move_order = function(poly){ # get the coefficient when polynomial multiple 1 times of x
  return(c(0,poly))
}

# we can get the plot from the coefficients
# input is n (the number of basis function)
# output is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
# always normalized

legendre_coeff = function(n){
  p_c = list()
  p_c[1] = c(1)
  if(n == 1){
    return(p_c[[1]])
  } else if(n==2){
    p_c[[2]] = c(0,sqrt(3))
    return(p_c)
    
  } else{
    for(i in 3:n){
      aux.i = i - 1
      p_c[[2]] = c(0, 1)
      p_c[[i]] = ((2*aux.i-1)/aux.i)*move_order(p_c[[i-1]]) - ((aux.i-1)/aux.i)*c(p_c[[i-2]], 0, 0)
    }
    p_n = list(p_c[[1]])
    for (i in 2:length(p_c)){
      p_n[[i]] = sqrt((2*i-1))*p_c[[i]]
    }
    # p_n[[1]] = 1/sqrt(2)
    return(p_n)
  }
}


# x \in [0,1]
# input is n (the number of basis function). coeffi is the list and contain the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
poly_val = function(coeffi, x){
  for(i in 1:length(coeffi)){
    aux = 0
    for(j in 1:length(coeffi[[i]])){
      aux = aux + coeffi[[i]][j]*(2*x-1)^(j - 1)
    }
    coeffi[[i]] = aux
  }
  return(coeffi)
}

# c is number of basis function, x is the inputs value.
Legendre_basis = function(c, x){
  aux_li = legendre_coeff(c)
  return(unlist(poly_val(aux_li, x)))
}


Legendre_kth_b = function(k, n){
  df = data.frame()
  aux_x = seq(0,1, length.out = n)
  coeffi = legendre_coeff(k)[k]
  for (i in aux_x){
    res = unlist(poly_val(coeffi, i))
    df = rbind(df, res)
  }
  colnames(df) = c("basis_value")
  return(df)
}


legendre_plot = function(n, title){
  coeffi = legendre_coeff(n)
  df = data.frame()
  aux_x = seq(0,1,0.005)
  
  for (i in aux_x){
    res = unlist(poly_val(coeffi, i))
    df = rbind(df, res)
  }
  
  new = c()
  for(i in 1:dim(df)[2]){
    new = c(new, df[, i])
  }
  f.df = as.data.frame(new)
  f.df$x = rep(aux_x, dim(df)[2])
  f.df$order = as.factor(rep(0:(dim(df)[2]-1), each = 201))
  
  theme_update(plot.title = element_text(hjust = 0.5))
  p1 <- ggplot(f.df, aes(x=x, y=new, group=order, colour = order))+ geom_line()  + ggtitle(title) +
    xlab("") + ylab("") + scale_colour_discrete(name  ="order")+theme(plot.title = element_text(size=18, face="bold"),
                                                                      legend.text=element_text(size=24, face = "bold"),
                                                                      axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                      axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                      axis.title.x=element_text(size=22,face='bold'),
                                                                      axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                      legend.title = element_text(face = "bold"))
  return(p1)
}



# ggtitle(TeX("L_1"))

# beta is the fitted estimators and Y is the design matrix. This Y is for cv.

beta_l = function(ts, c, b){
  n = length(ts)
  X = matrix(ts[(b+1):n], ncol = 1)
  aux = c(Legendre_basis(c, (b+1)/n))
  for(j in 1:b){
    aux = c(aux, Legendre_basis(c, (b+1)/n)*ts[b+1-j])
  }
  Y = matrix(aux, nrow = 1) # i =2 j >= 1
  
  for(i in (b+2):n){
    aux = c(Legendre_basis(c, i/n))
    for(j in 1:b){
      aux = c(aux, Legendre_basis(c, i/n)*ts[i-j])
    }
    
    aux_Y = matrix(aux, nrow = 1)
    Y = rbind(Y, aux_Y)
  }
  beta = solve(t(Y)%*%Y, tol = 1e-40)%*%t(Y)%*%X
  return(list(beta, Y))
}



# for each i/n, we can get the coefficient estimated.
phi_l = function(legendre, beta, b){
  c = length(legendre)
  b_res = list()
  for(i in 0:b){
    B.aux = matrix(c(rep(0, c*i), legendre, rep(0, c*(b-i))), ncol = 1)
    b_res[[i+1]] = as.numeric(t(beta)%*%B.aux)
  }
  
  return(b_res)
}


# ts is the time series data, c is the number of basis, b indicates AR(b_i)

alpha.legen = function(ts, c, b, m=500){
  
  l.alpha = list()
  aux.alpha = c()
  beta.es = beta_l(ts, c, b)
  
  for(j in 1:(b+1)){
    for(i in 1:m){
      aux.alpha[i] = phi_l(Legendre_basis(c, i/m), beta.es[[1]], b)[[j]]
    }
    l.alpha[[j]] = aux.alpha
    aux.alpha = c()
  }
  return(l.alpha)
}

# Leave one out
alpha.loocv.l = function(ts, c, b){
  n = length(ts)
  aux.true = ts[(b+1):n]
  aux.esti = c()
  leve.i = c()
  
  beta.es = beta_l(ts, c, b)
  hat = beta.es[[2]]%*%solve(t(beta.es[[2]])%*%beta.es[[2]])%*%t(beta.es[[2]])
  
  for(i in (b+1):n){
    aux.esti[i-b]  = matrix(beta.es[[2]][i-b,], nrow = 1)%*%beta.es[[1]]
    leve.i[i-b] = as.numeric(hat[i-b,i-b]) # hii is the diagonal of the hat matrix
  }
  error = (aux.true - aux.esti)^2
  lever = sum(error/((1-leve.i)^2))/(n-b)
  return(c(c,b,lever))
}

#  CV in the paper, the b incorrect when the data set is small.

alpha.cv.l = function(ts, c, b){
  n = length(ts)
  l = floor(3*log2(n))
  aux.train = ts[1:(n-l)]
  aux.vali = ts[(n-l+1):n]
  tt = fix.fit.legen(aux.train, c, b, length(aux.train))
  pre = predict.legen(aux.train, tt, length(aux.vali))
  error = sum((aux.vali - pre)^2)/l
  return(c(c,b,error))
}



# prediction

predict.legen = function(ts, esti.li, k){ # k indicates the number of predictions
  ts.pre = c()
  if(length(esti.li) == 3){
    phi.h = esti.li[[2]]
  }else{
    phi.h = esti.li[[1]]
  }
  
  n = length(phi.h[[1]])
  b = length(phi.h)
  for(h in 1:k){
    aux.pre = phi.h[[1]][n]
    for(j in 2:b){
      aux.pre = aux.pre + phi.h[[j]][n]*ts[n-h-j]
    }
    ts.pre[h] = aux.pre
  }
  return(ts.pre)
}


# the default number of c and b indicate the maximum number, fit.ts function will get the best c* and b* automatically as the returns.
# The return is the list contains estimate coefficients, the cross validation matrix, the estimate coefficients for basis functions, and the best c* and b* in order.
# ops = CV indicates cross validation and method = COOV means leave one out cross validation.(Elbow)

fix.fit.legen = function(ts, c, b, m){
  
  error.s = c()
  n = length(ts)
  es.alpha = alpha.legen(ts, c, b, n)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*ts[i-j+1]
    }
    error.s[i-b] = ts[i] - val.aux
  }
  
  return(list(ols.coef = beta_l(ts, c, b)[[1]], ts.coef = alpha.legen(ts, c, b, m), Residuals = error.s))
}
#timese = generate_AR1(786, 20)
#pp = fix.fit.legen(timese, 5, 1, 500)
auto.fit.legen = function(ts, c = 10, b = 3, m=500, method = "CV", threshold = 0){
  res.bc = matrix(ncol = 3, nrow = c*b)
  ind = 1
  for(i in 1:c){
    for(j in 1:b){
      if(method == "CV"){
        res.bc[ind, ] = alpha.cv.l(ts, i, j)
      } else{
        res.bc[ind, ] = alpha.loocv.l(ts, i, j)
      }
      
      ind = ind + 1
    }
  }
  colnames(res.bc) = c("c", "b", "cv")
  if(method == "Elbow"){
    
    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    res.bc = res.bc[which(res.bc[,2] == b.s), ]
    
    if(threshold == 0){
      c.s = 1 + which(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1) == max(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1)))
    } else{
      c.s = max(which(abs(res.bc[1:(length(res.bc[,3])-1),3]/res.bc[-1,3] - 1) >= threshold)) + 1
    }
    estimate = alpha.legen(ts, c.s, b.s, m)
    
  }else{
    b.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),2]
    c.s = res.bc[which(res.bc[,3] == min(res.bc[, 3])),1]
    estimate = alpha.legen(ts, c.s, b.s, m)
  }
  
  return(list(Estimate = estimate, CV = res.bc, Coefficients = beta_l(ts, c.s, b.s)[[1]], BC = c(c.s, b.s)))
}

# Testing
mv_method.legen = function(timese, c, b){
  h.0 = 3
  m.li = c(1:25)
  #library(Matrix)
  # Design matrix
  Y = beta_l(timese, c, b)[[2]]
  n = length(timese)
  # li.res = list()
  # m = 6
  
  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.legen(timese, c, b, n)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }
  
  Phi.li = list()
  for (m in m.li){
    aux_Phi=0
    Phi = 0
    for(i in (b+1):(n-m)){
      h = 0
      for(j in i:(i+m)){
        aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
        h = h + aux.h
      }
      B = matrix(Legendre_basis(c, i/n), ncol = 1)
      Phi = Phi + kronecker(h, B)
      aux_Phi = aux_Phi + Phi%*%t(Phi)
    }
    Phi.li[[m]] = 1/((n-m-b+1)*m)*aux_Phi
  }
  
  se.li = list()
  for(mj in (min(m.li)+h.0):(max(m.li)-h.0)){
    av.Phi = 0
    se = 0
    for (k in -3:3){
      av.Phi = av.Phi + Phi.li[[mj + k]]
    }
    av.Phi = av.Phi/7
    
    for(k in -3:3){
      se = se + norm(av.Phi - Phi.li[[mj + k]], "2")^2
    }
    se.li[[mj-3]] = sqrt(se/6)
  }
  return(m.op = which(unlist(se.li) == min(unlist(se.li))) + 3)
  # return(unlist(se.li))
}


fix.test.legen = function(timese, c, b, B.s, m){
  #library(Matrix)
  # Design matrix
  Y = beta_l(timese, c, b)[[2]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.legen(timese, c, b) #mv_method(timese, c, b)  #floor(n^(1/3))
  }
  
  # m = 6
  esti = alpha.legen(timese, c, b, 10000) # the estimate of coefficients
  
  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.legen(timese, c, b, n)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }
  length(error.s)
  
  
  # B
  inte = Legendre_basis(c, 1/10000)*(1/10000)
  for(i in 2:10000){
    inte = inte + Legendre_basis(c, i/10000)*(1/10000)
  }
  r.c = c # 2*(c-1)+1 only for tri basis.
  # I.bc
  I = matrix(rep(0, ((b+1)*r.c)^2), ncol = (b+1)*r.c)
  for(ind in 0:(b*r.c-1)){
    I[dim(I)[1] - ind, dim(I)[2] -ind] = 1
  }
  
  nT = 0
  for (i in 2:length(esti)){
    nT = nT + sum(((esti[[i]] - sum(esti[[i]]/10000))^2)/10000)
  }
  nT = n*nT
  
  Sigma = n*solve(t(Y)%*%Y, tol = 1e-40)
  inte = matrix(inte, ncol =1)
  W = diag(r.c) - inte%*%t(inte)
  W = matrix(bdiag(replicate(b+1,W,simplify=FALSE)), ncol = (b+1)*r.c)
  
  Tao = Sigma%*%I%*%W%*%Sigma
  
  
  # hist(unlist(Sta))
  # print(ite)
  
  Sta = list()
  Phi.li = list()
  
  for(k in 1:B.s){
    R = rnorm(n-m-b, 0, 1)
    Phi = 0
    for(i in (b+1):(n-m)){
      h = 0
      for(j in i:(i+m)){
        aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
        h = h + aux.h
      }
      B = matrix(Legendre_basis(c, i/n), ncol = 1)
      Phi = Phi + kronecker(h, B)*R[i-b]
    }
    
    Phi = (1/sqrt((n-m-b+1)*m))*Phi
    Phi.li[[k]] = Phi
  }
  # image(W)
  # W[(c+1):dim(W)[1], (c+1):dim(W)[2]] = 0
  for(k in 1:B.s){
    Sta[[k]] = t(Phi.li[[k]])%*%Tao%*%Phi.li[[k]]
  }
  
  # nT > sort(unlist(Sta))[950] if TRUE reject the null
  return(1 - sum(unlist(Sta) <= nT)/B.s) # P value
  
}




# Testing b

fit.testing.b.legen = function(timese, c, b.0 = 3, b = 8, B.s, m){
  # lag = 7, b = 10, we calculate, 1-10(reject), 2-10(not reject, statistics is small), 3-10 until 7-10 (not reject, statistics is small),
  if(b.0 >= b){return(FALSE)}
  #library(Matrix)
  # Design matrix
  Y = beta_l(timese, c, b)[[2]]
  n = length(timese)
  # li.res = list()
  if(m == 0){
    m = mv_method.legen(timese, c, b) #  mv_method(timese, c, b)  #floor(n^(1/3))
  }
  
  esti = alpha.legen(timese, c, b, 10000) # the estimate of coefficients
  
  # Error, i = b* + 1... n
  error.s = c()
  es.alpha = alpha.legen(timese, c, b, n)
  aux.len = length(es.alpha)
  for(i in (b+1):n){
    val.aux = es.alpha[[1]][i]
    for(j in 2:aux.len){
      val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
    }
    error.s[i-b] = timese[i] - val.aux
  }
  
  r.c = c # 2*(c-1)+1 only for tri basis.
  aux.pval = list()
  
  # testing checking.
  # B
  for(k.aux in 0:(b.0-1)){ # 0 ---> 1-15
    # 1 ---> 2-15
    nT = 0
    for (i in (2+k.aux):length(esti)){
      nT = nT + sum((esti[[i]]^2)/10000)
    }
    nT = n*nT
    # I.bc
    I = matrix(rep(0, ((b+1)*r.c)^2), ncol = (b+1)*r.c)
    for(ind in 0:((b-k.aux)*r.c-1)){
      I[dim(I)[1] - ind, dim(I)[2] -ind] = 1
    }
    
    Sigma = n*solve(t(Y)%*%Y, tol = 1e-40)
    Tao = Sigma%*%I%*%Sigma
    
    Sta = list()
    Phi.li = list()
    
    for(k in 1:B.s){
      R = rnorm(n-m-b, 0, 1)
      Phi = 0
      for(i in (b+1):(n-m)){
        h = 0
        for(j in i:(i+m)){
          aux.h = matrix(rev(c(timese[(j- b):(j - 1)],1)), ncol = 1)*error.s[j-b]
          h = h + aux.h
        }
        B = matrix(Legendre_basis(c, i/n), ncol = 1)
        Phi = Phi + kronecker(h, B)*R[i-b]
      }
      
      Phi = (1/sqrt((n-m-b+1)*m))*Phi
      Phi.li[[k]] = Phi
    }
    for(k.aux2 in 1:B.s){
      Sta[[k.aux2]] = t(Phi.li[[k.aux2]])%*%Tao%*%Phi.li[[k.aux2]]
      
    }
    #print(sort(unlist(Sta)))
    #print(nT)
    aux.pval[[k.aux+1]] = 1 - sum(unlist(Sta) <= nT)/B.s
  }
  
  # nT > sort(unlist(Sta))[950] if TRUE reject the null
  
  return(aux.pval)
}
