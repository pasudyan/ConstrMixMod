
# Required packages
req_packages <- c("mclust", "tmvtnorm", "invgamma",
                  "truncnorm", "mvnfast", "sp",
                  "matrixStats", "ggplot2", "MCMCpack",
                  "mvtnorm", "devtools")
if (!require(req_packages)){
  install.packages(req_packages)
  devtools::install_github("jakesherman/easypackages")
  library(easypackages)
  libraries(req_packages)
}

#'  Main function for inference
#'
#' This function evaluates the cumulative distribution of a univariate Gaussian distribution
#' @param mix Mixture type: one of "tmog" and "motg"
#' @param data_type one of "gauss", "beta", "bvgauss", "flowcyto", "crime"
#' @param thr vector of number of rejections per observations.
#' @return samples from posterior distribution
#' @export
#'
Infernce_ConstrMixMod <- function(mix, data_type, thr){

  # Check correct inputs
  if (!(mix %in% c("tmog", "motg"))){
    stop("Mixture type must be one of 'tmog' or 'motg'")
  }

  if (!(data_type %in% c("gauss", "beta", "bvgauss", "flowcyto", "crime"))){
    stop("Data type must be one of 'gauss', 'beta', 'bvgauss', 'flowcyto', 'crime'")
  }

  if (thr < 0){
    stop("Value of thr must be greater than 0")
  }

  # Sampling method for computation of the normalization function
  # "is" uses important sampling
  # package uses existing truncated Gaussian functions
  if (data_type == "crime"){
    samp_method <- "is"
  } else {
    samp_method <- "package"
  }


  K <- 50
  burnIn   <- 2000
  numIter  <- 5000

  # Maximum number of rejections
  ceil_rej <- 100000

  # Thinning value for the MCMC samples on the posterior
  thin_sc <- 1
  thr_l   <- length(thr)

  if (data_type == "gauss"){
    pr1 <- 0.5
    pr2 <- 0.1
    dm  <- 1
    load("sim_data/data_gauss_mean_0_500.RData")
  }

  if (data_type == "beta"){
    pr1 <- 0.1
    pr2 <- 0.1
    dm  <- 1
    load("sim_data/data_sim_beta_500.RData")
  }

  if (dm == 1){

    # Bounds for the domain
    lower <- 0
    upper <- 1

    # Function to eliminate points outside the hypercube
    constr <- function(X, u1 = lower, u2 = upper) {
      X >= u1 & X <= u2
    }

    ip_data <- ip_data_all[constr(ip_data_all)]

    num_obs <- length(ip_data)

    num_obs_sub <- sample.int(length(ip_data_all), num_obs, replace = F)
    ip_data     <- ip_data_all[num_obs_sub]
    ip_data     <- matrix(ip_data, nrow = 1, ncol = num_obs)

    # ********** Use only test data near the borders ********** #
    bound_data <- (ip_data - 0.5)*1.2 + 0.5
    bound_data_cs <- constr(bound_data)

    ind_out  <- which(bound_data_cs == FALSE)

    num_test <- 50
    test_ind <- sample(ind_out, num_test, replace = FALSE)

    # Prior parameters
    mu_pr <- mean(ip_data)
    sd_pr <- sd(ip_data)

    # Grid points for posterior predictive distribution
    delta   <- 0.005
    pts     <- seq(0, 1, delta)
  }

  if (data_type == "bvgauss"){

    # Load data
    load("sim_data/data_bvgauss_1000_2corner_modes.RData")

    # Function to eliminate points outside the hypercube
    genconst <- function(u1, u2, l1, l2) {
      in_set <- function(X) {
        if (nrow(X) == 1) X = t(X)
        X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2
      }
    }

    dm <- dim(ip_data_all)[1]
    l  <- 0
    u  <- 1

    constr  <- genconst(u, u, l, l)

    num_obs <- ncol(ip_data_all)

    smp <- sample.int(ncol(ip_data_all), num_obs, replace = F)
    ip_data <- ip_data_all[, smp]

    lower <- rep(l, dm)
    upper <- rep(u, dm)

    test_ind <- sample.int(num_obs, ceiling(0.2*num_obs), replace = F)
  }

  if (data_type == "flowcyto"){
    # Load data from library(mclust)
    # Use flowcytometry data
    data(GvHD)

    # Scale data
    sc  <- 1024
    ip_data_ctrl <- t(GvHD.control)/sc

    # Function to eliminate points outside the hypercube
    genconst <- function(u1, u2, l1, l2) {
      in_set <- function(X) {
        if (nrow(X) == 1) X = t(X)

        X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 & X[3,] > l1 &
          X[3,] < u1 & X[4,] > l2 & X[4,] < u2
      }
    }

    l1 <- 0; l2 <- 0; u  <- 1

    constr  <- genconst(u, u, l1, l2)
    num_obs <- ncol(ip_data_ctrl)
    ip_data <- ip_data_ctrl

    dm  <- dim(ip_data)[1]

    lower   <- rep(l1, dm)
    upper   <- rep(u, dm)

    # **************** Use only test data near the borders ***************** #
    bound_data <- apply(ip_data, 2, function(x) min(x) < 0.01 | max(x) > 0.99)

    ind_out  <- which(bound_data == TRUE)

    num_test <- 100
    test_ind <- sample(ind_out, num_test, replace = FALSE)
  }


  if (data_type == "crime"){

    # Load crime data and constraints for the data
    load("./crime_data/chi_murder_2012_2017.RData")
    load("./crime_data/Neigh_coord_list.RData")

    # Scaling data
    scale_crime   <- 5
    murder_loc_sc <- murder_loc_sc*scale_crime

    # Scale borders of Chicago
    coord_list_sc <- lapply(1:length(coord_list_sc_t), function(x) {
      tmp <- coord_list_sc_t[[x]]
      return(tmp*scale_crime)
    })

    num_neigh <- length(coord_list_sc)
    l_poly <- lapply(1:num_neigh, function(x) {
      temp <- Polygon(coord_list_sc[[x]])
      Polygons(list(temp), neigh_coord_df$COMMUNITY[x])
    })

    SpPol <- SpatialPolygons(l_poly)

    # Function to determine if point is in polygon
    cnstr_chicago  <- function(XX, l_poly){
      num_neigh <- length(l_poly)
      if (dim(XX)[1] == 1)
        XX <- t(XX)
      in_obs <- lapply(1:num_neigh, function(x){
        point.in.polygon(point.x = XX[1, ],
                         point.y = XX[2, ],
                         pol.x = l_poly[[x]]@Polygons[[1]]@coords[, 2],
                         pol.y = l_poly[[x]]@Polygons[[1]]@coords[, 1])
      })
      in_obs <- Reduce("+", in_obs)
      as.logical(in_obs)
    }

    constr  <- function(XX){
      cnstr_chicago(XX, l_poly)
    }

    # **************** Use only test data near the borders ***************** #
    bound_data <- t((t(murder_loc_sc) - c(1, 1.5))*1.2 + c(1, 1.5))
    bound_data_cs <- constr(t(bound_data))

    num_test <- 200
    ind_out  <- which(bound_data_cs == FALSE)
    test_ind <- sample(ind_out, num_test, replace = FALSE)

    ip_data <- t(murder_loc_sc)
    num_obs <- ncol(ip_data)

    dm <- dim(ip_data)[1]
  }

  # ======================= Main Inference ======================= #

  # Divide into test and train set with ratio 80:20
  sub <- seq(1, num_obs, 1)
  sub <- sub[-test_ind]

  train   <- ip_data[, sub]
  test    <- matrix(ip_data[, test_ind], nrow = dm)

  # Convert data to matrix
  if (dm == 1){
    train <- matrix(train, nrow = 1, ncol = length(train))
    test  <- matrix(test, nrow = 1, ncol = length(test))
  }

  # Initialize variables
  time.taken    <- c()
  test_pred     <- list()
  prod_pred_ts  <- c()
  prod_pred_pt  <- c()
  count_bins    <- list()
  rslt <- list()

  print(paste("Mix:", mix))
  print(paste("Data Type:", data_type))

  # Loop over all thresholds
  for (j in 1:length(thr)){
    print(paste("Threshold: ", thr[j]))

    # Function to run Gibbs sampling and timed
    start.time <- system.time({
      rslt[[j]] <- inf_gvhd(K, ip_data = train, axs,
                            thr = thr[j], burnIn, numIter,
                            constr, mix, lower, upper, data_type)
    })

    time.taken[j] <- start.time[3]

    # Functions to compute posterior predictive probabilities of the test set
    if (dm == 1){
      if (mix == "motg")
        test_pred[[j]]  <- post_pred_motgt(res = rslt[[j]]$res,
                                           test, K, pts, upper, lower,
                                           data_type, pr1, pr2)
      if (mix == "tmogt")
        test_pred[[j]]  <- post_pred_tmogt(res = rslt[[j]]$res,
                                           test, K, pts, upper, lower,
                                           data_type, pr1, pr2)

      prod_pred_ts[j]    <- sum(log(test_pred[[j]]$pred_ts))
      prod_pred_pt[j]    <- sum(log(test_pred[[j]]$pred_pt))

    } else {
      if (mix == "motg")
        test_pred[[j]]  <- post_pred_motgt_mv(res = rslt[[j]]$res,
                                              test, K, lower, upper, constr)

      if (mix == "tmogt")
        test_pred[[j]]  <- post_pred_tmogt_mv(res = rslt[[j]]$res, test, K,
                                              constr, thin_sc)

      prod_pred_ts[j]    <- sum(log(test_pred[[j]]$pred_ts))
    }

    count_bins[[j]] <- count_rej(res_list = rslt[[j]]$res)

    print(paste("Post_pred:", prod_pred_ts[j]))

    print("*******************************")
  }

  ppred_all <- prod_pred_ts
  time_all  <- time.taken

  return(list(result = rslt, ppred = ppred_all, time = time_all))
}

#' Rejection summary
#'
#' This function summarizes the value for the rejected samples
#' @param rejs vector of rejected samples
#' @param brks break points of the bins to summarize the rejections
#' @return number of samples in the bins
#' @export
#'
rej_summary <- function(rejs, brks){
  rej_all <- rejs
  bins    <- table(cut(rej_all, brks, include.lowest = T, right = F))
  return(bins)
}

#' Counting the rejections
#'
#' This function summarizes the value for the rejected samples
#' @param res_list list of the results of the gibbs sampling procedure
#' @return number of total rejections
#' @export
#'
count_rej <- function(res_list){
  unlist(lapply(1:length(res_list), function(x){
    sum(res_list[[x]]$cnt_rej)
  }))
}


