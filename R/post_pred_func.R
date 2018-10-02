#'  Cumulative distribution function
#'
#' This function evaluates the cumulative distribution of a univariate Gaussian distribution
#' @param mean mean of the distribution
#' @param sd standard deviation of the distribution
#' @param lower lower bound
#' @param upper upper bound
#' @return the area between the two bounds
#' @export
#' @examples
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' mc_inter(constr, 1/2, 0.1, 100)
#'
norm_cdf <- function(mean, sd, lower, upper){
  return(pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
}

#' Edit underflow of predictive probability
#'
#' This function edits the underflow of the predictive probability function for truncated Gaussian
#' @param X vector of probabilities
#' @export
edit_pred_prob <- function(X){
  # ================================= #
  # Function to edit overflows        #
  # Returns edited values             #
  # ================================= #

  ind_pinf <- which(X %in% c(Inf))
  ind_nan  <- which(is.nan(X) == T)
  if (length(ind_pinf) > 0)
    X[ind_pinf] <- 0
  if (length(ind_nan) > 0)
    X[ind_nan] <- 0
  return(X)
}

#' Posterior predictive for univariate MoTG
#'
#' This function evaluates the posterior predictive distribution for a univariate MoTG model
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test test set where columns are the individual observations and rows are the dimension
#' @param K number of clusters
#' @param pts grid points to draw the density distribution
#' @param upper upper bound of constraints, only use for unit square
#' @param lower lower bound of constraints, only use for unit square
#' @param data_type one of "gauss" and "beta"
#' @param pr1 first parameter for the density
#' @param pr2 second parameter for the density
#' @return list of prediction at points, prediction of tests, along with error bars
#' @export
#' @examples
#'
post_pred_motgt <- function(res, test, K, pts, upper, lower, data_type, pr1, pr2){
  iter = length(res)

  if (data_type == "beta"){
    p1_beta  <- pr1
    p2_beta  <- pr2
  }
  if (data_type == "gauss"){
    mu_r     <- pr1
    sd_r     <- pr2
  }

  grd_wt_pt   <- rep(0, length(pts))
  grd_wt_ts   <- rep(0, length(test))

  pred_pt   <- rep(list(), K)
  pred_ts   <- rep(list(), K)

  pred_pt_str  <- matrix(0, nrow = iter, ncol = length(pts))
  abs_err_pred <- matrix(0, nrow = iter, ncol = length(pts))
  sq_err_pred  <- matrix(0, nrow = iter, ncol = length(pts))

  probs_quant  <- c(0.10, 0.25, 0.5, 0.75, 0.90)

  pred_pt_summary <- matrix(0, nrow = length(probs_quant), ncol = length(pts))
  abs_err_summary <- matrix(0, nrow = length(probs_quant), ncol = length(pts))
  sq_err_summary  <- matrix(0, nrow = length(probs_quant), ncol = length(pts))

  if (data_type == "beta")
    real_dens <- dbeta(pts, p1_beta, p2_beta)
  if (data_type == "gauss")
    real_dens <- dnorm(pts, mu_r, sd_r)

  for (i in 1:iter){
    if(i %% 500 == 0) print(i)

    for(k in 1:K) {
      cv <- res[[i]]$cv[[k]][[1]]
      mn <- res[[i]]$mu[[k]]

      pred_pt[[k]] <- dtruncnorm(pts, a = lower, b = upper,
                                 mean = mn, sd = sqrt(cv))
      pred_ts[[k]] <- dtruncnorm(test, a = lower, b = upper,
                                 mean = mn, sd = sqrt(cv))

      pred_pt[[k]] <- edit_pred_prob(pred_pt[[k]])
      pred_ts[[k]] <- edit_pred_prob(pred_ts[[k]])

      grd_wt_pt   <- grd_wt_pt + res[[i]]$wt[k] * pred_pt[[k]]
      grd_wt_ts   <- grd_wt_ts + res[[i]]$wt[k] * pred_ts[[k]]

      pred_pt_str[i, ] <- pred_pt_str[i, ] +
        res[[i]]$wt[k] * pred_pt[[k]]
    }

    abs_err_pred[i, ] <- abs(pred_pt_str[i,] - real_dens)
    sq_err_pred[i, ]  <- (pred_pt_str[i,] - real_dens)^2

  }

  # Summary statistics of prediction point and errors
  pred_pt_summary <- apply(pred_pt_str, 2, function(x)
    quantile(x, probs = probs_quant))
  pred_pt_sd <- apply(pred_pt_str, 2, sd)

  abs_err_summary <- apply(abs_err_pred, 2, function(x)
    quantile(x, probs = probs_quant))
  sq_err_summary  <- apply(sq_err_pred, 2, function(x)
    quantile(x, probs = probs_quant))

  abs_err_sd <- apply(abs_err_pred, 2, sd)
  sq_err_sd  <- apply(sq_err_pred, 2, sd)

  grd_wt    <- list(pred_pt = grd_wt_pt/iter,
                    pred_ts = grd_wt_ts/iter,
                    pred_pt_summary = pred_pt_summary,
                    pred_pt_sd = pred_pt_sd,
                    abs_err_summary = abs_err_summary,
                    sq_err_summary  = sq_err_summary,
                    abs_err_sd = abs_err_sd,
                    sq_err_sd  = sq_err_sd)

  return(grd_wt)

}

#' Posterior predictive for univariate TMoG
#'
#' This function evaluates the posterior predictive distribution for a univariate TMoG model
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test test set where columns are the individual observations and rows are the dimension
#' @param K number of clusters
#' @param pts grid points to draw the density distribution
#' @param upper upper bound of constraints, only use for unit square
#' @param lower lower bound of constraints, only use for unit square
#' @param data_type one of "gauss" and "beta"
#' @param pr1 first parameter for the density
#' @param pr2 second parameter for the density
#' @return list of prediction at points, prediction of tests, along with error bars
#' @export
post_pred_tmogt <- function(res, test, K, pts, upper, lower, data_type, pr1, pr2){
  iter <- length(res)

  if (data_type == "beta"){
    p1_beta  <- pr1
    p2_beta  <- pr2
  }
  if (data_type == "gauss"){
    mu_r     <- pr1
    sd_r     <- pr2
  }

  grd_wt_pt   <- rep(0, length(pts))
  grd_wt_ts   <- rep(0, length(test))

  log_pred_pt_str <- matrix(0, ncol = length(pts), nrow = iter)
  log_pred_ts_str <- matrix(0, ncol = length(test), nrow = iter)

  pred_pt_str  <- matrix(0, ncol = length(pts), nrow = iter)
  pred_ts_str  <- matrix(0, ncol = length(test), nrow = iter)

  denom_pt    <- 0
  denom_ts    <- 0

  pred_pt_str  <- matrix(0, nrow = (iter), ncol = length(pts))
  abs_err_pred <- matrix(0, nrow = (iter), ncol = length(pts))
  sq_err_pred  <- matrix(0, nrow = (iter), ncol = length(pts))

  probs_quant  <- c(0.10, 0.25, 0.5, 0.75, 0.90)

  pred_pt_summary <- matrix(0, nrow = length(probs_quant), ncol = length(pts))
  abs_err_summary <- matrix(0, nrow = length(probs_quant), ncol = length(pts))
  sq_err_summary  <- matrix(0, nrow = length(probs_quant), ncol = length(pts))

  if (data_type == "beta")
    real_dens <- dbeta(pts, p1_beta, p2_beta)
  if (data_type == "gauss")
    real_dens <- dnorm(pts, mu_r, sd_r)

  for (i in 1:iter){
    if(i %% 500 == 0)
      print(i)

    for(k in 1:K) {
      cv <- res[[i]]$cv[[k]][[1]]
      mn <- res[[i]]$mu[[k]]

      log_pred_pt <- dnorm(pts, mean = mn, sd = sqrt(cv), log = T)
      log_pred_ts <- dnorm(test, mean = mn, sd = sqrt(cv), log = T)

      grd_wt_pt   <- grd_wt_pt + res[[i]]$wt[k] * exp(log_pred_pt)
      grd_wt_ts   <- grd_wt_ts + res[[i]]$wt[k] * exp(log_pred_ts)

      denom_pt    <- denom_pt + res[[i]]$wt[k] * norm_cdf(mn, sqrt(cv), lower, upper)
      denom_ts    <- denom_ts + res[[i]]$wt[k] * norm_cdf(mn, sqrt(cv), lower, upper)
    }

    if (log(denom_ts) == -Inf){
      log_pred_ts_str[i,] <- 0
    } else {
      log_pred_ts_str[i,]  <- log(grd_wt_ts) - log(denom_ts)
    }

    if (log(denom_pt) == -Inf){
      log_pred_pt_str[i,] <- 0
    } else {
      log_pred_pt_str[i,]  <- log(grd_wt_pt) - log(denom_pt)
    }

    abs_err_pred[i, ] <- abs(exp(log_pred_pt_str[i,]) - real_dens)
    sq_err_pred[i, ]  <- (exp(log_pred_pt_str[i,]) - real_dens)^2

    grd_wt_pt  <- rep(0, length(pts))
    grd_wt_ts  <- rep(0, length(test))

    denom_pt   <- 0
    denom_ts   <- 0

  }

  max_log_pred_pt <- apply(log_pred_pt_str, 2, max)
  max_log_pred_ts <- apply(log_pred_ts_str, 2, max)

  log_pred_pt     <- t(t(log_pred_pt_str) - max_log_pred_pt)
  log_pred_ts     <- t(t(log_pred_ts_str) - max_log_pred_ts)

  pred_pt  <- exp(log_pred_pt)
  pred_ts  <- exp(log_pred_ts)

  log_mn_pred_pt <- log(apply(pred_pt, 2, sum)) + max_log_pred_pt +
    log(1/(iter))
  log_mn_pred_ts <- log(apply(pred_ts, 2, sum)) + max_log_pred_ts +
    log(1/(iter))

  # Summary statistics of prediction point and errors
  pred_pt_summary <- apply(exp(log_pred_pt_str), 2, function(x)
    quantile(x, probs = probs_quant))
  pred_pt_sd <- apply(exp(log_pred_pt_str), 2, sd)

  abs_err_summary <- apply(abs_err_pred, 2, function(x)
    quantile(x, probs = probs_quant))
  sq_err_summary  <- apply(sq_err_pred, 2, function(x)
    quantile(x, probs = probs_quant))

  abs_err_sd <- apply(abs_err_pred, 2, sd)
  sq_err_sd  <- apply(sq_err_pred, 2, sd)

  grd_wt    <- list(pred_pt = exp(log_mn_pred_pt),
                    pred_ts = exp(log_mn_pred_ts),
                    pred_pt_summary = pred_pt_summary,
                    pred_pt_sd = pred_pt_sd,
                    abs_err_summary = abs_err_summary,
                    sq_err_summary  = sq_err_summary,
                    abs_err_sd = abs_err_sd,
                    sq_err_sd  = sq_err_sd)

  return(grd_wt)

}

#' Edit log probability
#'
#' This function edits the log probability for underflow and overflow. Substitute Inf to -Inf and NaN to -Inf.
#' @param X vector of probabilities.
#' @return vector of substituted log probabilities
#' @export
edit_log_prob <- function(X){
  X[X == -Inf] <- -Inf
  X[X == Inf]  <- -Inf
  X[is.nan(X) == T]  <- -Inf
  return(X)
}

#' Posterior predictive for multivariate MoTG
#'
#' This function evaluates the posterior predictive distribution for a multivariate MoTG model
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test test set where columns are the individual observations and rows are the dimension
#' @param K number of clusters
#' @param lower lower bound of constraints, only use for unit square
#' @param upper upper bound of constraints, only use for unit square
#' @param constr a function for the constraints
#' @return list of prediction of tests
#' @export
post_pred_motgt_mv <- function(res, test, K, lower, upper, constr){
  # ================================================================================ #
  # Posterior predictive distribution for Mixture of Truncated Gaussian Approximate  #
  # Only for Multi-dimensional data                                                  #
  # Returns posterior predictive of given points                                     #
  # ================================================================================ #

  iter <- length(res)
  dm   <- nrow(test)

  grd_wt_ts   <- rep(0, ncol(test))
  grd_wt   <- rep(0, ncol(test))

  pred_ts_str <- list(list())

  nonz_samp <- rep(0, iter)

  lower <- lower[1:dm]
  upper <- upper[1:dm]

  num_obs_prop <- 1e3

  for (i in 1:iter){
    if(i %% 100 == 0)
      print(i)

    pred_ts <- list()

    for(k in 1:K) {

      cv <- res[[i]]$cv[[k]][1:dm, 1:dm]
      mn <- res[[i]]$mu[[k]][1:dm]

      # Use of package tmvtnorm to compute Truncated Multivariate Gaussian density
      if (samp_method == "package")
        pred_ts[[k]] <- dtmvnorm_fast(t(test), mean = as.vector(mn), sigma = cv,
                                         lower = lower, upper = upper, log = T)
      if (samp_method == "is")
        pred_ts[[k]] <- dtmvnorm_impsamp_mv(t(test), mn = as.vector(mn), sig = cv,
                                            constr, num_obs_prop, log = T)

      pred_ts[[k]] <- edit_log_prob(pred_ts[[k]])

      grd_wt_ts   <- grd_wt_ts + res[[i]]$wt[k] * exp(pred_ts[[k]])

    }

    grd_wt <- grd_wt + grd_wt_ts
    grd_wt_ts <- rep(0, ncol(test))

  }

  return(list(pred_ts = grd_wt/iter))
}

#' Posterior predictive for multivariate TMoG
#'
#' This function evaluates the posterior predictive distribution for a multivariate TMoG model
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test test set where columns are the individual observations and rows are the dimension
#' @param K number of clusters
#' @param constr a function for the constraints
#' @param thin_sc thinning for MCMC chain (default no thinning)
#' @return list of prediction of tests
#' @export
post_pred_tmogt_mv <- function(res, test, K, constr, thin_sc = 1){
  # ================================================================================ #
  # Posterior predictive distribution for Truncated Mixtures of Gaussian             #
  # Only for Multi-dimensional data                                                  #
  # Returns posterior predictive of given points                                     #
  # ================================================================================ #

  iter <- length(res)/thin_sc

  dm   <- nrow(test)

  grd_wt_ts <- rep(0, ncol(test))
  denom_ts  <- 0

  log_pred_ts_str <- matrix(0, ncol = ncol(test), nrow = iter)
  nonz_samp <- rep(0, iter)

  if (samp_method == "package"){
    lower <- lower[1:dm]
    upper <- upper[1:dm]
  }

  num_obs_prop <- 1e3

  for (i in 1:iter){
    if(i %% 100 == 0)
      print(paste("i:", i))

    d <- i*thin_sc

    for(k in 1:K) {
      cv <- res[[d]]$cv[[k]][1:dm, 1:dm]
      mn <- res[[d]]$mu[[k]][1:dm]

      mu_prop  <- mn
      sig_prop <- cv*4

      # Compute Likelihood of test
      log_pred_ts <- mvnfast::dmvn(t(test), mu = mn, sigma = cv, log = TRUE)

      grd_wt_ts   <- grd_wt_ts + res[[d]]$wt[k] * exp(log_pred_ts)

      if (samp_method == "package")
        area_int  <- pmvnorm(lower, upper, mean = as.vector(mn), sigma = cv)[1]
      if (samp_method == "is")
        area_int  <- import_sampling_mv(constr, mn, cv, mu_prop, sig_prop, num_obs_prop)

      area_int    <- ifelse(area_int < 0, 0, area_int)
      denom_ts    <- denom_ts + res[[d]]$wt[k] * area_int

    }

    if (log(denom_ts) == -Inf){
      log_pred_ts_str[i, ] <- 0
    } else {
      log_pred_ts_str[i, ] <- log(grd_wt_ts) - log(denom_ts)
    }

    nonz_samp[i] <- ifelse(log(denom_ts) == -Inf, 0, 1)

    grd_wt_ts <- rep(0, ncol(test))
    denom_ts  <- 0
  }

  max_log_pred_ts <- apply(log_pred_ts_str, 2, max)
  log_pred_ts     <- t(t(log_pred_ts_str) - max_log_pred_ts)

  pred_ts  <- exp(log_pred_ts)
  eff_iter <- sum(nonz_samp)

  log_mn_pred_ts <- log(apply(pred_ts, 2, sum)) + max_log_pred_ts +
    log(1/eff_iter)

  return(list(pred_ts = exp(log_mn_pred_ts)))
}


