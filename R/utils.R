#' Rejection summary
#'
#' This function summarizes the value for the rejected samples
#' @param rejs vector of rejected samples
#' @param brks break points of the bins to summarize the rejections
#' @return number of samples in the bins
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
#'
count_rej <- function(res_list){
  unlist(lapply(1:length(res_list), function(x){
    sum(res_list[[x]]$cnt_rej)
  }))
}

#'  Cumulative distribution function
#'
#' This function evaluates the cumulative distribution of a univariate Gaussian distribution
#' @param mean mean of the distribution
#' @param sd standard deviation of the distribution
#' @param lower lower bound
#' @param upper upper bound
#' @return the area between the two bounds
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

#' Edit log probability
#'
#' This function edits the log probability for underflow and overflow. Substitute Inf to -Inf and NaN to -Inf.
#' @param X vector of probabilities.
#' @return vector of substituted log probabilities
edit_log_prob <- function(X){
  X[X == -Inf] <- -Inf
  X[X == Inf]  <- -Inf
  X[is.nan(X) == T]  <- -Inf
  return(X)
}

#' Mixture of Gaussians
#'
#' This function generates mixtures of gaussians samples
#' @param K number of clusters (default 2)
#' @param param_list list containing weight (wt), mean (mu), and covariance (cv)
#' @param N number of samples
#' @param dm dimensions of Gaussian
#' @return list with Gaussian samples and component assignment
gen_mog <- function(K = 2, param_list, N = 100, dm){
  wt  <- param_list$wt
  mu  <- param_list$mu
  cv  <- param_list$cv
  X   <- matrix(rep(0, N*dm), ncol = dm)
  z   <- sample(K, N, replace = T, prob = wt)
  for (i in 1:K) {
    indx  <- z == i
    if (sum(indx) > 0) {
      if (dm == 1){
        X[indx, ] <- mvnfast::rmvn(sum(indx), mu[[i]], cv[[i]])
      } else {
        X[indx, ] <- mvnfast::rmvn(sum(indx), mu[[i]], cv[[i]])
      }
    }
  }
  return(list(x = t(X), z = z))
}


#' Edit probability
#'
#' This function fixes the underflow and overflow of probability output with 0's
#' @param X vector of probabilities
#' @return vector of fixed probabilities
edit_prob <- function(X){
  X[is.nan(X) == T] <- 0
  X[X == Inf]       <- 0
  return(X)
}

#' Normalize probabilities
#'
#' This function normalizes a probability distribution
#' @param X vector of probabilities
#' @return normalize vector of probabilities
prob_norm <- function(X){
  norm <- X - logSumExp(X)
  return(norm)
}
