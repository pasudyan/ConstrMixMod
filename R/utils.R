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
  norm <- X - matrixStats::logSumExp(X)
  return(norm)
}

#' Index of rejections
#'
#' This function gives the rejected samples indices associated to the observation it belongs to
#' @param vl a vector of TRUE/FALSE indicating whether the samples are inside or outside the constraints
#' @param c_ind number of last accepted observations
give_ind <- function(vl, c_ind){
  indt   <- which(vl == TRUE)
  inda   <- rep(0, length(vl))
  inda[indt] <- Inf
  for (i in indt){
    inda[i:length(vl)] <- inda[i:length(vl)] + 1
  }

  return(list(ind_rej_obs = inda + c_ind, c_ind = length(indt)))
}

#' Index of last acceptance
#'
#' This function provides the index of the observations with that is last accepted
#' @param num_a total number of accepted samples
#' @param left number of observations without rejected samples
#' @param vl a vector of TRUE/FALSE indicating whether the samples are inside or outside the constraints
#'
ind_last_acc <- function(num_a, left, vl){
  # Function to determine index of last few rejections

  temp    <- num_a - (num_a - left)
  ind     <- which(vl == TRUE)[temp]
  t_ind   <- which(vl[1:ind] == FALSE)
  return(t_ind)
}

#' Imputation of rejected samples
#'
#' This function samples the rejected proposals for the Truncated Mixture of Gaussian (TMoG) model
#' @param K total number of clusters
#' @param params list of parameters which includes the weight (wt), mean (mu), and covariance matrix (cv)
#' @param N total number of observations
#' @param constr function for the constraint of the space
#' @param thr maximum number of rejections per observation
#' @param dm dimension of observations
#' @param ceil_rej max number of observations in total
#' @return rejected samples and indices to which observation it belongs to
#'
impute_rej <- function(K, params, N, constr, thr, dm, ceil_rej) {
  # Number of cumulative acceptance
  cum_acc <- 0

  # Number of rejected samples
  rej_acc <- 0

  rejs   <- matrix(rep(0, 0), nrow = length(params$mu[[1]]))
  z      <- c()

  left   <- N
  c_ind  <- 1
  ind_rej_obs <- c()

  while (left > 0) {
    if (thr == 0){
      left <- 0
      return(list(rejs = rejs, z = z, rej_acc = rej_acc))
    }

    # Sample from a mixture of gaussian
    smpls  <- gen_mog(K, params, N, dm)
    vl     <- constr(smpls$x)

    # Save accepted indices
    tmp     <- give_ind(vl, c_ind)
    c_ind   <- c_ind + tmp$c_ind

    # Number of accepted samples on this iteration
    num_a  <- sum(vl)

    # Number of rejected samples on this iteration
    num_r  <- N - num_a

    if (thr == Inf){
      if (num_a <= left){
        cum_acc <- cum_acc + num_a
        rej_acc <- rej_acc + num_r

        rejs    <- cbind(rejs, matrix(smpls$x[,!vl], nrow = dm))
        z       <- c(z, smpls$z[!vl])
        left    <- N - cum_acc

        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs)

        if (left == 0)
          return(list(rejs = rejs, z = z,
                      ind_rej_obs = ind_rej_obs[(ind_rej_obs != Inf) & (ind_rej_obs <= N)],
                      cum_acc = cum_acc, c_ind = c_ind, rej_acc = rej_acc))

        if (rej_acc >= ceil_rej)
          return(list(rejs = rejs, z = z,
                      ind_rej_obs = ind_rej_obs[(ind_rej_obs != Inf) & (ind_rej_obs <= N)],
                      cum_acc = cum_acc, c_ind = c_ind, rej_acc = rej_acc))


      } else {

        t_ind   <- ind_last_acc(num_a, left, vl)

        rejs    <- cbind(rejs, matrix(smpls$x[, t_ind], nrow = dm))
        z       <- c(z, smpls$z[t_ind])
        rej_acc <- ncol(rejs)

        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs[t_ind])

        cum_acc  <- min(cum_acc, N)

        return(list(rejs = rejs, z = z, ind_rej_obs = ind_rej_obs[ind_rej_obs != Inf],
                    cum_acc = cum_acc, c_ind = c_ind, rej_acc = rej_acc))
      }

    } else {

      cum_acc <- cum_acc + num_a
      rej_acc <- rej_acc + num_r

      if (cum_acc >= N & rej_acc < thr*N){
        if (rej_acc == 0){

          return(list(rejs = rejs, z = z, ind_rej_obs = ind_rej_obs[ind_rej_obs != Inf],
                      rej_acc = rej_acc, cum_acc = min(cum_acc, N)))
        } else {

          t_ind   <- ind_last_acc(num_a, left, vl)

          rejs    <- cbind(rejs, matrix(smpls$x[, t_ind], nrow = dm))
          z       <- c(z, smpls$z[t_ind])
          rej_acc <- ncol(rejs)
          ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs[t_ind])

          left    <- 0
          cum_acc <- min(cum_acc, N)

          return(list(rejs = rejs, z = z, ind_rej_obs = ind_rej_obs[ind_rej_obs != Inf],
                      rej_acc = rej_acc, cum_acc = cum_acc))
        }

      } else if (cum_acc < N & rej_acc >= thr*N){

        len_y   <- floor(thr*N)
        rejs    <- cbind(rejs, matrix(smpls$x[,!vl], nrow = dm))
        rejs    <- matrix(rejs[, 1:len_y], nrow = dm)
        z       <- c(z, smpls$z[!vl])

        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs)
        ind_rej_obs <- ind_rej_obs[ind_rej_obs != Inf]

        left    <- 0
        cum_acc <- min(cum_acc, N)

        return(list(rejs = rejs, z = z[1:len_y], ind_rej_obs = ind_rej_obs[1:len_y],
                    rej_acc = len_y, cum_acc = cum_acc))

      } else if (cum_acc >= N & rej_acc >= thr*N){

        t_ind   <- ind_last_acc(num_a, left, vl)
        rejs    <- cbind(rejs, matrix(smpls$x[,t_ind], nrow = dm))
        rej_acc <- ncol(rejs)
        z       <- c(z, smpls$z[t_ind])

        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs[t_ind])
        ind_rej_obs <- ind_rej_obs[ind_rej_obs != Inf]

        len_y   <- min(rej_acc, floor(thr*N))
        rejs    <- matrix(rejs[, 1:len_y], nrow = dm)

        left    <- 0

        cum_acc <- min(cum_acc, N)

        return(list(rejs = rejs, z = z[1:len_y], ind_rej_obs = ind_rej_obs[1:len_y],
                    rej_acc = len_y, cum_acc = cum_acc))
      } else {
        left    <- N - cum_acc
        rejs    <- cbind(rejs, matrix(smpls$x[,!vl], nrow = dm))
        z       <- c(z, smpls$z[!vl])
        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs)
      }
    }
  }
}

#' Imputation of rejected samples
#'
#' This function samples the rejected proposals for the Mixture of Truncated Gaussians (MoTG) model
#' @param mu mean parameters
#' @param cv covariance matrix
#' @param N total number of observations
#' @param constr function for the constraint of the space
#' @param thr maximum number of rejections per observation
#' @param dm dimension of observations
#' @param ceil_rej max number of observations in total
#' @return rejected samples and indices to which observation it belongs to
#'
impute_gauss <- function(mu, cv, N, constr, thr, dm, ceil_rej) {

  # Number of cumulative acceptance
  cum_acc <- 0

  # Number of rejected samples
  rej_acc <- 0

  # Store index of acceptance
  ind_rej_obs <- c()
  rejs        <- matrix(rep(0, 0), nrow = length(mu))

  left  <- N
  c_ind <- 1

  while (left > 0) {

    if (dm == 1){
      # Sample from a univariate normal distribution
      smpls   <- matrix(rnorm(N, mu, sqrt(cv[[1]])), nrow = dm, ncol = N)
    } else {
      # Sample from a multivariate normal distribution with the component parameters
      smpls   <- t(mvrnorm(N, mu, cv))
    }

    # Check if samples satisfies constraint
    vl      <- constr(smpls)

    # Save accepted indices
    tmp     <- give_ind(vl, c_ind)
    c_ind   <- c_ind + tmp$c_ind

    # Number of accepted samples on this iteration
    num_a  <- sum(vl)

    # Number of rejected samples on this iteration
    num_r  <- N - num_a

    if (thr == Inf){

      if (num_a <= left){
        # If the number of acceptance is still less than how much acceptance
        # needed
        cum_acc <- cum_acc + num_a # total number of accepted samples
        rej_acc <- rej_acc + num_r # total number of rejected samples

        rejs        <- cbind(rejs, matrix(smpls[,!vl], nrow = dm))
        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs)
        left        <- N - cum_acc

        if (left == 0)
          return(list(rejs = rejs,
                      ind_rej_obs = ind_rej_obs[(ind_rej_obs != Inf) & (ind_rej_obs <= N)],
                      cum_acc = cum_acc, c_ind = c_ind))

        if (rej_acc >= ceil_rej){
          return(list(rejs = rejs,
                      ind_rej_obs = ind_rej_obs[(ind_rej_obs != Inf) & (ind_rej_obs <= N)],
                      cum_acc = cum_acc, c_ind = c_ind))
        }

      } else {
        # If the number of acceptance is greater than how much acceptance is
        # still needed

        t_ind   <- ind_last_acc(num_a, left, vl)
        cum_acc <- cum_acc + num_a

        rejs        <- cbind(rejs, matrix(smpls[, t_ind], nrow = dm))
        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs[t_ind])
        rej_acc     <- ncol(rejs)

        cum_acc <- min(cum_acc, N)

        return(list(rejs = rejs, ind_rej_obs = ind_rej_obs[ind_rej_obs != Inf],
                    cum_acc = cum_acc, c_ind = c_ind))
      }

    } else {
      # If threshold is not infinity

      cum_acc <- cum_acc + num_a
      rej_acc <- rej_acc + num_r

      if (cum_acc >= N & rej_acc < thr*N){
        # If acceptance is more than total number of observations,
        # but number of rejections are less than threshold

        if (rej_acc == 0){

          # Keep all rejections
          return(list(rejs = rejs, ind_rej_obs = ind_rej_obs[ind_rej_obs != Inf],
                      cum_acc = min(cum_acc, N)))

        } else {

          t_ind   <- ind_last_acc(num_a, left, vl)

          rejs        <- cbind(rejs, matrix(smpls[, t_ind], nrow = dm))
          ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs[t_ind])
          rej_acc     <- ncol(rejs)

          left    <- 0
          cum_acc <- min(cum_acc, N)

          # Keep all rejections
          return(list(rejs = rejs, ind_rej_obs = ind_rej_obs[ind_rej_obs != Inf],
                      cum_acc = cum_acc))
        }

      } else if (cum_acc < N & rej_acc >= thr*N){
        # If acceptance is less than total number of observations,
        # but number of rejections are more than threshold

        len_y       <- floor(thr*N)
        rejs        <- cbind(rejs, matrix(smpls[, !vl], nrow = dm))
        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs)
        ind_rej_obs <- ind_rej_obs[ind_rej_obs != Inf]

        left    <- 0
        cum_acc <- min(cum_acc, N)

        return(list(rejs = matrix(rejs[, 1:len_y], nrow = dm),
                    ind_rej_obs = ind_rej_obs[1:len_y],
                    cum_acc = cum_acc))

      } else if (cum_acc >= N & rej_acc >= thr*N){
        # If acceptance is more than total number of observations,
        # and number of rejections are also more than threshold

        t_ind   <- ind_last_acc(num_a, left, vl)

        rejs    <- cbind(rejs, matrix(smpls[, t_ind], nrow = dm))
        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs[t_ind])
        ind_rej_obs <- ind_rej_obs[ind_rej_obs != Inf]
        rej_acc <- ncol(rejs)

        len_y   <- min(rej_acc, floor(thr*N))
        left    <- 0

        cum_acc <- min(cum_acc, N)

        return(list(rejs = matrix(rejs[, 1:len_y], nrow = dm),
                    ind_rej_obs = ind_rej_obs[1:len_y],
                    cum_acc = cum_acc))

      } else {
        left    <- N - cum_acc
        rejs    <- cbind(rejs, matrix(smpls[,!vl], nrow = dm))
        ind_rej_obs <- c(ind_rej_obs, tmp$ind_rej_obs)
      }
    }
  }
}


#' Monte Carlo integration function
#'
#' This function evaluates an integral using Monte Carlo integration.
#' @param constr a function that describes the domain or constraint of the observations
#' @param mn mean of the observations to be generated
#' @param cv covariance matrix of the observations to be generated
#' @param num_obs number of observations to be generated
#' @return value of integrals
#' @examples
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' mc_inter(constr, 1/2, 0.1, 100)
mc_inter <- function(constr, mn, cv, num_obs){
  # Monte Carlo integration:
  rsamp    <- matrix(mvnfast::rmvn(num_obs, mn, cv), nrow = 1)
  rsamp_in <- rsamp[constr(rsamp)]

  mc_int   <- (1/num_obs)*length(rsamp_in)
  return(mc_int)
}

#' Important Sampling function
#'
#' This function evaluates integrals using an important sampling method with a Gaussian distribution proposal
#' @param constr a function that describes the domain or constraint of the observations
#' @param mn mean of the observations to be generated
#' @param cv covariance matrix of the observations to be generated
#' @param num_obs number of observations to be generated
#' @return value of integrals
#' @examples
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' important_sampling(constr, 1/2, 0.1, 100)
#'
import_sampling <- function(constr, mn, cv, num_obs){
  # Parameters of proposal - q
  mu_prop <- 0
  cv_prop <- 4

  # Sample from standard normal multivariate distribution - q
  rsamp   <- t(mvnfast::rmvn(num_obs, mu_prop, cv_prop))

  # Compute log likelihoods for weight
  rsamp_pdf  <- mvnfast::dmvn(t(rsamp), mu = mn, sigma = cv, log = T)
  rsamp_prop <- mvnfast::dmvn(t(rsamp), mu = mu_prop, sigma = cv_prop, log = T)

  # Compute log weights
  log_is_w <- rsamp_pdf - rsamp_prop

  rsamp_in <- constr(rsamp)

  imps <- (1/num_obs)*sum(exp(log_is_w[rsamp_in]))

  return(imps)
}

#' Normalized Important Sampling function
#'
#' This function evaluates integrals using a normalized important sampling method with a Gaussian distribution proposal
#' @param constr a function that describes the domain or constraint of the observations
#' @param mn mean of the observations to be generated
#' @param cv covariance matrix of the observations to be generated
#' @param num_obs number of observations to be generated
#' @return value of integrals
#' @examples
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' important_sampling_norm(constr, 1/2, 0.1, 100)
import_sampling_norm <- function(constr, mn, cv, num_obs){
  # Normalized important sampling
  # Parameters of proposal - q
  mu_prop <- 0
  cv_prop <- 4

  # Sample from standard normal multivariate distribution - q
  rsamp   <- t(mvnfast::rmvn(num_obs, mu_prop, cv_prop))

  # Compute log likelihoods for weight
  rsamp_pdf  <- mvnfast::dmvn(t(rsamp), mu = mn, sigma = cv, log = T)
  rsamp_prop <- mvnfast::dmvn(t(rsamp), mu = mu_prop, sigma = cv_prop, log = T)

  # Compute log weights
  log_is_w <- rsamp_pdf - rsamp_prop

  rsamp_in <- constr(rsamp)

  norm <- (1/num_obs)*sum(exp(log_is_w))
  imps <- ((1/num_obs)*sum(exp(log_is_w[rsamp_in])))/norm

  return(imps)
}

#' Truncated univariate Gaussian density
#'
#' This function evaluates the truncated univariate Gaussian density using an important sampling procedure as the normalizing constant
#' @param X a vector of observations where each column is the independent observations
#' @param mn mean of the Gaussian density
#' @param sig covariance matrix of the Gaussian density
#' @param constr a function that describes the domain or constraint of the samples
#' @param num_samp number of samples in the important sampling procedure
#' @return value of the density
#' @examples
#' X <- rnorm(5, mean=0, sd=1)
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' dtmvnorm_impsamp(X, mn, sig, constr, num_samp=1e6)
dtmvnorm_impsamp <- function(X, mn, sig, constr, num_samp = 1e6){
  numer <- dmvn(X, mu = mn, sigma = sig, log = T)
  denom <- log(import_sampling(constr, mn = mn, cv = sig, num_obs = num_samp))
  val   <- numer - denom

  val    <- ifelse(constr(X), val, -Inf)

  return(exp(val))
}

#' Important Sampling for Multivariate variables
#'
#' This function evaluates integrals of a multivariate random variables using important sampling method with a multivariate Gaussian proposal
#' @param constr a function that describes the domain or constraint of the observations
#' @param mn mean of the observations to be generated
#' @param cv covariance matrix of the observations to be generated
#' @param mu_prop mean of the proposal distribution
#' @param sig_prop covariance matrix of the proposal distribution
#' @param num_obs number of observations to be generated
#' @return value of integrals
#' @examples
#' genconst <- function(u1, u2, l1, l2) { in_set <- function(X) {
#'  if (nrow(X) == 1) X = t(X)
#'  X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 & X[3,] > l1 &
#'  X[3,] < u1 & X[4,] > l2 & X[4,] < u2 }}
#' constr  <- genconst(u = 0, u = 0, l1 = 1, l2 = 1)
#' import_sampling_mv(constr, rep(1/2, 4), diag(0.1,4), rep(0,4), diag(0.3,4), 100)
import_sampling_mv <- function(constr, mn, cv, mu_prop, sig_prop, num_obs){

  dm <- length(mn)

  # Parameters of proposal - q
  cv_prop <- sig_prop

  # Sample from standard normal multivariate distribution - q
  rsamp   <- t(mvnfast::rmvn(num_obs, mu_prop, cv_prop))

  # Compute log likelihoods for weight
  rsamp_pdf  <- mvnfast::dmvn(t(rsamp), mu = mn, sigma = cv, log = T)
  rsamp_prop <- mvnfast::dmvn(t(rsamp), mu = mu_prop, sigma = cv_prop, log = T)

  # Compute log weights
  log_is_w <- rsamp_pdf - rsamp_prop

  rsamp_in <- constr(rsamp)

  imps <- (1/num_obs)*sum(exp(log_is_w[rsamp_in]))

  return(imps)
}

#' Truncated Multivariate Gaussian density
#'
#' This function evaluates the truncated multivariate Gaussian density using an important sampling procedure as the normalizing constant
#' @param X a vector of observations where each row is the dimension and each column is the independent observations
#' @param mn mean of the multivariate Gaussian density
#' @param sig covariance matrix of the multivariate Gaussian density
#' @param constr a function that describes the domain or constraint of the samples
#' @param num_samp number of samples in the important sampling procedure
#' @param log whether log value is returned (default TRUE)
#' @return value of density
#' @examples
#' genconst <- function(u1, u2, l1, l2) { in_set <- function(X) {
#'  if (nrow(X) == 1) X = t(X)
#'  X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 & X[3,] > l1 &
#'  X[3,] < u1 & X[4,] > l2 & X[4,] < u2 }}
#' constr  <- genconst(u1 = 0, u2 = 0, l1 = 1, l2 = 1)
#' dtmvnorm_impsamp_mv(X = mvnfast::rmvn(5, mu = rep(0,4), sigma = diag(0.1, 4)), mn = rep(0, 4), sig = diag(0.1,4), constr, num_samp = 1000, log = F)
dtmvnorm_impsamp_mv <- function(X, mn, sig, constr, num_samp = 10000, log = T){
  numer <- dmvn(X, mu = mn, sigma = sig, log = T)

  mn_imp  <- mn
  sig_imp <- sig*4

  denom <- log(import_sampling_mv(constr, mn = mn, cv = sig,
                                  mu_prop = mn_imp, sig_prop = sig_imp,
                                  num_obs = num_samp))
  val   <- numer - denom

  val    <- ifelse(constr(t(X)), val, -Inf)

  if (log)
    return(val)
  else
    return(exp(val))
}

#' Gibbs sampling algorithm
#'
#' This computes the Gibbs sampling procedure
#' @param K maximum number of clusters
#' @param ip_data data matrix where columns are individual observations and rows are dimensions
#' @param thr maximum number of rejections per observations
#' @param burnIn number of burn-in for the chain. Must be less than number of iterations
#' @param numIter total number of iterations
#' @param constr function for the constraints
#' @param mix one of "tmog" and "motg"
#' @param lower lower constraint for unit square
#' @param upper upper constriant for unit square
#' @param data_type one of "gauss", "beta", "flowcyto", "crime", "bvgauss"
#' @return a list of results consisting of weights, mean, and covariance matrix for every iterations, along with the rejected proposals
#' @export
#'
gibbs_sampling <- function(K, ip_data, thr = Inf, burnIn, numIter, constr, mix, lower, upper, data_type) {

  # Dimension and number of observations
  dm      <- dim(ip_data)[1]
  N       <- dim(ip_data)[2]

  if (dm == 1)
    mu_p    <- mean(ip_data)

  # Number of iteration and alpha for dirichlet
  alp      <- 1

  # NIW Prior
  if (dm == 1){
    if (data_type == "gauss"){
      niw_p  <- list(m0 = mu_p, v0 = 2,
                     a0 = 2,  b0 = 5/100)
    } else if (data_type == "beta"){
      niw_p  <- list(m0 = mu_p, v0 = 2,
                     a0 = 2,  b0 = 1/100)
    } else {
      stop("Data type for one-dimension is not one of Gaussian or Beta")
    }

  } else {
    if (data_type == "flowcyto"){
      niw_p  <- list(mu0 = rep(1/2, dm), lam = 1e-2,
                     Phi = diag(rep(1e-3, dm)),  nu = dm+1)
    } else if (data_type == "bvgauss"){
      niw_p  <- list(mu0 = 1/2, lam = 1/10,
                     Phi = diag(rep(1e-3, dm)),  nu = dm+1)
    } else if (data_type == "crime"){
      niw_p  <- list(mu0 = 1/2, lam = 1/10,
                     Phi = diag(rep(1e-3, dm)),  nu = dm+1)
    } else {
      stop("Data type for multi-dimensions is not one of Flowcytometry, BVgauss, or Crime")
    }
  }

  # Initialize clusters and parameters for each cluster
  z       <- kmeans(t(ip_data), 5, nstart = 10, iter.max = 50)$cluster
  params  <- updt_params(K, alp, ip_data, dm, niw_p, z)

  # Number of Rejections per Iteration
  rej_samp <- rep(list(list()), K)
  params$rej_samp <- rej_samp

  # Store result
  rslt    <- rep(list(params), (numIter-burnIn))

  brks  <- c(seq(-4, 0, 0.05), seq(1, 4, 0.05))
  count_bins_sum <- rej_summary(-5, brks)

  # *********************** Gibbs Sampling loop ************************* #
  for(iter in 1:numIter) {
    if(iter%%100 == 0)
      print(paste("Iter", iter))

    # Reorder observations
    new_order <- sample(ncol(ip_data), replace = F)
    ip_data   <- matrix(ip_data[, new_order], nrow = dm)
    z <- z[new_order]

    if (mix == "motg" | mix == "motgt"){
      # Update parameters for each cluster by first sampling the rejections
      for(i in 1:K) {
        pm <- updt_params_component(X = ip_data[, z == i],
                                    mu = params$mu[[i]],
                                    cv = params$cv[[i]],
                                    niw_p, constr, thr, dm)

        # Store parameters
        params$mu[[i]] <- pm$mu
        params$cv[[i]] <- pm$cv

        # Count number of rejections in partition
        if (mix == "motgt"){
          if (dm == 1){
            if (iter == numIter)
              rej_samp[[i]] <- pm$rej_samp
          } else {

            # For now, save rejections
            params$rej_samp[[i]] <- pm$rej_samp
          }
        } else {
          params$rej_samp[[i]] <- pm$rej_samp
        }

        params$cnt_rej[i]  <- pm$cnt_rej
        params$no[i]       <- pm$no

      }

      # Update assignments
      if (mix == "motgt"){
        z  <- updt_assgn_motgt(K, ip_data, N, dm, params, lower, upper)
      } else {
        z  <- updt_assgn_motg(K, ip_data, z, dm, rej = params$rej_samp, N, params, thr)
      }

      # Store new cluster parameters
      params$z <- z

      # Update weights
      new_wt <- updt_wghts(z, K, alp)
      params$wt <- new_wt

      if (iter > burnIn)
        rslt[[iter-burnIn]] <- params

      # Set rejection storage to empty again
      rej_samp <- rep(list(list()), K)

    } else if (mix == "tmog"){
      # Sample rejections
      smpls_rejs <- impute_rej(K, params, N, constr, thr, dm, ceil_rej)

      # Function to update params
      params <- updt_params_tmog(K, alp,
                                 x_rejs = cbind(ip_data, smpls_rejs$rejs),
                                 dm, params$wt,
                                 niw_p, z_rejs = c(z, smpls_rejs$z))
      # Update weights
      new_wt <- updt_wghts(c(z, smpls_rejs$z), K, alp)

      params$wt <- new_wt

      # Count number of rejections in partition
      if (dm == 1){
        # only save rejs summary for one dimensional data
        if (thr != 0){
          count_bins <- rej_summary(smpls_rejs$rejs, brks)
          count_bins_sum <- count_bins_sum + count_bins
        } else {
          count_bins_sum <- NULL
        }
      } else {
        # For now, save rejections
        if (iter > burnIn)
          rslt[[iter-burnIn]]$rej_smpls <- smpls_rejs$rejs
      }

      # Update assignments
      z <- updt_assgn_tmog(K, ip_data, N, dm, params)

      if (iter > burnIn){
        rslt[[iter-burnIn]] <- params
        rslt[[iter-burnIn]]$cnt_rej  <- smpls_rejs$rej_acc
      }
    } else {
      stop("Mixture is not one of motg, motgt, or tmog")
    }

  }

  return(list(res = rslt, rej_samp = rej_samp))

}









