#' Index of rejections
#'
#' This function gives the rejected samples indices associated to the observation it belongs to
#' @param vl a vector of TRUE/FALSE indicating whether the samples are inside or outside the constraints
#' @param c_ind number of last accepted observations
#' @export
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
#' @export
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
#' @export
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
#' @export
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





