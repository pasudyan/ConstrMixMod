#' Update components parameter
#'
#' This function update the parameters of each cluster for the MoTg model
#' @param X data matrix where each column is independent observations and each row is the different dimension
#' @param mu mean of each cluster
#' @param cv covariance matrix of each cluster
#' @param niw_p Normal-inverse-wishart prior
#' @param constr function for the constraint
#' @param thr maximum number of rejections per observation
#' @param dm dimension of data
#' @return list of the new mean, covariance matrix, number of observations in each cluster, number of rejections, and the rejected samples
#' @export
updt_params_component <- function(X, mu, cv, niw_p, constr, thr = Inf, dm) {

  if (dm == 1) {
    X    <- matrix(X, nrow = 1, ncol = length(X))
    nm   <- length(X)
  } else {
    X    <- matrix(X, nrow = dm)
    nm   <- ncol(X)
  }

  # If there are observations in X
  if(nm > 0) {

    # Draw rejected samples from a Gaussian distribution
    if (thr != 0){
      y  <- impute_gauss(mu, cv, nm, constr, thr, dm, ceil_rej)
      X  <- cbind(X, y$rejs)
    } else
      y  <- NULL

    nm2  <- ncol(X)
    mn   <- rowMeans(X)

    if (dm == 1){
      S  <- sum(X^2)
    } else {
      S  <- (X - mn) %*% t(X - mn)
    }

  } else {
    # If there are no observations in X (specific cluster)
    # Keep the number of observations
    nm2  <- nm

    ## Draw from a standard multivariate Gaussian
    # Update mean to 0
    mn   <- rep(0, dm)

    if (dm == 1){
      S   <- 0
    } else {
      # Update S matrix to 0
      S   <- diag(rep(0, dm))
    }

    y   <- NULL
  }

  SS    <- list(n = nm2, mn = mn, S = S)
  # print(paste("SS:", SS))

  # Draw from the posterior
  if (dm == 1){
    rslt  <- nig.post(niw_p$m0, niw_p$v0, niw_p$a0, niw_p$b0, SS)
  } else {
    rslt  <- niw.post(niw_p$mu0, niw_p$lam, niw_p$Phi, niw_p$nu, SS)
  }

  return(list(mu = rslt$mu, cv = rslt$cv, no = nm,
              cnt_rej = nm2-nm, rej_samp = y))
}

#' Update components weight
#'
#' This function update the weights associated with each cluster
#' @param z cluster assignment for every observations
#' @param K total number of clusters
#' @param alp hyperparameter of the dirichlet process prior
#' @return vector of weights for every components
#' @export
updt_wghts <- function(z, K, alp){
  count <- rep(0, K)

  # Count the number of observations in each component
  for (i in 1:K){
    count[i] <- sum(z == i)
  }

  # Draw weights from a beta distribution
  cum    <- rev(cumsum(rev(count)))
  bt     <- rbeta(K - 1, 1 + count, alp + cum[2:K])
  bt[K]  <- 1
  cmp_bt <- 1 - bt
  wt     <- bt

  # Compute the weights of the stick breaking process
  if (K > 1){
    for (i in 2:K)
      wt[i] <- prod(cmp_bt[1:(i-1)])*bt[i]
  }

  return(wt)
}

#' Update assignment for univariate MoTG model
#'
#' This function updates the cluster assignment of each observation in the MoTG model
#' @param K max number of clusters
#' @param ip_data data matrix where columns are individual observations and rows are number of dimensions
#' @param N number of observations
#' @param dm number of dimensions
#' @param params list of parameters which include weights (wt), mean (mu), and covariance matrix (cv) for all clusters
#' @param lower lower bound of constraint (for unit square only)
#' @param upper upper bound of constraint (for unit square only)
#' @return vector of new cluster assignments
#' @export
updt_assgn_motgt <- function(K, ip_data, N, dm, params, lower, upper) {

  # Store new parameters and weights
  wt <- params$wt
  mu <- params$mu
  cv <- params$cv

  # Matrix to store likelihood, row is cluster, col is obs
  lik <- matrix(rep(0, N*K), nrow = K)

  # Compute likelihood for each observations belonging in
  # each component. Use of the likelihood function
  for(i in 1:K) {
    if (dm == 1){
      # Truncated Gaussian density
      lik[i,] <- truncnorm::dtruncnorm(ip_data, lower, upper, mean = mu[[i]],
                                sd = sqrt(cv[[i]][[1]]))
      lik[i,] <- edit_prob(lik[i,])
    } else {
      # Truncated Multivariate Gaussian density
      lik[i,] <- dtmvnorm_fast(t(ip_data), mean = as.vector(mu[[i]]),
                               sigma = cv[[i]],
                               lower = lower,
                               upper = upper,
                               log = T)
      lik[i,] <- edit_log_prob(lik[i,])
    }
  }

  if (dm == 1){
    if (K > 1){
      log_lik  <- apply(log(lik), 2, edit_log_prob)
      lik_w    <- log_lik + log(wt)
      lik_norm <- apply(lik_w, 2, prob_norm)
    } else {
      log_lik  <- edit_log_prob(log(lik))
      lik_w    <- log_lik + log(wt)
      lik_norm <- prob_norm(lik_w)
    }

    z  <- sapply(1:N, function(x){
      sample(K, size = 1, replace = T, prob = exp(lik_norm[, x]))
    })

  } else {
    lik <- lik + log(wt)
    lik <- lik - apply(lik, 2, max)
    lik <- t(t(lik) - log(colSums(exp(lik))))
    lik <- exp(lik)
    lik <- apply(lik, 2, cumsum)
    z   <- rowSums(runif(N) > t(lik)) + 1L
  }

  return(z)

}

is.empty <- function(list_of){
  unlist(lapply(1:length(list_of), function(x){
    if (length(list_of[[x]]) == 0)
      return(TRUE)
    else
      return(FALSE)
  }))
}

edit_prob_lognan <- function(X){
  X[is.nan(X) == T] <- 1e-05
  return(X)
}

# # ------------- Funtion to update Cluster assignments -------------- #
#' Update assignment for univariate MoTG model
#'
#' This function updates the cluster assignment of each observation in the MoTG model
#' @param K max number of clusters
#' @param ip_data data matrix where columns are individual observations and rows are number of dimensions
#' @param z vector of cluster assignment
#' @param dm number of dimensions
#' @param rej list of rejected proposals
#' @param N total number of observations
#' @param params list of parameters which include weights (wt), mean (mu), and covariance matrix (cv) for all clusters
#' @param thr maximum number of rejected samples per observations
#' @return vector of new cluster assignments
#' @export
updt_assgn_motg <- function(K, ip_data, z, dm, rej, N, params, thr) {

  # Store new parameters and weights
  wt <- params$wt
  mu <- params$mu
  cv <- params$cv

  # Matrix to store likelihood, row is cluster, col is obs
  lik <- matrix(rep(0, N*K), nrow = K)

  # Compute likelihood for each observations belonging in
  # each component. Use of the likelihood function
  for(i in 1:K) {

    # First compute the likelihood of every observations in cluster i
    lik[i,] <- mvnfast::dmvn(t(ip_data), mu = as.vector(mu[[i]]),
                             sigma = cv[[i]], log = T)

    if (thr != 0){
      # Find clusters with nonempty rejections
      # vector gives number of rejections in each cluster
      non_empty_rej <- unlist(lapply(1:length(rej), function(x){
        ll <- ncol(rej[[x]]$rejs)
        if (is.null(ll) == T)
          return(0)
        else
          return(ll)
      }))

      non_empty_cl <- which(non_empty_rej != 0)

      for (j in non_empty_cl){
        obs_in <- which(z == j)

        rejs_sub <- rej[[j]]

        # Find observations with rejections
        rej_obs <- unique(rejs_sub$ind_rej_obs)

        for (d in rej_obs){
          rejs_num  <- which(rejs_sub$ind_rej_obs == d)
          rejs_obs  <- rejs_sub$rejs[, rejs_num]
          if (dm == 1){
            rejs_prob <- sum(mvnfast::dmvn(matrix(rejs_obs, ncol = dm), mu = as.vector(mu[[i]]),
                                           sigma = cv[[i]], log = T))
          } else {
            rejs_prob <- sum(mvnfast::dmvn(t(rejs_obs), mu = as.vector(mu[[i]]),
                                           sigma = cv[[i]], log = T))
          }
          lik[i, obs_in[d]] <- lik[i, obs_in[d]] + rejs_prob
        }
      }

    }
  }

  if (dm == 1){
    if (K > 1){
      lik_w    <- lik + log(wt)
      lik_norm <- apply(lik_w, 2, prob_norm)
    } else {
      lik_w    <- lik + log(wt)
      lik_norm <- prob_norm(lik_w)
    }

    z  <- sapply(1:N, function(x){
      sample(K, size = 1, replace = T, prob = exp(lik_norm[, x]))
    })

  } else {
    lik <- lik + log(wt)
    lik <- lik - apply(lik, 2, max)
    lik <- t(t(lik) - log(colSums(exp(lik))))
    lik <- exp(lik)
    lik <- apply(lik, 2, cumsum)
    z   <- rowSums(runif(N) > t(lik)) + 1L
  }

  return(z)
}

#' Update assignment for univariate TMoG model
#'
#' This function updates the cluster assignment of each observation in the univariate TMoG model
#' @param K max number of clusters
#' @param ip_data data matrix where columns are individual observations and rows are number of dimensions
#' @param N total number of observations
#' @param dm number of dimensions
#' @param params list of parameters which include weights (wt), mean (mu), and covariance matrix (cv) for all clusters
#' @return vector of new cluster assignments
#' @export
updt_assgn_tmog <- function(K, ip_data, N, dm, params) {

  # Store new parameters and weights
  wt <- params$wt
  mu <- params$mu
  cv <- params$cv

  # Matrix to store likelihood, row is cluster, col is obs
  lik <- matrix(rep(0, N*K), nrow = K)

  # Compute likelihood for each observations belonging in
  # each component. Use of the likelihood function
  for(i in 1:K) {
    # Gaussian density
    if (dm == 1){
      lik[i,] <- mvnfast::dmvn(t(ip_data), mu = mu[[i]],
                                sigma = cv[[i]][[1]], log = T)

      # lik[i,] <- edit_prob_lognan(lik[i,])
      # lik[i,] <- edit_log_prob(lik[i,])

    } else {
      lik[i, ]  <- mvnfast::dmvn(t(ip_data), mu[[i]], sigma = cv[[i]], log = T)

    }
  }

  # browser()
  if (dm == 1){
    lik_w    <- lik + log(wt)
      if (K > 1){
        lik_norm <- apply(lik_w, 2, prob_norm)
      } else {
        lik_norm <- prob_norm(lik_w)
      }
    z  <- sapply(1:N, function(x){
      sample(K, size = 1, replace = T, prob = exp(lik_norm[,x]))
    })

  } else {
    # Normalized then sample
    lik <- lik + log(wt)
    lik <- lik - apply(lik, 2, max)
    lik <- t(t(lik) - log(colSums(exp(lik))))
    lik <- exp(lik)
    lik <- apply(lik, 2, cumsum)
    z   <- rowSums(runif(N) > t(lik)) + 1L
  }

  return(z)

}



updt_params <- function(K, alp, ip_data, dm, niw_p, z) {
  # This function update the parameters for each component using
  # the stick breaking construction after observations are placed
  # in a cluster

  mu  <- list()
  cv  <- list()
  no  <- rep(0, K)
  cnt_rej <- rep(0, K)

  count   <- rep(0, K)

  # Count number of observations in each cluster
  for(i in 1:K)
    count[i] <- sum(z == i)

  # Draw weights from stick breaking construction
  cum    <- rev(cumsum(rev(count)))
  bt     <- rbeta(K - 1, 1 + count, alp + c(cum[2:K]))
  bt[K]  <- 1
  cmp_bt <- 1 - bt
  wt     <- bt

  if (K > 1){
    for(i in 2:K)
      wt[i] <- prod(cmp_bt[1:(i-1)])*bt[i]
  }

  # Loop update params for each component
  for(i in 1:K) {
    nm   <- count[i]
    X    <- ip_data[, z == i, drop = F]

    # Update likelihood mean and S for each component
    if(nm > 0) {
      mn     <- rowMeans(X)

      if (dm == 1){
        S      <- sum(X^2) #sum of squares
      } else {
        S      <- (X - mn) %*% t(X - mn)
      }

    } else {
      mn     <- rep(0, dm)

      if (dm == 1){
        S      <- 0
      } else {
        S      <- diag(rep(0, dm))
      }
    }

    SS   <- list(n = nm, mn = mn, S = S)

    # Update posterior for each component
    if (dm == 1){
      rslt <- nig.post(niw_p$m0, niw_p$v0, niw_p$a0, niw_p$b0, SS)
    } else {
      rslt <- niw.post(niw_p$mu0, niw_p$lam, niw_p$Phi, niw_p$nu, SS)
    }

    mu[[i]] <- rslt$mu
    cv[[i]] <- rslt$cv
  }

  return(list(wt = wt, mu = mu, cv = cv, cnt_rej = cnt_rej, no = no))
}

#' Update parameters for TMoG model
#'
#' This function updates the cluster parameters in the univariate TMoG model
#' @param K max number of clusters
#' @param alp parameter for dirichlet process prior
#' @param x_rejs data matrix of rejected proposals
#' @param dm number of dimensions
#' @param wt weights of each cluster
#' @param niw_p normal-inverse-wishart prior
#' @param z_rejs cluster assignments for rejected proposals
#' @return list of new weights, mean, covariance matrix, number of rejected proposals, and number of observations in each cluster
#' @export
updt_params_tmog <- function(K, alp, x_rejs, dm, wt, niw_p, z_rejs) {
  # This function update the parameters for each component using
  # the stick breaking construction after observations are placed
  # in a cluster

  mu  <- list()
  cv  <- list()
  no  <- rep(0, K)
  count   <- rep(0, K)

  # Count number of observations in each cluster
  for(i in 1:K)
    count[i] <- sum(z_rejs == i)

  # Loop update params for each component
  for(i in 1:K) {
    nm   <- count[i]
    X    <- matrix(x_rejs[, z_rejs == i], nrow = dm)

    # Update likelihood mean and S for each component
    if(nm > 0) {
      mn     <- rowMeans(X)
      if (dm == 1){
        S  <- sum(X^2) #sum of squares
        if (S == Inf)
          S  <- 1e05
      } else {
        S <- (X - mn) %*% t(X - mn)
      }

    } else {
      mn  <- rep(0, dm)
      if (dm == 1){
        S  <- 0
      } else {
        S  <- diag(rep(0, dm))
      }
    }

    SS   <- list(n = nm, mn = mn, S = S)

    # Update posterior for each component
    if (dm == 1){
      rslt <- nig.post(m0 = niw_p$m0, v0 = niw_p$v0,
                       a0 = niw_p$a0,
                       b0 = niw_p$b0, SS)
    } else {
      rslt <- niw.post(niw_p$mu0, niw_p$lam, niw_p$Phi, niw_p$nu, SS)
    }

    mu[[i]] <- rslt$mu
    cv[[i]] <- rslt$cv
  }

  return(list(wt = wt, mu = mu, cv = cv, no = count))
}

