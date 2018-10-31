#' Update components parameter
#'
#' This function update the parameters of each cluster for the MoTg model
#' @param X data matrix where each column is independent observations and each row is the different dimension
#' @param mu mean of each cluster
#' @param cv covariance matrix of each cluster
#' @param niwpr Normal-inverse-wishart prior
#' @param constr function for the constraint
#' @param thr maximum number of rejections per observation
#' @return list of the new mean, covariance matrix, number of observations in each cluster, number of rejections, and the rejected samples
updateParamsComponent <- function(X, mu, cv, niwpr, constr, thr = Inf) {

  dm <- dim(X)[1]

  if(is.null(dm))
    dm <- 1

  if (dm == 1) {
    X  <- matrix(X, nrow = 1, ncol = length(X))
    nm <- length(X)
  } else {
    X  <- matrix(X, nrow = dm)
    nm <- ncol(X)
  }

  if(nm > 0) {
    if (thr != 0){
      y <- imputeGauss(mu, cv, nm, constr, thr, ceilRej)
      if (!is.null(y$indRejObs)){
        X <- cbind(X, y$rejs)
      }
    } else{
      y <- NULL
    }

    nm2 <- ncol(X)
    mn  <- rowMeans(X)

    if (dm == 1) {
      S  <- sum(X^2)
    } else {
      S  <- (X - mn) %*% t(X - mn)
    }
  } else {
    nm2 <- nm
    mn  <- rep(0, dm)
    if (dm == 1){
      S <- 0
    } else {
      S <- diag(rep(0, dm))
    }
    y <- NULL
  }

  SS <- list(n = nm2, mn = mn, S = S)

  if (dm == 1){
    rslt  <- nigPost(niwpr$m0, niwpr$v0, niwpr$a0, niwpr$b0, SS)
  } else {
    rslt  <- niwPost(niwpr$mu0, niwpr$lam, niwpr$Phi, niwpr$nu, SS)
  }

  return(list(mu = rslt$mu, cv = rslt$cv, no = nm,
              cntRej = nm2 - nm, rejSamp = y))
}

#' Update components weight
#'
#' This function update the weights associated with each cluster
#' @param z cluster assignment for every observations
#' @param K total number of clusters
#' @param alp hyperparameter of the dirichlet process prior
#' @return vector of weights for every components
updateWeights <- function(z, K, alp){
  count <- rep(0, K)

  # Count the number of observations in each component
  for (i in 1:K){
    count[i] <- sum(z == i)
  }

  # Draw weights from a beta distribution
  cum    <- rev(cumsum(rev(count)))
  bt     <- rbeta(K - 1, 1 + count, alp + cum[2:K])
  bt[K]  <- 1
  cmpbt  <- 1 - bt
  wt     <- bt

  # Compute the weights of the stick breaking process
  if (K > 1){
    for (i in 2:K)
      wt[i] <- prod(cmpbt[1:(i-1)])*bt[i]
  }

  return(wt)
}


# # ------------- Funtion to update Cluster assignments -------------- #
#' Update assignment for univariate MoTG model
#'
#' This function updates the cluster assignment of each observation in the MoTG model
#' @param K max number of clusters
#' @param ipdata data matrix where columns are individual observations and rows are number of dimensions
#' @param z vector of cluster assignment
#' @param rej list of rejected proposals
#' @param N total number of observations
#' @param params list of parameters which include weights (wt), mean (mu), and covariance matrix (cv) for all clusters
#' @param thr maximum number of rejected samples per observations
#' @return vector of new cluster assignments
updateClustAssignMotg <- function(K, ipdata, z, rej, N, params, thr) {

  dm <- dim(ipdata)[1]

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
    lik[i,] <- mvnfast::dmvn(t(ipdata), mu = as.vector(mu[[i]]),
                             sigma = cv[[i]], log = T)

    if (thr != 0){
      # Find clusters with nonempty rejections
      # vector gives number of rejections in each cluster
      nonEmptyRej <- unlist(lapply(1:length(rej), function(x){
        ll <- ncol(rej[[x]]$rejs)
        if (is.null(ll) == T)
          return(0)
        else
          return(ll)
      }))

      nonEmptyCl <- which(nonEmptyRej != 0)

      for (j in nonEmptyCl){
        obsIn <- which(z == j)

        rejSub <- rej[[j]]

        # Find observations with rejections
        rejObs <- unique(rejSub$indRejObs)

        for (d in rejObs){
          rejsNum  <- which(rejSub$indRejObs == d)
          rejSubObs  <- rejSub$rejs[, rejsNum]
          if (dm == 1){
            rejsProb <- sum(mvnfast::dmvn(matrix(rejSubObs, ncol = dm), mu = as.vector(mu[[i]]),
                                           sigma = cv[[i]], log = T))
          } else {
            rejsProb <- sum(mvnfast::dmvn(t(rejSubObs), mu = as.vector(mu[[i]]),
                                           sigma = cv[[i]], log = T))
          }
          lik[i, obsIn[d]] <- lik[i, obsIn[d]] + rejsProb
        }
      }

    }
  }

  likW    <- lik + log(wt)
  likNorm <- apply(likW, 2, probNorm)
  z  <- sapply(1:N, function(x){
    sample(K, size = 1, replace = T, prob = exp(likNorm[, x]))
  })

  return(z)
}

#' Update assignment for univariate TMoG model
#'
#' This function updates the cluster assignment of each observation in the univariate TMoG model
#' @param K max number of clusters
#' @param ipdata data matrix of size p-by-N where N are number of observations and p are number of dimensions
#' @param N total number of observations
#' @param params list of parameters which include weights (wt), mean (mu), and covariance matrix (cv) for all clusters
#' @return vector of new cluster assignments
updateClustAssignTmog <- function(K, ipdata, N, params) {

  dm <- dim(ipdata)[1]

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
      lik[i,] <- mvnfast::dmvn(t(ipdata), mu = mu[[i]],
                                sigma = cv[[i]][[1]], log = T)
    } else {
      lik[i,] <- mvnfast::dmvn(t(ipdata), mu[[i]], sigma = cv[[i]], log = T)
    }
  }

  # Check equality here
  likW    <- lik + log(wt)
  likNorm <- apply(likW, 2, probNorm)
  z  <- sapply(1:N, function(x){
    sample(K, size = 1, replace = T, prob = exp(likNorm[,x]))
  })

  return(z)

}

updateParams <- function(K, alp, ipdata, niwpr, z) {
  # This function update the parameters for each component using
  # the stick breaking construction after observations are placed
  # in a cluster

  dm  <- dim(ipdata)[1]
  mu  <- list()
  cv  <- list()
  count  <- rep(0, K)

  # Count number of observations in each cluster
  for(i in 1:K)
    count[i] <- sum(z == i)

  # Draw weights from stick breaking construction
  cum    <- rev(cumsum(rev(count)))
  bt     <- rbeta(K - 1, 1 + count, alp + c(cum[2:K]))
  bt[K]  <- 1
  cmpbt  <- 1 - bt
  wt     <- bt

  if (K > 1){
    for(i in 2:K)
      wt[i] <- prod(cmpbt[1:(i-1)])*bt[i]
  }

  # Loop update params for each component
  for(i in 1:K) {
    nm <- count[i]
    X  <- ipdata[, z == i, drop = F]

    # Update likelihood mean and S for each component
    if(nm > 0) {
      mn <- rowMeans(X)
      if (dm == 1){
        S <- sum(X^2) #sum of squares
      } else {
        S <- (X - mn) %*% t(X - mn)
      }
    } else {
      mn <- rep(0, dm)
      if (dm == 1){
        S <- 0
      } else {
        S <- diag(rep(0, dm))
      }
    }

    SS <- list(n = nm, mn = mn, S = S)

    # Update posterior for each component
    if (dm == 1){
      rslt <- nigPost(niwpr$m0, niwpr$v0, niwpr$a0, niwpr$b0, SS)
    } else {
      rslt <- niwPost(niwpr$mu0, niwpr$lam, niwpr$Phi, niwpr$nu, SS)
    }
    mu[[i]] <- rslt$mu
    cv[[i]] <- rslt$cv
  }
  return(list(wt = wt, mu = mu, cv = cv))
}

#' Update parameters for TMoG model
#'
#' This function updates the cluster parameters in the univariate TMoG model
#' @param K max number of clusters
#' @param alp parameter for dirichlet process prior
#' @param xRejs data matrix of rejected proposals
#' @param wt weights of each cluster
#' @param niwpr normal-inverse-wishart prior
#' @param zRejs cluster assignments for rejected proposals
#' @return list of new weights, mean, covariance matrix, number of rejected proposals, and number of observations in each cluster
updateParamsTmog <- function(K, alp, xRejs, wt, niwpr, zRejs) {
  # This function update the parameters for each component using
  # the stick breaking construction after observations are placed
  # in a cluster

  dm  <- dim(xRejs)[1]
  mu  <- list()
  cv  <- list()
  no  <- rep(0, K)
  count <- rep(0, K)

  # Count number of observations in each cluster
  for(i in 1:K)
    count[i] <- sum(zRejs == i)

  # Loop update params for each component
  for(i in 1:K) {
    nm <- count[i]
    X  <- matrix(xRejs[, zRejs == i], nrow = dm)

    # Update likelihood mean and S for each component
    if(nm > 0) {
      mn <- rowMeans(X)
      if (dm == 1){
        S <- sum(X^2) #sum of squares
        if (S == Inf)
          S <- 1e05
      } else {
        S <- (X - mn) %*% t(X - mn)
      }
    } else {
      mn <- rep(0, dm)
      if (dm == 1){
        S <- 0
      } else {
        S <- diag(rep(0, dm))
      }
    }

    SS <- list(n = nm, mn = mn, S = S)

    # Update posterior for each component
    if (dm == 1){
      rslt <- nigPost(m0 = niwpr$m0, v0 = niwpr$v0,
                       a0 = niwpr$a0, b0 = niwpr$b0, SS)
    } else {
      rslt <- niwPost(niwpr$mu0, niwpr$lam, niwpr$Phi, niwpr$nu, SS)
    }

    mu[[i]] <- rslt$mu
    cv[[i]] <- rslt$cv
  }

  return(list(wt = wt, mu = mu, cv = cv, no = count))
}

