#'  Normal-Inverse-Gamma function
#'
#' This function generates samples from a Normal-Inverse-Gamma density
#' @param m0 location parameter for the Normal distribution (real)
#' @param v0 scale parameter for variance of the Normal distribution (real > 0)
#' @param a0 shape parameter for the Inverse-Gamma distribution (real > 0)
#' @param b0 rate parameter for the Inverse-Gamma distribution (real > 0)
#' @return list with mean and variance
#' @examples
#' nig(0, 2, 2, 1)
#'
nig <- function(m0, v0, a0, b0) {
  cv <- MCMCpack::rinvgamma(1, a0, b0) # Use rate model instead of scale
  mu <- mvnfast::rmvn(1, m0, cv*v0)

  return(list(mu = mu, cv = cv))
}

#'  Normal-Inverse-Gamma posterior density function
#'
#' This function generates samples from the posterior density of a Normal-Inverse-Gamma funcion
#' @param m0 location parameter for the Normal distribution (real)
#' @param v0 scale parameter for variance of the Normal distribution (real > 0)
#' @param a0 shape parameter for the Inverse-Gamma distribution (real > 0)
#' @param b0 rate parameter for the Inverse-Gamma distribution (real > 0)
#' @param SS a list which includes list(n = number of observations, mn = mean, S = covariance)
#' @return a list with mean and variance
#' @examples
#' nigPost(0, 1, 2, 3, list(n = 2, mn = 1, S = 0.1))
#'
nigPost <- function(m0, v0, a0, b0, SS) {
  n  <- SS$n
  mn <- SS$mn
  S  <- SS$S  #S = t(X - mn) %*% (X - mn)

  vn <- v0/(1 + n*v0)
  an <- a0 + n/2

  mu <- vn*((1/v0)*m0 + n*mn)
  bn <- b0 + 0.5*(m0^2/v0 + S - mu^2/vn)

  return(nig(mu, vn, an, bn))
}

#'  Normal-Inverse-Wishart function
#'
#' This function generates samples from a Normal-Inverse-Wishart density
#' @param m0 location parameter for the Gaussian distribution
#' @param lam scale parameter for the variance of the Gaussian distribution
#' @param Phi a positive definite scale matrix for the Inverse-Wishart distribution
#' @param nu degrees of freedom parameter for the Inverse-Whisart distribution
#' @return list with mean and covariance matrix
#' @examples
#' niw(rep(1/2,2), 1, diag(0.1, 2), 2)
#'
niw <- function(mu0, lam, Phi, nu) {

  cv <- MCMCpack::riwish(nu, Phi)
  mu <- mvnfast::rmvn(1, mu0, cv / lam)

  return(list(mu = mu, cv = cv))
}

#'  Normal-Inverse-Wishart posterior density function
#'
#' This function generates samples from the posterior density of a Normal-Inverse-Wishart funcion
#' @param mu0 location parameter for the Gaussian distribution
#' @param lam scale parameter for variance of the Gaussian distribution
#' @param Phi a positive definite scale matrix for the Inverse-Wishart distribution
#' @param nu degrees of freedom parameter for the Inverse-Whisart distribution
#' @param SS a list which includes list(n = number of observations, mn = mean, S = covariance matrix)
#' @return a list with mean and covariance matrix
#' @examples
#' SS = list(n = 2, mn = rep(0.7, 2), S = diag(0.1,2))
#' niwPost(rep(1/2,2), 1, diag(0.1, 2), 2, SS)
#'
niwPost <- function(mu0, lam, Phi, nu, SS) {
  n    <- SS$n
  mn   <- SS$mn
  S    <- SS$S     # S = t(X - mn) %*% (X - mn)
  muN  <- (lam*mu0 + n*mn)/(lam + n)
  Phin <- Phi + S + ((lam*n)/(lam + n))*outer(mn - mu0, mn - mu0)

  return(niw(muN, lam + n, Phin, nu + n))
}
