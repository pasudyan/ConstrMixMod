#'  Normal-Inverse-Gamma function
#' 
#' This function generates samples from a Normal-Inverse-Gamma density
#' @param m0 location parameter for the Normal distribution (real)
#' @param v0 scale parameter for variance of the Normal distribution (real > 0)
#' @param a0 shape parameter for the Inverse-Gamma distribution (real > 0)
#' @param b0 rate parameter for the Inverse-Gamma distribution (real > 0)
#' @return list with mean and variance 
#' @export
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
#' @export
#' @examples 
#' nig.post(0, 1, 2, 3, list(n = 2, mn = 1, S = 0.1))
#'
nig.post <- function(m0, v0, a0, b0, SS) {
  n     <- SS$n
  mn    <- SS$mn
  S     <- SS$S  #S = t(X - mn) %*% (X - mn)
  
  v_n   <- v0/(1 + n*v0)
  a_n   <- a0 + n/2
  
  mu_n  <- v_n*((1/v0)*m0 + n*mn)
  b_n   <- b0 + 0.5*(m0^2/v0 + S - mu_n^2/v_n)
  
  return(nig(mu_n, v_n, a_n, b_n))
}

#'  Normal-Inverse-Wishart function
#' 
#' This function generates samples from a Normal-Inverse-Wishart density 
#' @param m0 location parameter for the Gaussian distribution 
#' @param lam scale parameter for the variance of the Gaussian distribution
#' @param Phi a positive definite scale matrix for the Inverse-Wishart distribution 
#' @param nu degrees of freedom parameter for the Inverse-Whisart distribution
#' @return list with mean and covariance matrix
#' @export
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
#' @export
#' @examples 
#' SS = list(n = 2, mn = rep(0.7, 2), S = diag(0.1,2))
#' niw.post(rep(1/2,2), 1, diag(0.1, 2), 2, SS)
#'
niw.post <- function(mu0, lam, Phi, nu, SS) {
  n     <- SS$n
  mn    <- SS$mn
  S     <- SS$S     # S = t(X - mn) %*% (X - mn)
  mu_n  <- (lam*mu0 + n*mn)/(lam + n)
  Phi_n <- Phi + S + ((lam*n)/(lam + n))*outer(mn - mu0, mn - mu0)
  
  return(niw(mu_n, lam + n, Phi_n, nu + n))
}

# ************ Test for 1-dimensional ************ #
# m <- 0
# l <- 2
# a <- 2
# b <- 5/100
# num_obs = 10000
# temp_mu <- rep(0, num_obs)
# temp_cv <- rep(0, num_obs)
# 
# for (i in 1:num_obs){
#   temp_mu[i] <- nig(m, l, a, b)$mu
#   temp_cv[i] <- nig(m, l, a, b)$cv
# }
# 
# mean(temp_mu)
# sd(temp_cv)
# mean(sqrt(temp_cv))
# sd(sqrt(temp_cv))

# # ************ Test for Multi-dimensional ************ #
# dm <- 2
# mu0 <- rep(1, dm)
# lam <- 1e-2
# Phi <- diag(rep(1e-3, dm))
# nu  <- dm+2
# 
# numm <- 100000
# temp_mu <- matrix(0, ncol = dm, nrow = numm)
# temp_cv <- matrix(0, ncol = dm, nrow = numm)
# 
# for (i in 1:numm){
#   temp <- niw(mu0 = mu0, lam, Phi, nu)
#   temp_mu[i,] <- temp$mu
#   temp_cv[i,] <- c(temp$cv[1,1],
#                    temp$cv[2,2])
# }
# 
# # Mean of means
# colMeans(temp_mu)
# # Var of means
# apply(temp_mu, 2, var)
# # SD of means
# apply(temp_mu, 2, sd)
# 
# # Means of SD
# colMeans(sqrt(temp_cv))
# 
# # Means of Var
# colMeans(temp_cv)
# 
# # Var of SD
# apply(sqrt(temp_cv), 2, var)
# # SD of SD
# apply(sqrt(temp_cv), 2, sd)
# apply(temp_cv, 2, var)
# apply(temp_cv, 2, sd)

