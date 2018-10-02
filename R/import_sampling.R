#'  Monte Carlo integration function
#' 
#' This function evaluates an integral using Monte Carlo integration. 
#' @param constr a function that describes the domain or constraint of the observations
#' @param mn mean of the observations to be generated
#' @param cv covariance matrix of the observations to be generated
#' @param num_obs number of observations to be generated
#' @return value of integrals
#' @export
#' @examples 
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' mc_inter(constr, 1/2, 0.1, 100)
#'                                 
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
#' @export
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
#' @export
#' @examples 
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' important_sampling_norm(constr, 1/2, 0.1, 100)
#'
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
#' @export
#' @examples 
#' X <- rnorm(5, mean=0, sd=1)
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' dtmvnorm_impsamp(X, mn, sig, constr, num_samp=1e6)
#'
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
#' @export
#' @examples 
#' genconst <- function(u1, u2, l1, l2) { in_set <- function(X) {
#'  if (nrow(X) == 1) X = t(X) 
#'  X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 & X[3,] > l1 & 
#'  X[3,] < u1 & X[4,] > l2 & X[4,] < u2 }}
#' constr  <- genconst(u = 0, u = 0, l1 = 1, l2 = 1)
#' import_sampling_mv(constr, rep(1/2, 4), diag(0.1,4), rep(0,4), diag(0.3,4), 100)
#'  
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
#' @export
#' @examples 
#' genconst <- function(u1, u2, l1, l2) { in_set <- function(X) {
#'  if (nrow(X) == 1) X = t(X) 
#'  X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 & X[3,] > l1 & 
#'  X[3,] < u1 & X[4,] > l2 & X[4,] < u2 }}
#' constr  <- genconst(u1 = 0, u2 = 0, l1 = 1, l2 = 1)
#' dtmvnorm_impsamp_mv(X = mvnfast::rmvn(5, mu = rep(0,4), sigma = diag(0.1, 4)), mn = rep(0, 4), sig = diag(0.1,4), constr, num_samp = 1000, log = F)
#' 
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

# # =========== Test Important Sampling for 1 dimension ============ #
# dm <- 1
# mn <- 0
# cv <- 1
# 
# # Bounds for the domain
# lower <- 2
# upper <- 3
# 
# # Function to eliminate points outside the hypercube
# constr <- function(X, u1 = lower, u2 = upper) { 
#   X >= u1 & X <= u2 
# }
# 
# num_obs  <- 10000
# 
# imp_samp_res <- import_sampling(constr, mn, cv, num_obs)
# print(imp_samp_res)
# 
# # Check result
# pnorm(upper, mn, cv) - pnorm(lower, mn, cv)
# 
# # =========== Use important sampling to compute truncated dist. ============ #
# # Bounds for the domain
# lower <- 0
# upper <- 1
# 
# # Function to eliminate points outside the hypercube
# constr <- function(X, u1 = lower, u2 = upper) { 
#   X >= u1 & X <= u2 
# }
# 
# # Test using tmvtnorm function: 
# test_dat <- rmvn(5, mu = rep(1/2,1), sigma = 0.1)
# 
# test_dense <- dtmvnorm(test_dat, mean = rep(1/2,1), sigma = 0.1, 
#                        lower = rep(lower, 1), upper = rep(upper,1))
# print(test_dense)
# 
# test_is    <- dtmvnorm_impsamp(test_dat, mn = rep(1/2,1), sig = 0.1, 
#                                constr, num_samp = 1e6)
# print(test_is)
# 
# ============= General Important Sampling (Multidimensional) ============== #
# genconst <- function(u1, u2, l1, l2) {
#   in_set <- function(X) {
#     if (nrow(X) == 1) X = t(X)
#     X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2
#   }
# }
# 
# # Bounds for the domain
# lower <- 2
# upper <- 3
# 
# constr  <- genconst(upper, upper, lower, lower)
# 
# mn <- c(0,0)
# cv <- diag(1,2)
# 
# mu_prop  <- c(0,0)
# sig_prop <- diag(4,2)
# 
# num_obs  <- 1e5
# 
# imp_samp_res <- import_sampling_mv(constr, mn, cv, mu_prop, sig_prop, num_obs)
# print(imp_samp_res)
# 
# # Check result
# pmvnorm(rep(lower,2), rep(upper,2), mean = mn, sigma = cv)[1]

# =========== Use important sampling to compute truncated dist. - mv ============ #
# # Function to eliminate points outside the hypercube
# genconst_2 <- function(u1, u2, l1, l2) {
#   in_set <- function(X) {
#     if (nrow(X) == 1) X = t(X)
#     X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2
#   }
# }
# 
# genconst_3 <- function(u1, u2, l1, l2) {
#   in_set <- function(X) {
#     if (nrow(X) == 1) X = t(X)
#     X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 &
#       X[3,] > l1 & X[3,] < u1
#   }
# }
# 
# genconst_4 <- function(u1, u2, l1, l2) {
#   in_set <- function(X) {
#     if (nrow(X) == 1) X = t(X)
#     X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 &
#       X[3,] > l1 & X[3,] < u1 & X[4,] > l1 & X[4,] < u1
#   }
# }

# ================= Test against tmvtnorm function: ================ #
# dm <- 2
# 
# upper <- 1
# lower <- 0
# 
# if (dm == 2)
#   constr  <- genconst_2(upper, upper, lower, lower)
# 
# if (dm == 3)
#   constr  <- genconst_3(upper, upper, lower, lower)
# 
# if (dm == 4)
#   constr  <- genconst_4(upper, upper, lower, lower)
# 
# test_dat <- rmvn(15, mu = rep(1/2, dm), sigma = diag(rep(0.1, dm)))
# test_dense <- dtmvnorm(test_dat, mean = rep(1/2, dm), sigma = diag(rep(0.1, dm)),
#                        lower = rep(0, dm), upper = rep(1, dm))
# print(test_dense)
# 
# test_is    <- dtmvnorm_impsamp_mv(test_dat, mn = rep(1/2, dm), sig = diag(rep(0.1, dm)),
#                                constr, num_samp = 1e6)
# print(test_is)

  


