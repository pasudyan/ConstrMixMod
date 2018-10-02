checkSymmetricPositiveDefinite_fast <- function (x, name = "sigma") {
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be a symmetric matrix", name))
  }
  if (NROW(x) != NCOL(x)) {
    stop(sprintf("%s must be a square matrix", name))
  }
  if (any(diag(x) <= 0)) {
    stop(sprintf("%s all diagonal elements must be positive",
                 name))
  }
  if (det(x) <= 0) {
    stop(sprintf("%s must be positive definite", name))
  }
}

checkTmvArgs_fast <- function(mean, sigma, lower, upper) {
  if (is.null(lower) || any(is.na(lower)))
    stop(sQuote("lower"), " not specified or contains NA")
  if (is.null(upper) || any(is.na(upper)))
    stop(sQuote("upper"), " not specified or contains NA")
  if (!is.numeric(mean) || !is.vector(mean))
    stop(sQuote("mean"), " is not a numeric vector")
  if (is.null(sigma) || any(is.na(sigma)))
    stop(sQuote("sigma"), " not specified or contains NA")
  if (!is.matrix(sigma)) {
    sigma <- as.matrix(sigma)
  }
  if (NCOL(lower) != NCOL(upper)) {
    stop("lower and upper have non-conforming size")
  }
  checkSymmetricPositiveDefinite_fast(sigma)
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  if (length(lower) != length(mean) || length(upper) != length(mean)) {
    stop("mean, lower and upper must have the same length")
  }
  if (any(lower >= upper)) {
    stop("lower must be smaller than or equal to upper (lower<=upper)")
  }
  cargs <- list(mean = mean, sigma = sigma, lower = lower,
                upper = upper)
  return(cargs)
}

#' Fast density of truncated multivariate Gaussian
#'
#' This function computes the density of a truncated multivariate Gaussian faster than mvtnorm::dtmvnorm and use mvnfast to evaluate the density
#' @param x matrix of observations where the rows are individual observations and the columns are the dimension
#' @param mean mean of the density to be evaluated
#' @param sigma covariance matrix of the density to be evaluated
#' @param lower lower bounds of the truncation
#' @param upper upper bounds of the truncation
#' @param log whether to return log probabilities
#' @return value of density
#' @export
dtmvnorm_fast <- function(x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                           lower = rep(-Inf, length = length(mean)),
                           upper = rep(Inf, length = length(mean)),
                           log = FALSE) {
  cargs <- checkTmvArgs_fast(mean = mean, sigma = sigma, lower = lower,
                        upper = upper)
  mean <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper

  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  T <- nrow(x)
  insidesupportregion <- logical(T)

  for (i in 1:T) {
    insidesupportregion[i] = all(x[i, ] >= lower & x[i, ] <=
                                   upper & !any(is.infinite(x)))
  }

  if (log) {
    nomin   <- mvnfast::dmvn(x, mu = mean, sigma = sigma, log = TRUE)
    denomin <- round(mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma)[1], 5)
    denomin <- ifelse(denomin < 0, 0, denomin)
    log_denomin <- log(denomin)

    dvin  <-  nomin - log_denomin
    dvout <- -Inf

  }

  else {
    nomin   <- round(mvnfast::dmvn(x, mu = mean, sigma = sigma, log = FALSE), 10)
    denomin <- round(mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma), 10)
    denomin <- ifelse(denomin < 0, 0, denomin)

    dvin  <- nomin/denomin
    dvin  <- ifelse(is.nan(dvin), 0, dvin)
    dvin  <- ifelse(dvin == Inf, 0, dvin)

    dvout <- 0
  }

  f <- ifelse(insidesupportregion, dvin, dvout)
  return(f)
}

