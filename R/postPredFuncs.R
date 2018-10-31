#' Posterior predictive for univariate MoTG
#'
#' This function evaluates the posterior predictive distribution for a univariate MoTG model
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test test set where columns are the individual observations and rows are the dimension
#' @param pts grid points to draw the density distribution
#' @param upper upper bound of constraints, only use for unit square
#' @param lower lower bound of constraints, only use for unit square
#' @param realDens real density values at points
#' @return list of prediction at points, prediction of tests, along with error bars
#' @export
#' @examples
#'
postPredMotg <- function(res, test, pts, upper, lower, realDens){
  K <- length(res[[1]]$wt)
  iter = length(res)

  grdWtPt <- rep(0, length(pts))
  grdWtTs <- rep(0, length(test))

  predPt <- rep(list(), K)
  predTs <- rep(list(), K)

  predPtStr  <- matrix(0, nrow = iter, ncol = length(pts))
  absErrPreb <- matrix(0, nrow = iter, ncol = length(pts))
  sqErrPreb  <- matrix(0, nrow = iter, ncol = length(pts))

  probsQuant  <- c(0.10, 0.25, 0.5, 0.75, 0.90)

  predPtSumry <- matrix(0, nrow = length(probsQuant), ncol = length(pts))
  absErrSumry <- matrix(0, nrow = length(probsQuant), ncol = length(pts))
  sqErrSumry  <- matrix(0, nrow = length(probsQuant), ncol = length(pts))

  for (i in 1:iter){
    if(i %% 500 == 0) print(i)

    for(k in 1:K) {
      cv <- res[[i]]$cv[[k]][[1]]
      mn <- res[[i]]$mu[[k]]

      predPt[[k]] <- truncnorm::dtruncnorm(pts, a = lower, b = upper,
                                 mean = mn, sd = sqrt(cv))
      predTs[[k]] <- truncnorm::dtruncnorm(test, a = lower, b = upper,
                                 mean = mn, sd = sqrt(cv))

      predPt[[k]] <- editPredProb(predPt[[k]])
      predTs[[k]] <- editPredProb(predTs[[k]])

      grdWtPt <- grdWtPt + res[[i]]$wt[k] * predPt[[k]]
      grdWtTs <- grdWtTs + res[[i]]$wt[k] * predTs[[k]]

      predPtStr[i, ] <- predPtStr[i, ] + res[[i]]$wt[k] * predPt[[k]]
    }

    absErrPreb[i, ] <- abs(predPtStr[i,] - realDens)
    sqErrPreb[i, ]  <- (predPtStr[i,] - realDens)^2

  }

  # Summary statistics of prediction point and errors
  predPtSumry <- apply(predPtStr, 2, function(x) quantile(x, probs = probsQuant))
  predPtStd <- apply(predPtStr, 2, sd)

  absErrSumry <- apply(absErrPreb, 2, function(x) quantile(x, probs = probsQuant))
  sqErrSumry  <- apply(sqErrPreb, 2, function(x) quantile(x, probs = probsQuant))

  absErrStd <- apply(absErrPreb, 2, sd)
  sqErrStd  <- apply(sqErrPreb, 2, sd)

  gridVal    <- list(predPt = grdWtPt/iter,
                    predTs = grdWtTs/iter,
                    predPtSumry = predPtSumry,
                    predPtStd = predPtStd,
                    absErrSumry = absErrSumry,
                    sqErrSumry  = sqErrSumry,
                    absErrStd = absErrStd,
                    sqErrStd  = sqErrStd)

  return(gridVal)

}

#' Posterior predictive for univariate TMoG
#'
#' This function evaluates the posterior predictive distribution for a univariate TMoG model
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test test set where columns are the individual observations and rows are the dimension
#' @param pts grid points to draw the density distribution
#' @param upper upper bound of constraints, only use for unit square
#' @param lower lower bound of constraints, only use for unit square
#' @param realDens real density values at points
#' @return list of prediction at points, prediction of tests, along with error bars
#' @export
postPredTmog <- function(res, test, pts, upper, lower, realDens){
  K <- length(res[[1]]$wt)
  iter <- length(res)

  grdWtPt <- rep(0, length(pts))
  grdWtTs <- rep(0, length(test))

  logPredPtStr <- matrix(0, ncol = length(pts), nrow = iter)
  logPredTsStr <- matrix(0, ncol = length(test), nrow = iter)

  predPtStr <- matrix(0, ncol = length(pts), nrow = iter)
  predTsStr <- matrix(0, ncol = length(test), nrow = iter)

  denomPt <- 0
  denomTs <- 0

  predPtStr  <- matrix(0, nrow = (iter), ncol = length(pts))
  absErrPreb <- matrix(0, nrow = (iter), ncol = length(pts))
  sqErrPreb  <- matrix(0, nrow = (iter), ncol = length(pts))

  probsQuant  <- c(0.10, 0.25, 0.5, 0.75, 0.90)

  predPtSumry <- matrix(0, nrow = length(probsQuant), ncol = length(pts))
  absErrSumry <- matrix(0, nrow = length(probsQuant), ncol = length(pts))
  sqErrSumry  <- matrix(0, nrow = length(probsQuant), ncol = length(pts))

  for (i in 1:iter){
    if(i %% 500 == 0)
      print(i)

    for(k in 1:K) {
      cv <- res[[i]]$cv[[k]][[1]]
      mn <- res[[i]]$mu[[k]]

      logPredPt <- dnorm(pts, mean = mn, sd = sqrt(cv), log = T)
      logPredTs <- dnorm(test, mean = mn, sd = sqrt(cv), log = T)

      grdWtPt <- grdWtPt + res[[i]]$wt[k] * exp(logPredPt)
      grdWtTs <- grdWtTs + res[[i]]$wt[k] * exp(logPredTs)

      denomPt <- denomPt + res[[i]]$wt[k] * normCDF(mn, sqrt(cv), lower, upper)
      denomTs <- denomTs + res[[i]]$wt[k] * normCDF(mn, sqrt(cv), lower, upper)
    }

    if (log(denomTs) == -Inf){
      logPredTsStr[i,] <- 0
    } else {
      logPredTsStr[i,]  <- log(grdWtTs) - log(denomTs)
    }

    if (log(denomPt) == -Inf){
      logPredPtStr[i,] <- 0
    } else {
      logPredPtStr[i,] <- log(grdWtPt) - log(denomPt)
    }

    absErrPreb[i, ] <- abs(exp(logPredPtStr[i,]) - realDens)
    sqErrPreb[i, ]  <- (exp(logPredPtStr[i,]) - realDens)^2

    grdWtPt <- rep(0, length(pts))
    grdWtTs <- rep(0, length(test))

    denomPt <- 0
    denomTs <- 0

  }

  maxLogPredPt <- apply(logPredPtStr, 2, max)
  maxLogPredTs <- apply(logPredTsStr, 2, max)

  logPredPt <- t(t(logPredPtStr) - maxLogPredPt)
  logPredTs <- t(t(logPredTsStr) - maxLogPredTs)

  predPt  <- exp(logPredPt)
  predTs  <- exp(logPredTs)

  logMeanPredPt <- log(apply(predPt, 2, sum)) + maxLogPredPt +
    log(1/(iter))
  logMeanPredTs <- log(apply(predTs, 2, sum)) + maxLogPredTs +
    log(1/(iter))

  # Summary statistics of prediction point and errors
  predPtSumry <- apply(exp(logPredPtStr), 2, function(x) quantile(x, probs = probsQuant))
  predPtStd   <- apply(exp(logPredPtStr), 2, sd)

  absErrSumry <- apply(absErrPreb, 2, function(x) quantile(x, probs = probsQuant))
  sqErrSumry  <- apply(sqErrPreb, 2, function(x) quantile(x, probs = probsQuant))

  absErrStd <- apply(absErrPreb, 2, sd)
  sqErrStd  <- apply(sqErrPreb, 2, sd)

  gridVal <- list(predPt = exp(logMeanPredPt),
                    predTs = exp(logMeanPredTs),
                    predPtSumry = predPtSumry,
                    predPtStd = predPtStd,
                    absErrSumry = absErrSumry,
                    sqErrSumry  = sqErrSumry,
                    absErrStd = absErrStd,
                    sqErrStd  = sqErrStd)

  return(gridVal)

}

#' Posterior predictive for univariate data
#'
#' This function evaluates the posterior predictive distribution for univariate models
#' @param mix one of "tmog" or "motg"
#' @param constrType one of "simple" and "complex". If "simple" include lower and upper bounds.
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test matrix of p-by-N where N are number of observations and p are number of dimension
#' @param pts grid points to draw the density distribution
#' @param realDens real density values at grid points
#' @return list of prediction at points, prediction of tests, along with error bars
#' @export
#' @examples
#'
postPredDist <- function(mix, constrType, res, test, pts, realDens, ...){
  if (constrType == "simple"){
    if (hasArg(upper) == T & hasArg(lower) == T){
      if (mix == "tmog"){
        ppres <- postPredTmog(res, test, pts, upper, lower, realDens)
      } else if (mix == "motg"){
        ppres <- postPredMotg(res, test, pts, upper, lower, realDens)
      } else {
        stop("Mix must be one of 'tmog' and 'motg'")
      }
    } else {
      stop("Must include upper and lower bounds for simple constraint type")
    }
  } else if (constrType == "complex"){
    stop("Constraint type for complex in one dimension is not yet supported")
  } else {
    stop("Constraint type must be one of 'simple' or 'complex'")
  }
  return(ppres)
}

#' Posterior predictive for multivariate MoTG
#'
#' This function evaluates the posterior predictive distribution for a multivariate MoTG model
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test test set where columns are the individual observations and rows are the dimension
#' @param constr a function for the constraints
#' @param thinsc thinning for MCMC chain (default no thinning)
#' @param constrType one of "simple" or "complex". Follow with lower and upper bounds if "simple".
#' @return list of prediction of tests
#' @export
postPredMotgMv <- function(res, test, constr, thinsc = 1, constrType, ...){
  # ================================================================================ #
  # Posterior predictive distribution for Mixture of Truncated Gaussian Approximate  #
  # Only for Multi-dimensional data                                                  #
  # Returns posterior predictive of given points                                     #
  # ================================================================================ #

  K <- length(res[[1]]$wt)
  iter <- length(res)
  dm   <- nrow(test)

  grdWtTs   <- rep(0, ncol(test))
  gridVal   <- rep(0, ncol(test))

  predTsStr <- list(list())
  nonzSamp  <- rep(0, iter)

  if (constrType == "simple"){
    sampMethod = "package"
  } else if (constrType == "complex"){
    sampMethod = "is"
  }

  numObsProp <- 1e3

  for (i in 1:iter){
    if(i %% 100 == 0)
      print(i)

    predTs <- list()

    for(k in 1:K) {

      cv <- res[[i]]$cv[[k]][1:dm, 1:dm]
      mn <- res[[i]]$mu[[k]][1:dm]

      # Use of package tmvtnorm to compute Truncated Multivariate Gaussian density
      if (sampMethod == "package")
        predTs[[k]] <- dtmvnormFast(t(test), mean = as.vector(mn), sigma = cv,
                                         lower = lower, upper = upper, log = T)
      if (sampMethod == "is")
        predTs[[k]] <- dtmvnormMvImportSampling(t(test), mn = as.vector(mn), sig = cv,
                                            constr, numObsProp, log = T)

      predTs[[k]] <- editLogProb(predTs[[k]])

      grdWtTs   <- grdWtTs + res[[i]]$wt[k] * exp(predTs[[k]])

    }

    gridVal <- gridVal + grdWtTs
    grdWtTs <- rep(0, ncol(test))

  }

  return(list(predTs = gridVal/iter))
}

#' Posterior predictive for multivariate TMoG
#'
#' This function evaluates the posterior predictive distribution for a multivariate TMoG model
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test test set where columns are the individual observations and rows are the dimension
#' @param constr a function for the constraints
#' @param thinsc thinning for MCMC chain (default no thinning)
#' @param constrType one of "simple" or "complex". Follow with lower and upper bounds if "simple".
#' @return list of prediction of tests
#' @export
postPredTmogMv <- function(res, test, constr, thinsc = 1, constrType, ...){
  # ================================================================================ #
  # Posterior predictive distribution for Truncated Mixtures of Gaussian             #
  # Only for Multi-dimensional data                                                  #
  # Returns posterior predictive of given points                                     #
  # ================================================================================ #

  K <- length(res[[1]]$wt)
  iter <- length(res)/thinsc
  dm   <- nrow(test)

  grdWtTs <- rep(0, ncol(test))
  denomTs <- 0

  logPredTsStr <- matrix(0, ncol = ncol(test), nrow = iter)
  nonzSamp <- rep(0, iter)

  if (constrType == "simple"){
    sampMethod = "package"
  } else if (constrType == "complex"){
    sampMethod = "is"
  }

  numObsProp <- 1e3

  for (i in 1:iter){
    if(i %% 100 == 0)
      print(paste("i:", i))

    d <- i*thinsc

    for(k in 1:K) {
      cv <- res[[d]]$cv[[k]][1:dm, 1:dm]
      mn <- res[[d]]$mu[[k]][1:dm]

      muProp  <- mn
      sigProp <- cv*4

      # Compute Likelihood of test
      logPredTs <- mvnfast::dmvn(t(test), mu = mn, sigma = cv, log = TRUE)
      grdWtTs   <- grdWtTs + res[[d]]$wt[k] * exp(logPredTs)

      if (sampMethod == "package")
        areaInt <- mvtnorm::pmvnorm(lower, upper, mean = as.vector(mn), sigma = cv)[1]
      if (sampMethod == "is")
        areaInt <- mvImportanceSampling(constr, mn, cv, muProp, sigProp, numObsProp)

      areaInt <- ifelse(areaInt < 0, 0, areaInt)
      denomTs <- denomTs + res[[d]]$wt[k] * areaInt
    }

    if (log(denomTs) == -Inf){
      logPredTsStr[i, ] <- 0
    } else {
      logPredTsStr[i, ] <- log(grdWtTs) - log(denomTs)
    }

    nonzSamp[i] <- ifelse(log(denomTs) == -Inf, 0, 1)
    grdWtTs <- rep(0, ncol(test))
    denomTs <- 0
  }

  maxLogPredTs <- apply(logPredTsStr, 2, max)
  logPredTs    <- t(t(logPredTsStr) - maxLogPredTs)

  predTs  <- exp(logPredTs)
  effIter <- sum(nonzSamp)

  logMeanPredTs <- log(apply(predTs, 2, sum)) + maxLogPredTs +
    log(1/effIter)

  return(list(predTs = exp(logMeanPredTs)))
}

#' Posterior predictive for multivariate data
#'
#' This function evaluates the posterior predictive distribution for univariate models
#' @param mix type of model. One of "tmog" or "motg"
#' @param constrType type of bounds. One of "simple" or "complex". If "simple", include lower and upper bounds for all dimension
#' @param res list of results containing weight (wt), mean (mn), and covariance (cv) matrix for each cluster
#' @param test matrix of p-by-N where p is the number of dimensions and N is the total number of observations
#' @param constr function which inputs test and outputs TRUE/FALSE
#' @param thinsc thinning of chain (default is 1, no thinning)
#' @param ... input lower and upper bounds for "simple" constraint type
#' @return list of prediction at points, prediction of tests, along with error bars
#' @export
#' @examples
#'
postPredDistMv <- function(mix, constrType, res, test, constr, thinsc, ...){
  if (constrType == "simple"){
    if (hasArg(lower) == T & hasArg(upper) == T){
      if (length(lower) != dm | length(upper) != dm){
        stop("Length of lower and upper must equal to dimension of data")
      }
    } else {
      stop("Must specify lower and upper bounds for 'simple' constraint type")
    }

    if (mix == "tmog"){
      ppres <- postPredTmogMv(res, test, constr, thinsc, constrType, lower, upper)
    } else if (mix == "motg"){
      ppres <- postPredMotgMv(res, test, constr, thinsc, constrType, lower, upper)
    } else {
      stop("Mix must be one of tmog or motg")
    }
  } else if (constrType == "complex"){
    if (mix == "tmog"){
      ppres <- postPredTmogMv(res, test, constr, thinsc, constrType, ...)
    } else if (mix == "motg"){
      ppres <- postPredMotgMv(res, test, constr, thinsc, constrType, ...)
    } else {
      stop("Mix must be one of tmog or motg")
    }
  } else {
    stop("Constraint type must be one of 'simple' or 'complex'")
  }

  return(ppres)
}
