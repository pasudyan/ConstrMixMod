#' Rejection summary
#'
#' This function summarizes the value for the rejected samples
#' @param rejs vector of rejected samples
#' @param brks break points of the bins to summarize the rejections
#' @return number of samples in the bins
#'
rejSummary <- function(rejs, brks){
  rejAll <- rejs
  bins    <- table(cut(rejAll, brks, include.lowest = T, right = F))
  return(bins)
}

editProbLogNan <- function(X){
  X[is.nan(X) == T] <- 1e-05
  return(X)
}

is.empty <- function(listOf){
  unlist(lapply(1:length(listOf), function(x){
    if (length(listOf[[x]]) == 0)
      return(TRUE)
    else
      return(FALSE)
  }))
}

#' Counting the rejections
#'
#' This function summarizes the value for the rejected samples
#' @param resList list of the results of the gibbs sampling procedure
#' @return number of total rejections
#'
countRejs <- function(resList){
  unlist(lapply(1:length(resList), function(x){
    sum(resList[[x]]$cnt_rej)
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
#' mcInter(constr, 1/2, 0.1, 100)
#'
normCDF <- function(mean, sd, lower, upper){
  return(pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
}

#' Edit underflow of predictive probability
#'
#' This function edits the underflow of the predictive probability function for truncated Gaussian
#' @param X vector of probabilities
editPredProb <- function(X){
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
editLogProb <- function(X){
  X[X == -Inf] <- -Inf
  X[X == Inf]  <- -Inf
  X[is.nan(X) == T]  <- -Inf
  return(X)
}

#' Mixture of Gaussians
#'
#' This function generates mixtures of gaussians samples
#' @param K number of clusters (default 2)
#' @param paramList list containing weight (wt), mean (mu), and covariance (cv)
#' @param N number of samples
#' @param dm dimensions of Gaussian
#' @return list with Gaussian samples and component assignment
genMog <- function(K = 2, paramList, N = 100, dm){
  wt  <- paramList$wt
  mu  <- paramList$mu
  cv  <- paramList$cv

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
editProb <- function(X){
  X[is.nan(X) == T] <- 0
  X[X == Inf]       <- 0
  return(X)
}

#' Normalize probabilities
#'
#' This function normalizes a probability distribution
#' @param X vector of probabilities
#' @return normalize vector of probabilities
probNorm <- function(X){
  norm <- X - matrixStats::logSumExp(X)
  return(norm)
}

#' Index of rejections
#'
#' This function gives the rejected samples indices associated to the observation it belongs to
#' @param vl a vector of TRUE/FALSE indicating whether the samples are inside or outside the constraints
#' @param cInd number of last accepted observations
giveInd <- function(vl, cind){
  indt   <- which(vl == TRUE)
  inda   <- rep(0, length(vl))
  inda[indt] <- Inf
  for (i in indt){
    inda[i:length(vl)] <- inda[i:length(vl)] + 1
  }

  return(list(indRejObs = inda + cind, cind = length(indt)))
}

#' Index of last acceptance
#'
#' This function provides the index of the observations with that is last accepted
#' @param numA total number of accepted samples
#' @param left number of observations without rejected samples
#' @param vl a vector of TRUE/FALSE indicating whether the samples are inside or outside the constraints
#'
indLastAcc <- function(numA, left, vl){
  # Function to determine index of last few rejections

  temp <- numA - (numA - left)
  ind  <- which(vl == TRUE)[temp]
  tInd <- which(vl[1:ind] == FALSE)
  return(tInd)
}

#' Imputation of rejected samples
#'
#' This function samples the rejected proposals for the Truncated Mixture of Gaussian (TMoG) model
#' @param K total number of clusters
#' @param params list of parameters which includes the weight (wt), mean (mu), and covariance matrix (cv)
#' @param N total number of observations
#' @param constr function for the constraint of the space
#' @param thr maximum number of rejections per observation
#' @param ceilRej max number of observations in total
#' @return rejected samples and indices to which observation it belongs to
#'
imputeRejections <- function(K, params, N, constr, thr, ceilRej) {

  dm <- dim(params$cv[[1]])[1]

  if(is.null(dm))
    dm <- 1

  # Number of cumulative acceptance
  cumAcc <- 0

  # Number of rejected samples
  rejAcc <- 0

  rejs <- matrix(rep(0, 0), nrow = length(params$mu[[1]]))
  z    <- c()
  left <- N
  cind <- 1
  indRejObs <- c()

  while (left > 0) {
    if (thr == 0){
      left <- 0
      return(list(rejs = rejs, z = z, rejAcc = rejAcc))
    }

    # Sample from a mixture of gaussian
    smpls  <- genMog(K, params, N, dm)
    vl     <- constr(smpls$x)

    # Save accepted indices
    tmp  <- giveInd(vl, cind)
    cind <- cind + tmp$cind

    # Number of accepted samples on this iteration
    numA <- sum(vl)

    # Number of rejected samples on this iteration
    numR <- N - numA

    if (thr == Inf){
      if (numA <= left){
        cumAcc <- cumAcc + numA
        rejAcc <- rejAcc + numR

        rejs <- cbind(rejs, matrix(smpls$x[,!vl], nrow = dm))
        z    <- c(z, smpls$z[!vl])
        left <- N - cumAcc

        indRejObs <- c(indRejObs, tmp$indRejObs)

        if (left == 0)
          return(list(rejs = rejs, z = z,
                      indRejObs = indRejObs[(indRejObs != Inf) & (indRejObs <= N)],
                      cumAcc = cumAcc, cind = cind, rejAcc = rejAcc))

        if (rejAcc >= ceilRej)
          return(list(rejs = rejs, z = z,
                      indRejObs = indRejObs[(indRejObs != Inf) & (indRejObs <= N)],
                      cumAcc = cumAcc, cind = cind, rejAcc = rejAcc))


      } else {

        tInd <- indLastAcc(numA, left, vl)

        rejs <- cbind(rejs, matrix(smpls$x[, tInd], nrow = dm))
        z    <- c(z, smpls$z[tInd])
        rejAcc <- ncol(rejs)

        indRejObs <- c(indRejObs, tmp$indRejObs[tInd])

        cumAcc <- min(cumAcc, N)

        return(list(rejs = rejs, z = z, indRejObs = indRejObs[indRejObs != Inf],
                    cumAcc = cumAcc, cind = cind, rejAcc = rejAcc))
      }

    } else {

      cumAcc <- cumAcc + numA
      rejAcc <- rejAcc + numR

      if (cumAcc >= N & rejAcc < thr*N){
        if (rejAcc == 0){

          return(list(rejs = rejs, z = z, indRejObs = indRejObs[indRejObs != Inf],
                      rejAcc = rejAcc, cumAcc = min(cumAcc, N)))
        } else {

          tInd <- indLastAcc(numA, left, vl)

          rejs <- cbind(rejs, matrix(smpls$x[, tInd], nrow = dm))
          z    <- c(z, smpls$z[tInd])
          rejAcc <- ncol(rejs)
          indRejObs <- c(indRejObs, tmp$indRejObs[tInd])

          left   <- 0
          cumAcc <- min(cumAcc, N)

          return(list(rejs = rejs, z = z, indRejObs = indRejObs[indRejObs != Inf],
                      rejAcc = rejAcc, cumAcc = cumAcc))
        }

      } else if (cumAcc < N & rejAcc >= thr*N){

        lenY <- floor(thr*N)
        rejs <- cbind(rejs, matrix(smpls$x[,!vl], nrow = dm))
        rejs <- matrix(rejs[, 1:lenY], nrow = dm)
        z    <- c(z, smpls$z[!vl])

        indRejObs <- c(indRejObs, tmp$indRejObs)
        indRejObs <- indRejObs[indRejObs != Inf]

        left   <- 0
        cumAcc <- min(cumAcc, N)

        return(list(rejs = rejs, z = z[1:lenY], indRejObs = indRejObs[1:lenY],
                    rejAcc = lenY, cumAcc = cumAcc))

      } else if (cumAcc >= N & rejAcc >= thr*N){

        tInd   <- indLastAcc(numA, left, vl)
        rejs   <- cbind(rejs, matrix(smpls$x[,tInd], nrow = dm))
        rejAcc <- ncol(rejs)
        z      <- c(z, smpls$z[tInd])

        indRejObs <- c(indRejObs, tmp$indRejObs[tInd])
        indRejObs <- indRejObs[indRejObs != Inf]

        lenY <- min(rejAcc, floor(thr*N))
        rejs <- matrix(rejs[, 1:lenY], nrow = dm)

        left <- 0
        cumAcc <- min(cumAcc, N)

        return(list(rejs = rejs, z = z[1:lenY], indRejObs = indRejObs[1:lenY],
                    rejAcc = lenY, cumAcc = cumAcc))
      } else {
        left <- N - cumAcc
        rejs <- cbind(rejs, matrix(smpls$x[,!vl], nrow = dm))
        z    <- c(z, smpls$z[!vl])
        indRejObs <- c(indRejObs, tmp$indRejObs)
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
#' @param ceilRej max number of observations in total
#' @return rejected samples and indices to which observation it belongs to
#'
imputeGauss <- function(mu, cv, N, constr, thr, ceilRej) {

  dm <- dim(cv)[1]

  if(is.null(dm))
    dm <- 1

  # Number of cumulative acceptance
  cumAcc <- 0

  # Number of rejected samples
  rejAcc <- 0

  # Store index of acceptance
  indRejObs <- c()
  rejs      <- matrix(rep(0, 0), nrow = length(mu))

  left <- N
  cind <- 1

  while (left > 0) {

    if (dm == 1){
      # Sample from a univariate normal distribution
      smpls <- matrix(rnorm(N, mu, sqrt(cv[[1]])), nrow = dm, ncol = N)
    } else {
      # Sample from a multivariate normal distribution with the component parameters
      smpls <- t(mvnfast::rmvn(N, mu, cv))
    }

    # Check if samples satisfies constraint
    vl <- constr(smpls)

    # Save accepted indices
    tmp  <- giveInd(vl, cind)
    cind <- cind + tmp$cind

    # Number of accepted samples on this iteration
    numAcc  <- sum(vl)

    # Number of rejected samples on this iteration
    numRejs <- N - numAcc

    if (thr == Inf){

      if (numAcc <= left){
        # If the number of acceptance is still less than how much acceptance
        # needed
        cumAcc <- cumAcc + numAcc # total number of accepted samples
        rejAcc <- rejAcc + numRejs # total number of rejected samples
        rejs      <- cbind(rejs, matrix(smpls[ , !vl], nrow = dm))
        indRejObs <- c(indRejObs, tmp$indRejObs)
        left      <- N - cumAcc

        if (left == 0)
          return(list(rejs = rejs,
                      indRejObs = indRejObs[(indRejObs != Inf) & (indRejObs <= N)],
                      cumAcc = cumAcc, cind = cind))

        if (rejAcc >= ceilRej){
          return(list(rejs = rejs,
                      indRejObs = indRejObs[(indRejObs != Inf) & (indRejObs <= N)],
                      cumAcc = cumAcc, cind = cind))
        }

      } else {
        # If the number of acceptance is greater than how much acceptance is
        # still needed

        tInd   <- indLastAcc(numAcc, left, vl)
        cumAcc <- cumAcc + numAcc

        rejs      <- cbind(rejs, matrix(smpls[, tInd], nrow = dm))
        indRejObs <- c(indRejObs, tmp$indRejObs[tInd])
        rejAcc    <- ncol(rejs)

        cumAcc <- min(cumAcc, N)

        return(list(rejs = rejs, indRejObs = indRejObs[indRejObs != Inf],
                    cumAcc = cumAcc, cind = cind))
      }

    } else {
      # If threshold is not infinity

      cumAcc <- cumAcc + numAcc
      rejAcc <- rejAcc + numRejs

      if (cumAcc >= N & rejAcc < thr*N){
        # If acceptance is more than total number of observations,
        # but number of rejections are less than threshold

        if (rejAcc == 0){
          # Keep all rejections
          return(list(rejs = rejs, indRejObs = indRejObs[indRejObs != Inf],
                      cumAcc = min(cumAcc, N)))
        } else {

          tInd  <- indLastAcc(numAcc, left, vl)
          rejs  <- cbind(rejs, matrix(smpls[, tInd], nrow = dm))
          indRejObs <- c(indRejObs, tmp$indRejObs[tInd])
          rejAcc    <- ncol(rejs)

          left   <- 0
          cumAcc <- min(cumAcc, N)

          # Keep all rejections
          return(list(rejs = rejs, indRejObs = indRejObs[indRejObs != Inf],
                      cumAcc = cumAcc))
        }

      } else if (cumAcc < N & rejAcc >= thr*N){
        # If acceptance is less than total number of observations,
        # but number of rejections are more than threshold

        lenY  <- floor(thr*N)
        rejs  <- cbind(rejs, matrix(smpls[, !vl], nrow = dm))
        indRejObs <- c(indRejObs, tmp$indRejObs)
        indRejObs <- indRejObs[indRejObs != Inf]

        left   <- 0
        cumAcc <- min(cumAcc, N)

        return(list(rejs = matrix(rejs[, 1:lenY], nrow = dm),
                    indRejObs = indRejObs[1:lenY],
                    cumAcc = cumAcc))

      } else if (cumAcc >= N & rejAcc >= thr*N){
        # If acceptance is more than total number of observations,
        # and number of rejections are also more than threshold

        tInd <- indLastAcc(numAcc, left, vl)
        rejs <- cbind(rejs, matrix(smpls[, tInd], nrow = dm))

        indRejObs <- c(indRejObs, tmp$indRejObs[tInd])
        indRejObs <- indRejObs[indRejObs != Inf]
        rejAcc <- ncol(rejs)

        lenY <- min(rejAcc, floor(thr*N))
        left <- 0

        cumAcc <- min(cumAcc, N)

        return(list(rejs = matrix(rejs[, 1:lenY], nrow = dm),
                    indRejObs = indRejObs[1:lenY],
                    cumAcc = cumAcc))

      } else {
        left <- N - cumAcc
        rejs <- cbind(rejs, matrix(smpls[ , !vl], nrow = dm))
        indRejObs <- c(indRejObs, tmp$indRejObs)
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
#' @param numObs number of observations to be generated
#' @return value of integrals
#' @examples
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' mcInter(constr, 1/2, 0.1, 100)
mcInter <- function(constr, mn, cv, numObs){
  # Monte Carlo integration:
  rsamp    <- matrix(mvnfast::rmvn(numObs, mn, cv), nrow = 1)
  rsampIn <- rsamp[constr(rsamp)]

  mc_int   <- (1/numObs)*length(rsampIn)
  return(mc_int)
}

#' Important Sampling function
#'
#' This function evaluates integrals using an important sampling method with a Gaussian distribution proposal
#' @param constr a function that describes the domain or constraint of the observations
#' @param mn mean of the observations to be generated
#' @param cv covariance matrix of the observations to be generated
#' @param numObs number of observations to be generated
#' @return value of integrals
#' @examples
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' important_sampling(constr, 1/2, 0.1, 100)
#'
importanceSampling <- function(constr, mn, cv, numObs){
  # Parameters of proposal - q
  muProp <- 0
  cvProp <- 4

  # Sample from standard normal multivariate distribution - q
  rsamp   <- t(mvnfast::rmvn(numObs, muProp, cvProp))

  # Compute log likelihoods for weight
  rsampPdf  <- mvnfast::dmvn(t(rsamp), mu = mn, sigma = cv, log = T)
  rsampProp <- mvnfast::dmvn(t(rsamp), mu = muProp, sigma = cvProp, log = T)

  # Compute log weights
  logWghts <- rsampPdf - rsampProp

  rsampIn <- constr(rsamp)

  imps <- (1/numObs)*sum(exp(logWghts[rsampIn]))

  return(imps)
}

#' Normalized Important Sampling function
#'
#' This function evaluates integrals using a normalized important sampling method with a Gaussian distribution proposal
#' @param constr a function that describes the domain or constraint of the observations
#' @param mn mean of the observations to be generated
#' @param cv covariance matrix of the observations to be generated
#' @param numObs number of observations to be generated
#' @return value of integrals
#' @examples
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' important_sampling_norm(constr, 1/2, 0.1, 100)
importSamplingNorm <- function(constr, mn, cv, numObs){
  # Normalized important sampling
  # Parameters of proposal - q
  muProp <- 0
  cvProp <- 4

  # Sample from standard normal multivariate distribution - q
  rsamp   <- t(mvnfast::rmvn(numObs, muProp, cvProp))

  # Compute log likelihoods for weight
  rsampPdf  <- mvnfast::dmvn(t(rsamp), mu = mn, sigma = cv, log = T)
  rsampProp <- mvnfast::dmvn(t(rsamp), mu = muProp, sigma = cvProp, log = T)

  # Compute log weights
  logWghts <- rsampPdf - rsampProp

  rsampIn <- constr(rsamp)

  norm <- (1/numObs)*sum(exp(logWghts))
  imps <- ((1/numObs)*sum(exp(logWghts[rsampIn])))/norm

  return(imps)
}

#' Truncated univariate Gaussian density
#'
#' This function evaluates the truncated univariate Gaussian density using an important sampling procedure as the normalizing constant
#' @param X a vector of observations where each column is the independent observations
#' @param mn mean of the Gaussian density
#' @param sig covariance matrix of the Gaussian density
#' @param constr a function that describes the domain or constraint of the samples
#' @param numSamp number of samples in the important sampling procedure
#' @return value of the density
#' @examples
#' X <- rnorm(5, mean=0, sd=1)
#' constr <- function(X, u1 = lower, u2 = upper) { X >= u1 & X <= u2 }
#' dtmvnormImportSampling(X, mn, sig, constr, numSamp=1e6)
dtmvnormImportSampling <- function(X, mn, sig, constr, numSamp = 1e6){
  numer <- dmvn(X, mu = mn, sigma = sig, log = T)
  denom <- log(importanceSampling(constr, mn = mn, cv = sig, numObs = numSamp))
  val   <- numer - denom

  val   <- ifelse(constr(X), val, -Inf)

  return(exp(val))
}

#' Important Sampling for Multivariate variables
#'
#' This function evaluates integrals of a multivariate random variables using important sampling method with a multivariate Gaussian proposal
#' @param constr a function that describes the domain or constraint of the observations
#' @param mn mean of the observations to be generated
#' @param cv covariance matrix of the observations to be generated
#' @param muProp mean of the proposal distribution
#' @param sigProp covariance matrix of the proposal distribution
#' @param numObs number of observations to be generated
#' @return value of integrals
#' @examples
#' genconst <- function(u1, u2, l1, l2) { in_set <- function(X) {
#'  if (nrow(X) == 1) X = t(X)
#'  X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 & X[3,] > l1 &
#'  X[3,] < u1 & X[4,] > l2 & X[4,] < u2 }}
#' constr  <- genconst(u = 0, u = 0, l1 = 1, l2 = 1)
#' mvImportanceSampling(constr, rep(1/2, 4), diag(0.1,4), rep(0,4), diag(0.3,4), 100)
mvImportanceSampling <- function(constr, mn, cv, muProp, sigProp, numObs){

  dm <- length(mn)

  # Parameters of proposal - q
  cvProp <- sigProp

  # Sample from standard normal multivariate distribution - q
  rsamp   <- t(mvnfast::rmvn(numObs, muProp, cvProp))

  # Compute log likelihoods for weight
  rsampPdf  <- mvnfast::dmvn(t(rsamp), mu = mn, sigma = cv, log = T)
  rsampProp <- mvnfast::dmvn(t(rsamp), mu = muProp, sigma = cvProp, log = T)

  # Compute log weights
  logWghts <- rsampPdf - rsampProp

  rsampIn <- constr(rsamp)

  imps <- (1/numObs)*sum(exp(logWghts[rsampIn]))

  return(imps)
}

#' Truncated Multivariate Gaussian density
#'
#' This function evaluates the truncated multivariate Gaussian density using an important sampling procedure as the normalizing constant
#' @param X a vector of observations where each row is the dimension and each column is the independent observations
#' @param mn mean of the multivariate Gaussian density
#' @param sig covariance matrix of the multivariate Gaussian density
#' @param constr a function that describes the domain or constraint of the samples
#' @param numSamp number of samples in the important sampling procedure
#' @param log whether log value is returned (default TRUE)
#' @return value of density
#' @examples
#' genconst <- function(u1, u2, l1, l2) { in_set <- function(X) {
#'  if (nrow(X) == 1) X = t(X)
#'  X[1,] > l1 & X[1,] < u1 & X[2,] > l2 & X[2,] < u2 & X[3,] > l1 &
#'  X[3,] < u1 & X[4,] > l2 & X[4,] < u2 }}
#' constr  <- genconst(u1 = 0, u2 = 0, l1 = 1, l2 = 1)
#' dtmvnormMvImportSampling(X = mvnfast::rmvn(5, mu = rep(0,4), sigma = diag(0.1, 4)), mn = rep(0, 4), sig = diag(0.1,4), constr, numSamp = 1000, log = F)
dtmvnormMvImportSampling <- function(X, mn, sig, constr, numSamp = 10000, log = T){
  numer <- dmvn(X, mu = mn, sigma = sig, log = T)

  meanImp <- mn
  sigImp  <- sig*4

  denom <- log(mvImportanceSampling(constr, mn = mn, cv = sig, muProp = meanImp,
                                    sigProp = sigImp, numObs = numSamp))
  val <- numer - denom
  val <- ifelse(constr(t(X)), val, -Inf)

  if (log)
    return(val)
  else
    return(exp(val))
}

#' Gibbs sampling algorithm
#'
#' This computes the Gibbs sampling procedure
#' @param K maximum number of clusters
#' @param ipdata data matrix where columns are individual observations and rows are dimensions
#' @param thr maximum number of rejections per observations
#' @param burnIn number of burn-in for the chain. Must be less than number of iterations
#' @param numIter total number of iterations
#' @param constr function for the constraints
#' @param mix one of "tmog" and "motg"
#' @param niwpr list of NIW prior parameters
#' @return a list of results consisting of weights, mean, and covariance matrix for every iterations, along with the rejected proposals
#' @export
#'
gibbSampling <- function(K, ipdata, thr = Inf, burnIn, numIter, constr, mix, niwpr) {

  # Dimension and number of observations
  dm  <- dim(ipdata)[1]
  N   <- dim(ipdata)[2]

  # Number of iteration and alpha for dirichlet
  alp      <- 1

  # Initialize clusters and parameters for each cluster
  z       <- kmeans(t(ipdata), 5, nstart = 10, iter.max = 50)$cluster
  params  <- updateParams(K, alp, ipdata, niwpr, z)

  # Number of Rejections per Iteration
  rejSamp <- rep(list(list()), K)
  params$rejSamp <- rejSamp

  # Store result
  rslt    <- rep(list(params), (numIter-burnIn))

  # *********************** Gibbs Sampling loop ************************* #
  for(iter in 1:numIter) {
    if(iter%%100 == 0)
      print(paste("Iter", iter))

    # Reorder observations
    neworder <- sample(ncol(ipdata), replace = F)
    ipdata   <- matrix(ipdata[, neworder], nrow = dm)
    z <- z[neworder]

    if (mix == "motg"){
      # Update parameters for each cluster by first sampling the rejections
      for(i in 1:K) {
        pm <- updateParamsComponent(X  = matrix(ipdata[, z == i], nrow = dm),
                                    mu = params$mu[[i]], cv = params$cv[[i]],
                                    niwpr, constr, thr)

        # Store parameters
        params$mu[[i]] <- pm$mu
        params$cv[[i]] <- pm$cv
        params$rejSamp[[i]] <- pm$rejSamp
        params$cntRej[i]    <- pm$cntRej
        params$no[i]        <- pm$no
      }

      # Update assignments
      z  <- updateClustAssignMotg(K, ipdata, z, rej = params$rejSamp, N, params, thr)

      # Store new cluster parameters
      params$z  <- z

      # Update weights
      newWghts  <- updateWeights(z, K, alp)
      params$wt <- newWghts

      if (iter > burnIn)
        rslt[[iter - burnIn]] <- params

      # Set rejection storage to empty again
      rejSamp <- rep(list(list()), K)

    } else if (mix == "tmog"){
      # Sample rejections
      smplRejs <- imputeRejections(K, params, N, constr, thr, ceilRej)

      # Function to update params
      params <- updateParamsTmog(K, alp, xRejs = cbind(ipdata, smplRejs$rejs),
                                  params$wt, niwpr, zRejs = c(z, smplRejs$z))
      # Update weights
      params$wt <- updateWeights(c(z, smplRejs$z), K, alp)

      # Count number of rejections in partition
      if (iter > burnIn)
        rslt[[iter - burnIn]]$rejSamp <- smplRejs$rejs

      # Update assignments
      z <- updateClustAssignTmog(K, ipdata, N, params)

      if (iter > burnIn){
        rslt[[iter - burnIn]] <- params
        rslt[[iter - burnIn]]$cntRej <- smplRejs$rejAcc
      }
    }

  }

  return(result = rslt)

}

checkSymmetricPositiveDefiniteFast <- function (x, name = "sigma") {
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

checkTmvArgsFast <- function(mean, sigma, lower, upper) {
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
  checkSymmetricPositiveDefiniteFast(sigma)
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

dtmvnormFast <- function(x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                          lower = rep(-Inf, length = length(mean)),
                          upper = rep(Inf, length = length(mean)),
                          log = FALSE) {
  cargs <- checkTmvArgsFast(mean = mean, sigma = sigma, lower = lower,
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





