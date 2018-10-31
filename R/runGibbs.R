#' Main function for inference
#'
#' This function evaluates the cumulative distribution of a univariate Gaussian distribution
#' @param mix Mixture type: one of "tmog" and "motg"
#' @param Xdata p-by-N matrix where N is the total number of observations and p is the dimension
#' @param constr function to test whether observation belongs in a space. Must input Xdata and returns vector of TRUE/FALSE
#' @param K maximum number of clusters
#' @param thr vector of number of rejections per observations.
#' @param constrType one of "simple" and "complex". If "simple", follow by lower and upper bounds.
#' @param niwpr list of niw prior. For one dimensional, must include m0, v0, a0, b0. For more than one dimension, must include mu0, lam, Phi, and nu.
#' @param ... pass lower and upper bounds of the constraint if type is "simple"
#' @return samples from posterior distribution
#' @export
#'
infConstrMixMod <- function(mix, Xdata, constr, K = 100, thr = Inf, constrType, niwpr, ...){

  # Check correct inputs
  if (!(mix %in% c("tmog", "motg"))){
    stop("Mixture type must be one of 'tmog' or 'motg'")
  }

  if (thr < 0){
    stop("Value of thr must be greater than 0")
  }

  if (K < 0){
    stop("Value of K must be greater than 0")
  }

  if (!is.matrix(Xdata)){
    stop("Data must be in a p-by-N matrix format, where p is the number of dimensions and
         N is the number of observations")
  }

  dm <- dim(Xdata)[1]

  if (constrType == "simple"){
    sampMethod <- "package"
    if (hasArg(lower) == T & hasArg(upper) == T){
      if (length(lower) != dm & length(upper) != dm){
        stop("Length of lower and upper must equal dimension of data")
      }
    } else {
      stop("Must input lower and upper arguments for simple constraint type")
    }
  } else {
    sampMethod <- "is"
  }

  burnIn  <- 2000
  numIter <- 5000

  # Maximum number of rejections
  ceilRej <- 100000

  # Thinning value for the MCMC samples on the posterior
  thinsc <- 1
  thrl   <- length(thr)

  # Initialize variables
  rslt <- rep(list(list()), thrl)

  print(paste("Mix:", mix))

  # ============= Run Gibbs Sampling and Posterior Probabilities ============= #
  # Loop over all thresholds
  for (j in 1:thrl){
    print(paste("Running Gibbs Sampling for Threshold: ", thr[j]))

    # Function to run Gibbs sampling and timed
    rslt[[j]] <- gibbSampling(K, ipdata = Xdata, thr = thr[j], burnIn, numIter,
                              constr, mix, niwpr)

    print("*******************************")
  }

  return(result = rslt)
}







