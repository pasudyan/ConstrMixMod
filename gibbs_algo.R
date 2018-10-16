#' Gibbs sampling algorithm
#'
#' This computes the Gibbs sampling procedure
#' @param K maximum number of clusters
#' @param ip_data data matrix where columns are individual observations and rows are dimensions
#' @param thr maximum number of rejections per observations
#' @param burnIn number of burn-in for the chain. Must be less than number of iterations
#' @param numIter total number of iterations
#' @param constr function for the constraints
#' @param mix one of "tmog" and "motg"
#' @param lower lower constraint for unit square
#' @param upper upper constriant for unit square
#' @param data_type one of "gauss", "beta", "flowcyto", "crime", "bvgauss"
#' @return a list of results consisting of weights, mean, and covariance matrix for every iterations, along with the rejected proposals
#' @export
#'
gibbs_sampling <- function(K, ip_data, thr = Inf, burnIn, numIter, constr, mix, lower, upper, data_type) {

  # Dimension and number of observations
  dm      <- dim(ip_data)[1]
  N       <- dim(ip_data)[2]

  if (dm == 1)
    mu_p    <- mean(ip_data)

  # Number of iteration and alpha for dirichlet
  alp      <- 1

  # NIW Prior
  if (dm == 1){
    if (data_type == "gauss"){
      niw_p  <- list(m0 = mu_p, v0 = 2,
                     a0 = 2,  b0 = 5/100)
    } else if (data_type == "beta"){
      niw_p  <- list(m0 = mu_p, v0 = 2,
                     a0 = 2,  b0 = 1/100)
    } else {
      stop("Data type for one-dimension is not one of Gaussian or Beta")
    }

  } else {
    if (data_type == "flowcyto"){
      niw_p  <- list(mu0 = rep(1/2, dm), lam = 1e-2,
                     Phi = diag(rep(1e-3, dm)),  nu = dm+1)
    } else if (data_type == "bvgauss"){
      niw_p  <- list(mu0 = 1/2, lam = 1/10,
                     Phi = diag(rep(1e-3, dm)),  nu = dm+1)
    } else if (data_type == "crime"){
      niw_p  <- list(mu0 = 1/2, lam = 1/10,
                     Phi = diag(rep(1e-3, dm)),  nu = dm+1)
    } else {
      stop("Data type for multi-dimensions is not one of Flowcytometry, BVgauss, or Crime")
    }
  }

  # Initialize clusters and parameters for each cluster
  z       <- kmeans(t(ip_data), 5, nstart = 10, iter.max = 50)$cluster
  params  <- updt_params(K, alp, ip_data, dm, niw_p, z)

  # Number of Rejections per Iteration
  rej_samp <- rep(list(list()), K)
  params$rej_samp <- rej_samp

  # Store result
  rslt    <- rep(list(params), (numIter-burnIn))

  brks  <- c(seq(-4, 0, 0.05), seq(1, 4, 0.05))
  count_bins_sum <- rej_summary(-5, brks)

  # *********************** Gibbs Sampling loop ************************* #
  for(iter in 1:numIter) {
    if(iter%%100 == 0)
      print(paste("Iter", iter))

    # Reorder observations
    new_order <- sample(ncol(ip_data), replace = F)
    ip_data   <- matrix(ip_data[, new_order], nrow = dm)
    z <- z[new_order]

    if (mix == "motg" | mix == "motgt"){
      # Update parameters for each cluster by first sampling the rejections
      for(i in 1:K) {
        pm <- updt_params_component(X = ip_data[, z == i],
                                    mu = params$mu[[i]],
                                    cv = params$cv[[i]],
                                    niw_p, constr, thr, dm)

        # Store parameters
        params$mu[[i]] <- pm$mu
        params$cv[[i]] <- pm$cv

        # Count number of rejections in partition
        if (mix == "motgt"){
          if (dm == 1){
            if (iter == numIter)
              rej_samp[[i]] <- pm$rej_samp
          } else {

            # For now, save rejections
            params$rej_samp[[i]] <- pm$rej_samp
          }
        } else {
          params$rej_samp[[i]] <- pm$rej_samp
        }

        params$cnt_rej[i]  <- pm$cnt_rej
        params$no[i]       <- pm$no

      }

      # Update assignments
      if (mix == "motgt"){
        z  <- updt_assgn_motgt(K, ip_data, N, dm, params, lower, upper)
      } else {
        z  <- updt_assgn_motg(K, ip_data, z, dm, rej = params$rej_samp, N, params, thr)
      }

      # Store new cluster parameters
      params$z <- z

      # Update weights
      new_wt <- updt_wghts(z, K, alp)
      params$wt <- new_wt

      if (iter > burnIn)
        rslt[[iter-burnIn]] <- params

      # Set rejection storage to empty again
      rej_samp <- rep(list(list()), K)

    } else if (mix == "tmog"){
      # Sample rejections
      smpls_rejs <- impute_rej(K, params, N, constr, thr, dm, ceil_rej)

      # Function to update params
      params <- updt_params_tmog(K, alp,
                                  x_rejs = cbind(ip_data, smpls_rejs$rejs),
                                  dm, params$wt,
                                  niw_p, z_rejs = c(z, smpls_rejs$z))
      # Update weights
      new_wt <- updt_wghts(c(z, smpls_rejs$z), K, alp)

      params$wt <- new_wt

      # Count number of rejections in partition
      if (dm == 1){
        # only save rejs summary for one dimensional data
        if (thr != 0){
          count_bins <- rej_summary(smpls_rejs$rejs, brks)
          count_bins_sum <- count_bins_sum + count_bins
        } else {
          count_bins_sum <- NULL
        }
      } else {
        # For now, save rejections
        if (iter > burnIn)
          rslt[[iter-burnIn]]$rej_smpls <- smpls_rejs$rejs
      }

      # Update assignments
      z <- updt_assgn_tmog(K, ip_data, N, dm, params)

      if (iter > burnIn){
        rslt[[iter-burnIn]] <- params
        rslt[[iter-burnIn]]$cnt_rej  <- smpls_rejs$rej_acc
      }
    } else {
      stop("Mixture is not one of motg, motgt, or tmog")
    }

  }

  return(list(res = rslt, rej_samp = rej_samp))

}



