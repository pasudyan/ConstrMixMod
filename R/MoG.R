#' Mixture of Gaussian function
#' 
#' This function generates a mixture of Gaussians samples. 
#' @param K number of clusters (default 2)
#' @param params list with names: wt = weights, mu = mean, and cv = covariance
#' @param N number of random samples (default 100)
#' @param dm dimension
#' @return list of the vector of random samples and the associated cluster assignment
#' @keywords Mixture of Gaussians
#' @export
#' @examples 
#' params = list(wt = c(0.5, 0.5), mu=c(0.1, 0.5), cv = c(0.1, 0.7))
#' gen_mog(K = 2, params, N = 100, dm = 2)

gen_mog <- function(K = 2, params,  N = 100, dm) {
  
  wt  <- params$wt
  mu  <- params$mu
  cv  <- params$cv

  X <- matrix(rep(0, N*dm), ncol = dm)
  z <- sample(K, N, replace = T, prob = wt)

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
  
  list(x = t(X), z = z)
}

