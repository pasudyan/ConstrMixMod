% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updateFuncs.R
\name{updateClustAssignTmog}
\alias{updateClustAssignTmog}
\title{Update assignment for univariate TMoG model}
\usage{
updateClustAssignTmog(K, ipdata, N, params)
}
\arguments{
\item{K}{max number of clusters}

\item{ipdata}{data matrix of size p-by-N where N are number of observations and p are number of dimensions}

\item{N}{total number of observations}

\item{params}{list of parameters which include weights (wt), mean (mu), and covariance matrix (cv) for all clusters}
}
\value{
vector of new cluster assignments
}
\description{
This function updates the cluster assignment of each observation in the univariate TMoG model
}
