% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimpleVAR2.r
\name{simulate.simpleVAR2}
\alias{simulate.simpleVAR2}
\title{Simulate a VAR model}
\usage{
\method{simulate}{simpleVAR2}(object, nsim = 10, seed = NULL, colnames = NULL)
}
\arguments{
\item{nsim}{number of periods for the simulated series.}

\item{seed}{either \code{NULL} or an integer that will be used in a
call to \code{\link{set.seed}} before simulating the time series.
The default, \code{NULL} will not change the random generator state.}

\item{colnames}{Name of the columns (time periods)}

\item{object:}{An object of class \code{"simpleVAR2"}.}
}
\description{
Returns one simulated path of the VAR model in \code{object}.
}
