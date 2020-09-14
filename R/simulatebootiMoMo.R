#' Simulate future sample paths from a Bootstrapped Improvement
#' Rate Model
#'
#' Simulate mortaility improvment rates and mortality rates using
#' the a bootstrapped a mortality improvement rate model.
#' The period indexes are \eqn{\kappa_t^{(i)}, i = 1,..N,} are forecasted
#' using integrated vector autoregressive model. The cohort index
#' \eqn{\gamma_{t-x}} is forecasted using an ARIMA\eqn{(p, d, q)}.
#' By default an ARIMA\eqn{(1, 1, 0)} with a constant is used.
#'
#' @param object an object of class \code{"bootiMoMo"} with the bootstrapped
#' parameters of an improvement rate model.
#' @param nsim number of sample paths to simulate from each bootstrapped
#' sample. Thus if there are \code{nBoot} bootstrapped samples the total
#' number of paths will be \code{nsim * nBoot}.
#' @inheritParams simulate.simpleVAR2
#' @inheritParams forecast.fitiMoMo
#'
#' @return A list of class \code{"simiMoMo"} with components:
#'
#' \item{improvements}{a three dimensional array with the future simulated improvement rates.}
#' \item{rates}{a three dimensional array with the future simulated mortality rates.}
#' \item{ages}{vector of ages corresponding to the first dimension of \code{improvements}.}
#' \item{years}{vector of years for which a forecast has been produced. This
#'  corresponds to the second dimension of \code{improvements}.}#'
#' \item{kt.s}{ information on the simulated paths of the period indexes of the model.
#' This is a list with the \code{model} fitted to \eqn{\kappa_t}; the simulated paths
#' \code{sim}; and the \code{years} for which a forecast was produced. If the
#' model does not have any age-period terms (i.e. \eqn{N=0}) this is set to
#' \code{NULL}.}
#' \item{gc.s}{ information on the simulated paths of the cohort index of
#' the model. This is a list with the \code{model} fitted to \eqn{\gamma_c};
#' the simulated paths (\code{sim}); and the \code{cohorts} for which
#' simulations were produced. If the mortality model does not have a cohort
#' effect this is set to \code{NULL}.}
#' \item{fittedImprovements}{ a three dimensional array with the in-sample
#' improvements of the model for the years for which the improvement rate model
#' was fitted.}
#'
#' @details
#' For further details see \code{\link{simulate.fitiMoMo}}.
#'
#' @seealso \code{\link{bootstrap.fitiMoMo}}, \code{\link{simulate.fitiMoMo}}
#'
#' @export
simulate.bootiMoMo <-function(object, nsim = 1, seed = NULL, h = 50,
                              kt.order = c(1, 0, 0),
                              kt.include.constant = TRUE,
                              kt.include.trend = FALSE,
                              gc.order = c(1, 0, 0),
                               gc.include.constant = TRUE,
                              jumpRates = NULL,
                               kt.lookback = NULL, gc.lookback = NULL,
                               ...) {

  ages <- object$model$ages
  nAges <- length(ages)
  nBoot <- length(object$bootParameters)
  nPath <- nsim * nBoot

  ## Handle generator seed
  if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)
  if (is.null(seed))
    RNGstate <- .Random.seed
  else {
    R.seed <- .Random.seed
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  ## Do simulations

  #Initialse output variables
  modelFit <- object$bootParameters[[1]]
  modelFit$Dxt <- object$model$Dxt
  modelFit$Ext <- object$model$Ext
  tempSim <- simulate.fitiMoMo(object = modelFit,
                                nsim = nsim, seed = NULL, h = h,
                               kt.order = kt.order,
                               kt.include.constant = kt.include.constant,
                               kt.include.trend = kt.include.trend,
                               gc.order = gc.order,
                                gc.include.constant = gc.include.constant,
                               jumpRates = jumpRates,
                                kt.lookback = kt.lookback,
                                gc.lookback = gc.lookback)
  improvements <- array(NA, c(dim(tempSim$improvements)[1:2], nPath),
                 list(dimnames(tempSim$improvements)[[1]],
                      dimnames(tempSim$improvements)[[2]], 1:nPath))
  rates <- array(NA, c(dim(tempSim$rates)[1:2], nPath),
                 list(dimnames(tempSim$rates)[[1]],
                      dimnames(tempSim$rates)[[2]], 1:nPath))
  fittedImprovements <- array(NA, c(dim(tempSim$fittedImprovements)[1:2], nPath),
                  list(dimnames(tempSim$fittedImprovements)[[1]],
                       dimnames(tempSim$fittedImprovements)[[2]], 1:nPath))
  if (is.null(tempSim$kt.s)) {
    kt.s <- NULL
  } else {
    kt.s <- list(sim = NULL, years = tempSim$kt.s$years)
    kt.s$sim <- array(NA, c(dim(tempSim$kt.s$sim)[1:2], nPath),
                      list(dimnames(tempSim$kt.s$sim)[[1]],
                           dimnames(tempSim$kt.s$sim)[[2]], 1:nPath))
  }
  if (is.null(tempSim$gc.s)) {
    gc.s <- NULL
  } else {
    gc.s <- list(sim = NULL, cohorts = tempSim$gc.s$cohorts)
    gc.s$sim <- array(NA, c(dim(tempSim$gc.s$sim)[1], nPath),
                      list(dimnames(tempSim$gc.s$sim)[[1]], 1:nPath))
  }

  # Do simulations for each bootstrap sample
  i <- 1
  while (i <= nBoot) {
    improvements[, , ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$improvements
    rates[, , ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$rates
    fittedImprovements[, , ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$fittedImprovements
    if (!is.null(kt.s))
      kt.s$sim[, , ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$kt.s$sim
    if (!is.null(gc.s))
      gc.s$sim[, ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$gc.s$sim
    i <- i + 1
    if ( i <= nBoot) {
      modelFit <- object$bootParameters[[i]]
      modelFit$Dxt <- object$model$Dxt
      modelFit$Ext <- object$model$Ext
      tempSim <- simulate.fitiMoMo(object = modelFit,
                                   nsim = nsim, seed = NULL, h = h,
                                   kt.order = kt.order,
                                   kt.include.constant = kt.include.constant,
                                   kt.include.trend = kt.include.trend,
                                   gc.order = gc.order,
                                   gc.include.constant = gc.include.constant,
                                   jumpRates = jumpRates,
                                   kt.lookback = kt.lookback,
                                   gc.lookback = gc.lookback)

    }
  }
  structure(list(improvements = improvements,
                 rates = rates, ages = tempSim$ages,
                 years = tempSim$years, kt.s = kt.s, gc.s = gc.s,
                 fittedImprovements = fittedImprovements,
                 model = object, call = match.call()),
            class = "simiMoMo")
}

