#' Simulate future sample paths from an improvement  Model
#'
#' Simulate mortaility improvment rates and mortality rates using
#' the fit from a mortality improvement rate model.
#' The period indexes are \eqn{\kappa_t^{(i)}, i = 1,..N,} are forecasted
#' using integrated vector autoregressive model. The cohort index
#' \eqn{\gamma_{t-x}} is forecasted using an ARIMA\eqn{(p, d, q)}.
#' By default an ARIMA\eqn{(1, 1, 0)} with a constant is used.
#'
#' @param object an object of class \code{"fitiMoMo"} with the fitted
#' parameters of an improvement rate mortality model.
#' @param nsim number of sample paths to simulate.
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
#'
#' The modelling of the period indexes \eqn{kappa_t} is done using a integrated vector
#' autoregressive model of differencing order \eqn{d} and autorregressive order \eqn{p}:
#'
#' \deqn{\Delta^d k_t = C+Dt+\sum_{i=1}^p A_i \Delta^d k_{t-i} + \epsilon_t}
#'
#' where \eqn{C} and  \eqn{D} are \eqn{N}-dimensional vectors for parameters and
#' \eqn{A_1,...,A_p} are \eqn{N\times N} matrices of autoregressive parampeters.
#' If \code{kt.include.constant} is \code{TRUE} then \eqn{C} is included in the
#' equation. Similarly, if \code{kt.include.trend} is \code{TRUE} then \eqn{D}
#' is included in the equation.
#'
#' Fitting and simulating of the VAR model is done using the fucntion
#' \code{\link{simpleVAR2}}.
#'
#' Fitting and simulating of the ARIMA model for the cohort index
#' is done with function \code{\link[forecast]{Arima}} from package
#' \pkg{forecast}. See the latter function for further details on
#' input arguments \code{gc.order} and \code{gc.include.constant}.
#'
#'@export
simulate.fitiMoMo <- function(object, nsim = 1000, seed = NULL, h = 50,
                               kt.order = c(1, 0),
                               kt.include.constant = TRUE,
                               kt.include.trend = FALSE,
                               gc.order = c(1, 0, 0),
                               gc.include.constant = TRUE,
                               jumpRates = NULL,
                               kt.lookback = NULL, gc.lookback = NULL,
                               ...) {

  #Handle generato seed
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

  #Fit model to kt
  N <- object$model$N
  kt <- object$kt
  years <- object$years
  nYears <- length(years)
  if (is.null(kt.lookback))
    kt.lookback <- nYears
  if (kt.lookback <= 0)
    stop("kt.lookback must be positive")
  kt.lookback <- min(c(kt.lookback, nYears))
  yearsSim <- (years[nYears] + 1):(years[nYears] + h)
  ages <- object$ages
  nAges <- length(ages)
  kt.h <- kt
  kt.path <- NULL
  kt.model <- NULL
  years.h <- years
  years.s <- yearsSim
  if (N > 0) {
    kt.nNA <- max(which(!is.na(kt[1, ])))
    kt.hNA <- nYears - kt.nNA
    if (kt.include.constant & kt.include.trend) type <- "both"
    if (kt.include.constant & !kt.include.trend) type <- "const"
    if (!kt.include.constant & kt.include.trend) type <- "trend"
    if (!kt.include.constant & !kt.include.trend) type <- "none"
    if (object$model$N > 1) {
      kt.obs <- t(kt[, (1 + nYears - kt.lookback):kt.nNA])
    } else {
      kt.obs <- (kt[, (1 + nYears - kt.lookback):kt.nNA])
    }
    kt.model <- simpleVAR2(kt.obs,
                           p = kt.order[1], d = kt.order[2],
                           type = type)
    if (kt.hNA > 0) {
      years.h <- years[-((kt.nNA+1):nYears)]
      years.s <- c(years[(kt.nNA+1):nYears], years.s)
      kt.h <- array(kt.h[, 1:kt.nNA], c(nrow(kt), kt.nNA))
      dimnames(kt.h)[[2]] <- years.h
    }
    kt.sim <- array(NA,c(N, length(years.s), nsim), list(1:N, years.s, 1:nsim))
  }

  #fit model to gc
  gc <- object$gc
  cohorts <- object$cohorts
  nCohorts <- length(cohorts)
  if (is.null(gc.lookback)) gc.lookback <- nCohorts
  if (gc.lookback <= 0)
    stop("gc.lookback must be positive")
  gc.lookback <- min(c(gc.lookback, nCohorts))
  gc.h <- gc
  cohorts.h <- cohorts
  gc.model <- NULL
  gc.path <- NULL
  cohorts.s <- (cohorts[nCohorts] + 1):(cohorts[nCohorts] + h)
  if (!is.null(object$model$cohortAgeFun)) {
    gc.nNA <- max(which(!is.na(gc)))
    gc.hNA <- nCohorts - gc.nNA
    gc.model <- forecast::Arima(gc[(1 + nCohorts - gc.lookback):gc.nNA],
                                order = gc.order,
                                include.constant = gc.include.constant)
    if (gc.hNA > 0) {
      gc.h <- gc[-((gc.nNA+1):nCohorts)]
      cohorts.h <- cohorts[-((gc.nNA + 1):nCohorts)]
      cohorts.s <- c(cohorts[(gc.nNA + 1):nCohorts], cohorts.s)
    }
    gc.sim <- array(NA, c(length(cohorts.s), nsim), list(cohorts.s, 1:nsim))
  }

  #Do simulations
  forcastImprovements <- array(NA, c(nAges, h, nsim),
                        list(ages, yearsSim, 1:nsim))
  fittedImprovements  <- array(NA, c(nAges, nYears, nsim),
                        list(ages, years, 1:nsim))
  rates <- forcastImprovements
  if (is.null(jumpRates)) {
    jumpRates <- object$Dxt[,nYears + 1] / object$Ext[,nYears + 1]
  }

  for (i in 1:nsim) {
    if (!is.null(kt.model)) {
      kt.path <- t(simulate(kt.model,h + kt.hNA))
      kt.sim[, , i] <- kt.path
    }
    if (!is.null(gc.model)) {
      gc.path <- as.vector(simulate(gc.model, h + gc.hNA))
      gc.sim[, i] <- gc.path
    }
    improvementsi <- StMoMo:::predictLink(ax = object$ax, bx = object$bx,
                                         kt = cbind(kt.h, kt.path),
                                         b0x = object$b0x, gc = c(gc.h, gc.path),
                                         oxt = NULL, ages = object$ages,
                                         years = c(years.h, years.s))

    forcastImprovements[, , i] <- improvementsi[, (nYears + 1):(nYears + h)]
    fittedImprovements[, , i] <- improvementsi[, 1:nYears]
    RF <- exp(-t(apply(forcastImprovements[, , i], 1, cumsum)))
    rates[, , i] <- array(jumpRates, dim = dim(RF)) * RF
  }

  if (is.null(kt.model)) {
    kt.s <- NULL
  } else {
    kt.s <- list(sim = kt.sim, model = kt.model, years = years.s)
  }
  if (is.null(gc.model)) {
    gc.s <- NULL
  } else {
    gc.s <- list(sim = gc.sim, model = gc.model, cohorts = cohorts.s)
  }
  structure(list(improvements = forcastImprovements,  rates = rates, ages = ages,
                 years = yearsSim, kt.s = kt.s, gc.s = gc.s,
                 fittedImprovements = fittedImprovements,
                 model = object, call = match.call()),
            class ="simiMoMo")

}
