#' Forecast mortality rates using an improvement rate model
#'
#' Forecast mortaility improvment rates and mortality rates using
#' the fit from a mortality improvement rate model.
#' The period indexes are \eqn{\kappa_t^{(i)}, i = 1,..N,} are forecasted
#' using integrated vector autoregressive model. The cohort index
#' \eqn{\gamma_{t-x}} is forecasted using an ARIMA\eqn{(p, d, q)}.
#' By default an ARIMA\eqn{(1, 1, 0)} with a constant is used.
#'
#' @param object an object of class \code{"fitiMoMo"} with the fitted
#' parameters of an improvement rate mortality model.
#' @param h number of years ahead to forecast.
#' @param kt.order  an optional vector indicating the order of
#' autorregression and of differetiation of the VARI model. The two
#' components \eqn{(p, d)} are the AR order and the degree of  differencing.
#' @param kt.include.constant a logical value indicating if the VARI model
#' should include a constant value. The default is \code{TRUE}.
#' @param kt.include.trend a logical value indicating if the VARI model
#' should have a linear trend. The default is \code{FALSE}.
#' @param gc.order a specification of the ARIMA model for the cohort effect:
#' the three components \eqn{(p, d, q)} are the AR order, the degree of
#' differencing, and the MA. The default is an ARIMA\eqn{(1, 0, 0)}.
#' @param gc.include.constant a logical value indicating if the ARIMA model
#' should include a constant value. The default is \code{TRUE}.
#' @param jumpRates optional vector of moratlity rates for the last year used as starting rates.
#' for the projection. If it is not provided the rates from the  the actual rates from
#' the final year are used.
#' @param kt.lookback optional argument to specify the look-back window to use
#' in the estimation of the time series model for the period indexes. By
#' default all the estimated values are used. If
#' \code{kt.lookback} is provided then the last \code{kt.lookback}
#' years of \eqn{\kappa_t^{(i)}, i = 1,..N,} are used.
#' @param gc.lookback optional argument to specify the look-back window to use
#'        in the estimation of the ARIMA model for the cohort effect. By
#'        default all the estimated values are used in estimating the ARIMA
#'        model. If \code{gc.lookback} is provided then the last
#'        \code{gc.lookback} years of \eqn{\gamma_{t-x}} are used.
#' @param ... other arguments.
#'
#' @return A list of class \code{"foriMoMo"} with components:
#'
#' \item{improvements}{a matrix with the point forecast of the improvement rates.}
#' \item{rates}{a matrix with the point forecast of the rates.}
#' \item{ages}{ vector of ages corresponding to the rows of \code{improvements}.}
#' \item{years}{vector of years for which a forecast has been produced. This
#'  corresponds to the columns of \code{improvements}.}#'
#' \item{kt.f}{ forecasts of period indexes of the model. This is a list with
#' the \code{model} fitted to \eqn{\kappa_t}; the \code{mean}(central)
#' forecast; and the \code{years} for which a forecast was produced. If the
#' model does not have any age-period terms (i.e. \eqn{N=0}) this is set to
#' \code{NULL}.}#'
#' \item{gc.f}{ forecasts of cohort index of the model. This is a list with
#' the \code{model} fitted to \eqn{\gamma_c}; the \code{mean}(point) forecast;
#' and the \code{cohorts} for which a forecast was produced. If the model
#' does not have a cohort effect this is set to \code{NULL}.} #'
#' \item{fittedImprovements}{ a matrix with the fitted in-sample improvements of
#' the model for  the years for which the improvement rate model was fitted.}#'
#' \item{model}{the model fit from which the forecast was produced.}
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
#' Fitting and forecasting of the VAR model is done using the fucntion
#' \code{\link{simpleVAR2}}.
#'
#' Fitting and forecasting of the ARIMA model for the cohort index
#' is done with function \code{\link[forecast]{Arima}} from package
#' \pkg{forecast}. See the latter function for further details on
#' input arguments \code{gc.order} and \code{gc.include.constant}.
#'
#' @export
forecast.fitiMoMo <-function(object, h = 50,
                             kt.order = c(1, 0),
                             kt.include.constant = TRUE,
                             kt.include.trend = FALSE,
                              gc.order = c(1, 0, 0),
                              gc.include.constant = TRUE,
                             jumpRates = NULL,
                              kt.lookback = NULL, gc.lookback = NULL,
                              ...) {

  #forecast kt
  kt <- object$kt
  years <- object$years
  nYears <- length(years)
  if (is.null(kt.lookback)) kt.lookback <- nYears
  if (kt.lookback <= 0)
    stop("kt.lookback must be positive")
  kt.lookback <- min(c(kt.lookback, nYears))
  yearsFor <- (years[nYears] + 1):(years[nYears] + h)
  agesFor <- object$ages
  nAges <- length(object$ages)
  kt.h <- kt
  kt.f <- NULL
  kt.model <- NULL
  years.h <- years
  years.f <- yearsFor
  if (object$model$N > 0) {
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
    kt.for <- forecast(kt.model, h = h + kt.hNA)
    if (kt.hNA > 0) {
      years.h <- years[-((kt.nNA+1):nYears)]
      years.f <- c(years[(kt.nNA+1):nYears], years.f)
      kt.h <-array(kt.h[, 1:kt.nNA], c(nrow(kt), kt.nNA))
      dimnames(kt.h)[[2]] <- years.h
    }
    kt.f <- list(mean = t(kt.for), model = kt.model, years = years.f)
  }
  #forecast gc
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
  gc.f <- NULL
  cohorts.f <- (cohorts[nCohorts] + 1):(cohorts[nCohorts] + h)
  if (!is.null(object$model$cohortAgeFun)) {
    gc.nNA <- max(which(!is.na(gc)))
    gc.hNA <- nCohorts - gc.nNA
    gc.model <- forecast::Arima(gc[(1 + nCohorts - gc.lookback):gc.nNA],
                                order = gc.order,
                                include.constant = gc.include.constant)
    gc.for <- forecast(gc.model, h = h + gc.hNA)

    if (gc.hNA > 0) {
      gc.h <- gc[-((gc.nNA+1):nCohorts)]
      cohorts.h <- cohorts[-((gc.nNA + 1):nCohorts)]
      cohorts.f <- c(cohorts[(gc.nNA + 1):nCohorts], cohorts.f)
    }
    gc.f <- list(mean = as.vector(gc.for$mean),
                 model = gc.model, cohorts = cohorts.f)

    names(gc.f$mean)  <- cohorts.f
  }
  #predict rates
  improvements <- StMoMo:::predictLink(ax = object$ax, bx = object$bx,
                                kt = cbind(kt.h, kt.f$mean),
                                b0x = object$b0x, gc = c(gc.h, gc.f$mean),
                                oxt = NULL, ages = object$ages,
                                years = c(years.h, years.f))

  forcastImprovements <- improvements[, (nYears + 1):(nYears + h)]
  fittedImprovements <- improvements[, 1:nYears]

  #Apply jump-off option
  RF <- exp(-t(apply(forcastImprovements, 1, cumsum)))
  if (is.null(jumpRates)) {
      jumpRates <- object$Dxt[,nYears + 1] / object$Ext[,nYears + 1]
  }
  rates <- array(jumpRates, dim = dim(RF)) * RF

  #prepare output
  structure(list(improvements = forcastImprovements, rates = rates,
                 ages = agesFor, years = yearsFor,
                 kt.f = kt.f, gc.f = gc.f,
                 fittedImprovements = fittedImprovements, model = object,
                 call = match.call()),
            class = "foriMoMo")
}
