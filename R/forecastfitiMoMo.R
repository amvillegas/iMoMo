#' Forecast mortality rates using an improvement rate Model
#'
#' @export
forecast.fitiMoMo <-function(object, h = 50,
                             kt.order = c(1, 0, 0),
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
  # forcastRates <- rates[, (nYears + 1):(nYears + h)]
  # fittedRates <- rates[, 1:nYears]
  # if (jumpchoice == "actual") {
  #   jumpoffRates <- (object$Dxt / object$Ext)[, nYears]
  #   forcastRates <- forcastRates * jumpoffRates / fittedRates[ , nYears]
  # }
  #
  #prepare output
  structure(list(improvements = forcastImprovements, rates = rates,
                 ages = agesFor, years = yearsFor,
                 kt.f = kt.f, gc.f = gc.f,
                 fittedImprovements = fittedImprovements, model = object,
                 call = match.call()),
            class = "foriMoMo")
}
