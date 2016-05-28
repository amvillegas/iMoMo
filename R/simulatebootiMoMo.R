
#' Simulate future sample paths from a Bootstrapped Improvement
#' Rate Model
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

