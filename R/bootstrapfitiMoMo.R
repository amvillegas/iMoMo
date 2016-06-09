#' Bootstrap an Improvement Rate Mortality Model of type fitted
#' @export
bootstrap.fitiMoMoF <- function(object, nBoot = 1,
                                type = c("semiparametric", "residual"),
                                deathType = c("observed", "fitted"), ...) {


  out <- bootstrap(object$fittingModel, nBoot = nBoot, type = type,
                   deathType = deathType)

  out$model <- object
  out$bootParameters <- lapply(out$bootParameters, transParamiMoMoF)
  class(out) <- c("bootiMoMoF", "bootiMoMo", class(out))
  out
}

#' Bootstrap an Improvement Rate Mortality Model of type fitted
#' with constant improvement rates
#' @export
bootstrap.fitiMoMoCIF <- function(object, nBoot = 1,
                                  type = c("semiparametric", "residual"),
                                  deathType = c("observed", "fitted"), ...) {
  out <- bootstrap(object$fittingModel, nBoot = nBoot, type = type,
                   deathType = deathType)

  out$model <- object
  out$bootParameters <- lapply(out$bootParameters, transParamiMoMoCIF)
  class(out) <- c("bootiMoMoCIF", "bootiMoMo", class(out))
  out
}


# #' Bootstrap an Improvement Rate Mortality Model of type observed
# #' @export
# bootstrap.fitiMoMoO <- function(object, nBoot = 1,
#                                 type = c("semiparametric", "residual"),
#                                 deathType = c("observed", "fitted"), ...) {
#
#
#   out <- bootstrap(object$fittingModel, nBoot = nBoot, type = type,
#                    deathType = deathType)
#
#   out$model <- object
#   out$bootParameters <- lapply(out$bootParameters, transParamiMoMoO)
#   class(out) <- c("bootiMoMoO", "bootiMoMo", class(out))
#   out
# }



#' Bootstrap an Improvement Rate Mortality Model of type observed
#' @export
bootstrap.fitiMoMoO <- function(object, nBoot = 1,
                                type = c("semiparametric", "residual"),
                                deathType = c("observed", "fitted"), ...) {

  type <- match.arg(type)
  deathType <- match.arg(deathType)
  link <- object$model$link

  if (type == "residual") {
    bootSamples <- StMoMo:::genPoissonResBootSamples(
      devRes = residuals(object$fittingModel, scale = FALSE),
      dhat = fitted(object$fittingModel, type = "deaths"), nBoot)
    bootSamples <- lapply(bootSamples, function(x) cbind(object$Dxt[,1], x))

  } else if (type == "semiparametric") {
    if (deathType == "observed") {
      D <- object$Dxt
    } else {
      D <- fitted(object$fittingModel, type = "deaths")
      D <- cbind(object$Dxt[,1], D)
    }
    bootSamples <- StMoMo:::genPoissonSemiparametricBootSamples(D = D, nBoot)

  }
  #Fit the model to each of
  refit <- function(Dxt) {
    getMinimalFitiMoMo(fit(object = object$model, Dxt = Dxt, Ext = object$Ext,
                           ages = object$ages,
                           years = c(object$years[1]-1, object$years),
                           wxt = object$wxt, start.ax = object$ax,
                           start.bx = object$bx, start.kt = object$kt,
                           start.b0x = object$b0x, start.gc = object$gc,
                           verbose = FALSE))
  }
  bootParameters <- lapply(bootSamples, FUN = refit)
  structure(list(bootParameters = bootParameters, model = object, type = type,
                 deathType = deathType, call = match.call()),
            class = c("bootiMoMoO", "bootiMoMo","bootStMoMo"))


}

#' Extract a lighter version of a  Improvment Rate Mortality Model
#' @keywords internal
getMinimalFitiMoMo <- function(object) {
  structure(list(model = object$model, ax = object$ax, bx = object$bx,
                 kt = object$kt, b0x = object$b0x, gc = object$gc,
                 oxt = object$oxt, ages = object$ages, years = object$years,
                 cohorts = object$cohorts), class = class(object))
}

