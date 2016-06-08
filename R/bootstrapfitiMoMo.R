#' Bootstrap an Improvement Rate Mortality Model of type observed
#' @export
bootstrap.fitiMoMoO <- function(object, nBoot = 1,
                                type = c("semiparametric", "residual"),
                                deathType = c("observed", "fitted"), ...) {


out <- bootstrap(object$fittingModel, nBoot = nBoot, type = type,
                 deathType = deathType)

out$model <- object
out$bootParameters <- lapply(out$bootParameters, transParamiMoMoO)
class(out) <- c("bootiMoMoO", "bootiMoMo", class(out))
out
}


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