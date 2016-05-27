#' Bootstrap an Improvement Rate Mortality Model of type observed
#' @export
bootstrap.fitiMoMoO <- function(object, nBoot = 1,
                                type = c("semiparametric", "residual"),
                                deathType = c("observed", "fitted"), ...) {


out <- bootstrap(object$fittingModel, nBoot = nBoot, type = type,
                 deathType = deathType)

out$model <- object
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
  class(out) <- c("bootiMoMoF", "bootiMoMo", class(out))
  out
}
