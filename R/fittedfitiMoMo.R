#' Compute fitted values for an Improvement Rate Model
#'
#' Returns fitted values for the data used in fitting the
#' Improvement Rate Stochastic  Mortality Model.
#'
#' @param object an object of class \code{"fitiMoMo"} with the fitted
#' parameters of an improvment rate  model.
#' @param type the type of the fitted values that should be returned. The
#' alternatives are \code{"improvementes"}(default), \code{"rates"}, and
#' \code{"deaths"}.
#'
#' @return A matrix with the fitted values.
#'
#' @export
fitted.fitiMoMo <- function(object, type = c("improvements", "rates", "deaths"), ...) {

  type <- match.arg(type)
  improvements <- with(object, StMoMo:::predictLink(ax = ax, bx = bx, kt = kt, b0x = b0x,
                                           gc = gc, oxt = NULL, ages = ages,
                                           years = years))
  n <- dim(object$Dxt)[2]
  if (object$model$type == "indirect"){
    RF <- exp(-t(apply(cbind(0, improvements), 1, cumsum)))
    rates <- array(exp(object$Ax), dim = dim(RF)) * RF
    dimnames(rates) <- dimnames(object$Ext)
    deaths <- object$Ext * rates
  } else {
    rates <- exp(-improvements) * object$Dxt[, 1:(n - 1)] / object$Ext[, 1:(n - 1)]
    deaths <- object$Ext[, 1:(n - 1)] * rates
  }

  switch(type, rates = rates, deaths = deaths, improvements = improvements)
}
