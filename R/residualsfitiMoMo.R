#' Extract deviance residuals of an Improvement Rate  Model
#'
#' Compute deviance residuals of a fitted improvement rate  Model.
#' These residuals can be plotted using \code{\link[StMoMo]{plot.resStMoMo}}.
#'
#' @param object an object of class \code{"fitiMoMo"} with the fitted
#' parameters of an improvment mortality model.
#' @param scale logical indicating whether the residuals should be scaled
#' or not by dividing the deviance by the  overdispersion of the model.
#' Default is \code{TRUE}.
#' @param ... other arguments.
#'
#' @return An object of class \code{"resStMoMo"} with the residuals. This
#' object has components:
#'   \item{residuals}{ a matrix with the residuals.}
#'   \item{ages}{ ages corresponding to the rows in \code{residuals}.}
#'   \item{years}{ years corresponding to the columns in \code{residuals}.}
#'
#' @seealso \code{\link[StMoMo]{plot.resStMoMo}}
#'
#' @export
residuals.fitiMoMoI <- function(object, scale = TRUE, ...) {
  D <- object$Dxt + 0.00001 #Add a small amount to compensate for the
  #possibility of 0 deaths
  W <- object$wxt
  ind <- (W > 0)
  res <- array(NA, dim(W), dimnames = dimnames(D))
  Dhat <- fitted(object, type = "deaths")
  res[ind] <- 2 * W[ind] * (D[ind] * log(D[ind] / Dhat[ind]) - (D[ind] - Dhat[ind]))
  signRes <- sign(D - Dhat)

  if (scale)
    phi <- sum(res[ind]) / (object$nobs - object$npar)
  else
    phi <- 1
  res <- signRes * sqrt(abs(res) / phi)
    structure(list(residuals = res, ages = object$ages,
                   years = c(object$years[1]-1, object$years)),
            class = "resStMoMo")
}
