#' Transform a mortality rate model to an improvement rate model
#'
#' Transforms a StMoMo object representing a mortality rate model
#' into an iMoMo object representing an improvement rate model
#'
#' @param model A StMoMo object representing the mortality rate model
#'   to be transformed into an improvement rate model. This StMoMo
#'   model must include an static age function \eqn{\alpha_x} and must
#'   use a \code{"log"} link.
#'
#' @param type defines the type of estimation method to be used.
#'   \code{"indirect"} would estimate the equivalent mortality rate model
#'   and the transform the model into improvement rates. \code{"direct"}
#'   would estimate the model directly on the improvement rate model.
#'
#' @param  constFun function defining the identifiability constraints of the
#'   model. It must be a function of the form
#'   \code{constFun <- function(ax, bx, kt, b0x, gc, wxt, ages)} taking a set
#'   of fitted model parameters and returning a list
#'   \code{list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
#'   of the model parameters with the identifiability constraints applied. If
#'   omitted no identifiability constraints are applied to the model.
#'
#' @export
StMoMo2iMoMo <- function(model, type = c("indirect", "direct"),
                         constFun = function(ax, bx, kt, b0x, gc, wxt, ages)
                           list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc) ){

  #---------------------------------------------------------------------
  # Check inputs
  #---------------------------------------------------------------------

  #StMoMo model
  type <- match.arg(type)
  if (sum("StMoMo" %in% class(model)) == 0) {
    stop("The model argument needs to be of class StMoMo")
  }
  if (model$link != "log")
    stop("The base StMoMo model needs to use a log link")
  if (type == "fitted" && model$staticAgeFun == FALSE) {
    stop("For models of type fitted the StMoMo model needs
         to include an static age function")
  }

  #Estimation model
  if (type == "direct"){
    model <- StMoMo(link = model$link, staticAgeFun = FALSE,
                    periodAgeFun = model$periodAgeFun,
                    cohortAgeFun = model$cohortAgeFun)
  }

  #Structure of the model is the same but without a static age
  textFormula <- sub("log m\\[x,t\\] = a\\[x\\] \\+ ", "eta\\[x,t\\] = ",
                     model$textFormula)
  out <- list(staticAgeFun = FALSE,
              periodAgeFun = model$periodAgeFun,
              cohortAgeFun = model$cohortAgeFun,
              N = model$N,
              textFormula = textFormula,
              type = type,
              constFun = constFun,
              model = model)

  #Estimation type
  if (type == "indirect")
    class(out) <-  c("iMoMo", "iMoMoI", "StMoMo")
  else
    class(out) <-  c("iMoMo", "iMoMoD", "StMoMo")
  out

}
