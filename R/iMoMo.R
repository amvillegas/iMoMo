#' Create a new Improvement Rate Mortality Model
#'
#' Initialises an iMoMo oject which represents a generalised
#' Age-Period-Cohort improvement rate mortality model
#'
#' @param staticAgeFun logical value indicating if a static age function
#'   \eqn{\alpha_x} is to be included.
#'
#' @param periodAgeFun  a list of length \eqn{N} with the definitions of the
#'   period age modulating parameters \eqn{\beta_x^{(i)}}. Each entry can take
#'   values: \code{"NP"} for non-parametric age terms, \code{"1"} for
#'   \eqn{\beta_x^{(i)}=1} or a predefined parametric function of
#'   age (see details). Set this to \code{NULL} if there are no period terms
#'   in the model.
#'
#' @param cohortAgeFun defines the cohort age modulating parameter
#'   \eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for non-parametric
#'   age terms, \code{"1"} for \eqn{\beta_x^{(0)}=1}, a predefined parametric
#'   function of age (see details) or \code{NULL} if there is no cohort effect.
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
#' @param  constFunEst function defining the identifiability constraints for the
#'   equivalent mortality rate model. It must be a function of the form
#'   \code{constFunEst <- function(Ax, ax, bx, kt, b0x, gc, wxt, ages)} taking a set
#'   of fitted model parameters and returning a list
#'   \code{list(Ax = Ax, ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
#'   of the model parameters with the identifiability constraints applied. If
#'   omitted no identifiability constraints are applied to the estimaion model.
#'
#' @export
iMoMo  <- function(staticAgeFun = TRUE, periodAgeFun = 'NP',
                    cohortAgeFun = NULL, type = c("indirect", "direct"),
                    constFun = function(ax, bx, kt, b0x, gc, wxt, ages)
                      list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc),
                    constFunEst = function(Ax, ax, bx, kt, b0x, gc, wxt, ages)
                      list(Ax = Ax, ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)) {

  type <- match.arg(type)
  #Create estimation model
  model <- StMoMo(link = "log",
                  staticAgeFun = ifelse(type == "indirect", TRUE, staticAgeFun),
                  periodAgeFun = periodAgeFun, cohortAgeFun = cohortAgeFun)

  #Add constant improvement rates if necessary
  if (type == "indirect" && staticAgeFun){
    model$gnmFormula <- paste(model$gnmFormula, "factor(x):t", sep = " + ")
  }


  #Structure of the model is the same but on improvement rates
  textFormula <- sub("log m\\[x,t\\]", "eta\\[x,t\\]", model$textFormula)
  if (!staticAgeFun)
    textFormula <- sub("a\\[x\\] \\+ ", "", textFormula)


  out <- list(staticAgeFun = staticAgeFun,
              periodAgeFun = model$periodAgeFun,
              cohortAgeFun = model$cohortAgeFun,
              N = model$N,
              textFormula = textFormula,
              type = type,
              constFun = constFun,
              constFunEst = constFunEst,
              model = model)


  #Estimation type
  if (type == "indirect")
    class(out) <-  c("iMoMo", "iMoMoI", "StMoMo")
  else
    class(out) <-  c("iMoMo", "iMoMoD", "StMoMo")
  out

}


#' @export
print.iMoMo <- function(x, ...) {
  cat(paste(x$type, "model with predictor: "))
  cat("")
  cat(x$textFormula)
}


