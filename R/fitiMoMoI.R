#' Fit an Improvement Rate Model indirectly
#'
#' Fit an Improvement Rate Model indirectly via fiting the equivalent mortality
#' rate model and then trasforming the output to improvement rates
#'
#' @param object an object of class \code{"iMoMoI"} defining the improvement rate
#' model.
#' @param data an optional object of type StMoMoData containing information on
#' deaths and exposures to be used for fitting the model. This is typically created
#' with  function \code{\link{StMoMoData}}. If this is not provided then the fitting
#' data is taken from arguments, \code{Dxt}, \code{Ext}, \code{ages}, \code{years}.
#' @param Dxt optional matrix of deaths data.
#' @param Ext optional matrix of observed exposures of the same dimension of
#' \code{Dxt}.
#' @param ages optional vector of ages corresponding to rows of \code{Dxt} and
#' \code{Ext}.
#' @param years optional vector of years corresponding to rows of \code{Dxt} and
#' \code{Ext}.
#' @param ages.fit optional vector of ages to include in the fit. Must be a
#' subset of \code{ages}.
#' @param years.fit optional vector of years to include in the fit. Must be a
#' subset of \code{years}.
#' @param wxt optional matrix of 0-1 weights to be used in the fitting process.
#' This can be used, for instance, to zero weight some cohorts in the data.
#' See \code{\link{genWeightMat}} which is a helper function for defining
#' weighting matrices.
#' @param start.ax optional vector with starting values for \eqn{\alpha_x}.
#' @param start.bx optional matrix with starting values for \eqn{\beta_x^{(i)}}.
#' @param start.kt optional matrix with starting values for \eqn{\kappa_t^{(i)}}.
#' @param start.b0x optional vector with starting values for \eqn{\beta_x^{(0)}}.
#' @param start.gc optional vector with starting values for \eqn{\gamma_c}.
#' @param verbose a logical value. If \code{TRUE} progress indicators are
#' printed as the model is fitted. Set \code{verbose = FALSE} to silent the
#' fitting and avoid progress messages.
#' @param ... arguments to be passed to or from other methods. This can be
#' used to control the fitting parameters of \code{gnm}. See
#' \code{\link[gnm]{gnm}}.
#'
#' @export
fit.iMoMoI <- function(object, data = NULL, Dxt = NULL, Ext = NULL,
                       ages = NULL, years = NULL, ages.fit = NULL,
                       years.fit = NULL, wxt = NULL,
                       start.ax = NULL, start.bx = NULL, start.kt = NULL,
                       start.b0x = NULL, start.gc = NULL, verbose = TRUE,
                       ...) {


  # Select data from data or from Dxt, Ext, ages, years
  if(!is.null(data)) {
    if (class(data) != "StMoMoData")
      stop("Argument data needs to be of class StMoMoData.")
    Dxt <- data$Dxt
    Ext <- data$Ext
    ages <- data$ages
    years <- data$years
  } else {
    if (is.null(Dxt) || is.null(Ext))
      stop("Either argument data or arguments Dxt and Ext need to be provided.")
    if (is.null(ages)) ages <- 1:nrow(Dxt)
    if (is.null(years)) years <- 1:ncol(Dxt)
    data <- structure(list(Dxt = Dxt, Ext = Ext, ages = ages, years = years,
                           type = ifelse(object$link == "log", "central", "initial"),
                           series = "unknown", label = "unknown"), class = "StMoMoData")
  }
  if (is.null(ages.fit)) ages.fit <- ages
  if (is.null(years.fit)) years.fit <- years



  out <- fit(object$model, Dxt = Dxt, Ext = Ext, ages = ages,
             years = years, ages.fit = ages.fit,
             years.fit = years.fit, wxt = wxt,
             start.ax = start.ax, start.bx = start.bx, start.kt = start.kt,
             start.b0x = start.b0x, start.gc = start.gc, verbose = verbose,
             ...)
  out$fittingModel <- out

  if(object$staticAgeFun == FALSE){
    out <- transParamiMoMoI(out)
  } else {

  }

  #Apply identfiability constraints
  constPar <- object$constFun(out$ax, out$bx, out$kt,
                              out$b0x, out$gc, wxt, ages.fit)
  out$ax <- constPar$ax
  out$bx <- constPar$bx
  out$kt <- constPar$kt
  out$b0x <- constPar$b0x
  out$gc <- constPar$gc

  #Prepare output
  out$model <- object
  class(out) <- c("fitiMoMoI", "fitiMoMo", "fitStMoMo")
  out
}
