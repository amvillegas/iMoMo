#' Fit an Improvement Rate Mortality Model of type fitted
#'
#' @export
fit.iMoMoF <- function(object, Dxt, Ext, ages = 1:nrow(Dxt),
                       years = 1:ncol(Dxt), ages.fit = ages,
                       years.fit = years, wxt = NULL,
                       start.ax = NULL, start.bx = NULL, start.kt = NULL,
                       start.b0x = NULL, start.gc = NULL, verbose = TRUE,
                       ...) {

  #Fit the equivalent mortality rate model
  tempObject <- object
  class(tempObject) <- "StMoMo"
  out <- fit(tempObject, Dxt, Ext, ages = ages,
             years = years, ages.fit = ages.fit,
             years.fit = years.fit, wxt = wxt,
             start.ax = start.ax, start.bx = start.bx, start.kt = start.kt,
             start.b0x = start.b0x, start.gc = start.gc, verbose = verbose,
             ...)
  out$fittingModel <- out
  #Transform the parameters
#   out$years <- out$years[2:length(out$years)]
#   out$cohorts <- out$cohorts[2:length(out$cohorts)]
#   if (object$N > 0) {
#     for (i in 1:object$N) {
#       ci <- out$kt[i, 1]
#       out$kt[i, ] <- out$kt[i, ] - ci
#       out$kt[i, ] <- c(0, diff(out$kt[i, ]))
#       out$ax <- out$ax + ci * out$bx[, i]
#
#     }
#     out$kt <- matrix(out$kt[,-1], nrow = object$N, ncol = length(out$years),
#                      dimnames = list(1:object$N, out$years))
#   }
#
#   if (!is.null(out$gc)) {
#       c0 <- out$gc[1]
#       out$gc <- out$gc - c0
#       out$ax <- out$ax + c0 * out$b0x
#       out$gc <- diff(out$gc)
#   }
#   out$Ax <-  out$ax
#   out$ax <- NULL
#   out <- c(out, list(ax=NULL))
#   if(!is.null(out$kt)) out$kt <- out$kt * -1
#   if(!is.null(out$gc)) out$gc <- out$gc * -1
  out <- transParamiMoMoF(out)
  class(out$model) <- class(object)
  class(out) <- c("fitiMoMoF", "fitiMoMo", "fitStMoMo")
  out
}
