#' Fit an Improvement Rate Mortality Model of type observed
#'
#' @export
fit.iMoMoO <- function(object, Dxt, Ext, ages = 1:nrow(Dxt),
                       years = 1:ncol(Dxt), ages.fit = ages,
                       years.fit = years, wxt = NULL,
                       start.ax = NULL, start.bx = NULL, start.kt = NULL,
                       start.b0x = NULL, start.gc = NULL, verbose = TRUE,
                       ...) {



  n <- length(years.fit)
  mxt <- Dxt/Ext
  oxt <- log(mxt[which(ages %in% ages.fit),
                 which(years %in% years.fit[1:(n-1)])])
  tempObject <- object
  class(tempObject) <- "StMoMo"
  out <- fit(tempObject, Dxt, Ext, ages = ages,
           years = years, ages.fit = ages.fit,
           years.fit = years.fit[-1], oxt = oxt, wxt = wxt,
           start.ax = start.ax, start.bx = start.bx, start.kt = start.kt,
           start.b0x = start.b0x, start.gc = start.gc, verbose = verbose,
           ...)
  class(out$model) <- class(object)
  class(out) <- c("fitiMoMoO", class(out))
  out
}
