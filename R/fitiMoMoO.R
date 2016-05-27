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
  out$fittingModel <- out

  out <- transParamiMoMoO(out)

#   if(!is.null(out$kt)) out$kt <- out$kt * -1
#   if(!is.null(out$gc)) out$gc <- out$gc * -1
#   if(!is.null(out$ax)) out$ax <- out$ax * -1
  out$Dxt <- Dxt[which(ages %in% ages.fit),
                 which(years %in% years.fit)]
  out$Ext <- Ext[which(ages %in% ages.fit),
                 which(years %in% years.fit)]
  class(out$model) <- class(object)
  class(out) <- c("fitiMoMoO", "fitiMoMo", class(out))
  out
}
