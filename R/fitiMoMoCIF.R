#' Fit an Improvement Rate Mortality Model of type fitted
#'
#' @export
fit.iMoMoCIF <- function(object, Dxt, Ext, ages = 1:nrow(Dxt),
                       years = 1:ncol(Dxt), ages.fit = ages,
                       years.fit = years, wxt = NULL,
                       start.ax = NULL, start.bx = NULL, start.kt = NULL,
                       start.b0x = NULL, start.gc = NULL, verbose = TRUE,
                       ...) {

  #Fit the equivalent mortality rate model
  tempObject <- object
  class(tempObject) <- "StMoMo"
  out <- fit(tempObject, Dxt = Dxt, Ex = Ext, ages = ages,
             years = years, ages.fit = ages.fit,
             years.fit = years.fit, wxt = wxt,
             start.ax = start.ax, start.bx = start.bx, start.kt = start.kt,
             start.b0x = start.b0x, start.gc = start.gc, verbose = verbose,
             ...)
  out$fittingModel <- out
  out$fittingModel$ax <- out$fittingModel$ax[1:length(ages.fit)]
  #Get the model parameters
  # ages <- ages.fit
  # nAges <- length(ages)
  # kt <- out$kt
  # bx <- out$bx
  # b0x <- out$b0x
  # gc <- out$gc
  # #Get static age function and constant improvements
  # axTemp <- out$ax[grep(pattern = "^[-]?[[:digit:]]+$",
  #                                          names(out$ax))]
  # dxTemp <- out$ax[grep(pattern = ":t",
  #                                          names(out$ax))]
  # names(dxTemp) <- sub(pattern = ":t", replacement = "" ,
  #                      x = names(dxTemp))
  # ax <- rep(0,nAges)
  # names(ax) <- ages
  # ax[names(axTemp)] <- axTemp
  #
  # dx <- rep(0,nAges)
  # names(dx) <- ages
  # dx[names(dxTemp)] <- dxTemp
  #
  # #Apply trasnformations
  # t0 <- years.fit[1]
  # ax <- ax +  t0 * dx
  #
  # constPar <- object$constFun2(ax = ax, bx = bx, kt = kt,
  #                           b0x = b0x, gc = gc, wxt = wxt,
  #                           ages = ages.fit)
  #
  # out$ax <- constPar$ax
  # out$bx <- constPar$bx
  # out$kt <- constPar$kt
  # out$b0x <- constPar$b0x
  # out$gc <- constPar$gc
  # if (object$N > 0){
  #   constPar <- constRemoveTrends(ax = constPar$ax, bx = constPar$bx,
  #                                 kt = constPar$kt, dx = dx)
  #
  #   out$ax <- constPar$ax
  #   out$bx <- constPar$bx
  #   out$kt <- constPar$kt
  #   dx <- constPar$dx
  # }
  # #Transform to improvement rate setting
  # out <- transParamiMoMoF(out)
  # out$ax <- -dx

  out <- transParamiMoMoCIF(out)
  class(out$model) <- class(object)
  class(out) <- c("fitiMoMoCIF", "fitiMoMo", "fitStMoMo")
  out
}
