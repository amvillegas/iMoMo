
#' Plot fitted parameters from two stochastic mortality model
#'
#' @param x an object of class \code{"fitStMoMo"} with the fitted
#' parameters of a stochastic mortality model.
#' @param y an object of class \code{"fitStMoMo"} with the fitted
#' parameters of a stochastic mortality model.
#' @param nCol number of columns to use in the plot.
#' @param parametricbx if \code{FALSE} parametric age-modulating terms,
#' which don't need to be estimated, are not plotted.
#' @param ... additional arguments to control graphical appearance.
#' See \code{\link[graphics]{plot}}.
#' @param type what type of plot should be drawn. See
#' \code{\link[graphics]{plot}}.
#'
#' @export
plotTwoModels <- function(x, y, nCol = 2, parametricbx = TRUE, type = "l",
                          legend = NULL,  ...) {

  years <- x$years
  ages <- x$ages
  cohorts <- x$cohorts
  ax <- x$ax
  bx <- x$bx
  kt <- x$kt
  b0x <- x$b0x
  gc <- x$gc
  N <- x$model$N

  #Calculate number of plots and rows
  nPlots <- 2 * N + (!is.null(ax)) + 2 * (!is.null(gc))
  is.nonparametric <- function(bx) {
    is.character(bx) && bx == "NP"
  }
  if (parametricbx == FALSE) {  #substract the parametric plots
    nParametricbx <- ifelse( is.null(x$model$periodAgeFun), 0,
                             sum(sapply(x$model$periodAgeFun,
                                        FUN = function(x) !is.nonparametric(x)))) +
      (!is.null(x$model$cohortAgeFun) && !is.nonparametric(x$model$cohortAgeFun))

    nPlots <- nPlots - nParametricbx
  }
  if (!is.null(ax) && nCol == 2 && parametricbx == TRUE) {
    nPlots <- nPlots + 1  # add and empty plot to the rigth of ax
  }

  nRow <- ceiling(nPlots / nCol)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(nRow, nCol))

  #ax
  if (!is.null(ax)) {
    plot(x = ages, y = ax, ylab = "", xlab = "age",
         main = expression(paste(alpha[x], " vs. x", "")), type = type, ...)
  }

  if (!is.null(ax) && nCol == 2 && parametricbx == TRUE) {
    frame()  # add and empty plot to the rigth of ax
  }

  # bx, kt
  if (N > 0) {
    for (i in 1:N) {
      #bx
      if (parametricbx == TRUE || is.nonparametric(x$model$periodAgeFun[[i]])) {
        matplot(x = ages, y = cbind(bx[, i], y$bx[, i]), ylab = "", xlab = "age",
             main = substitute(paste(beta[x]^{(i)}, " vs. x", ""),
                               list(i = i)), type = type, ...)
      }
      #kt
      matplot(x = years,y = t(rbind(kt[i, ], y$kt[i, ])), ylab = "", xlab = "year",
           main = substitute(paste(kappa[t]^{(i)}, " vs. t", ""),
                             list(i = i)), type = type, ...)
    }
  }


  if (!is.null(gc) == TRUE) {
    #bx0
    if (parametricbx == TRUE || is.nonparametric(x$model$cohortAgeFun)) {
      matplot(x = ages, y = t(rbind(b0x, y$b0x)), ylab = "", xlab = "age",
           main = substitute(paste(beta[x]^{(i)}, " vs. x", ""),
                             list(i = 0)), type = type, ...)
    }
    #gc
    matplot(x = cohorts, y = t(rbind(gc, y$gc)), ylab = "", xlab = "cohort",
         main = expression(paste(gamma[t-x], " vs. t-x","")),
         type = type, ...)
  }
  par(oldpar)
}
