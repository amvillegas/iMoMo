#' Plot bootstrapped parameters of a Stochastic Mortality Model
#'
#' @export
plotTwoBootstraps <- function(x, y, nCol = 2, parametricbx = TRUE,
                            colour = rgb(1, 0, 0), colour.y = rgb(0, 0, 0),
                            probs = c(2.5, 10, 25, 50, 75, 90, 97.5), ...) {
  years <- x$model$years
  ages <- x$model$ages
  cohorts <- x$model$cohorts
  ax <- x$model$ax
  bx <- x$model$bx
  kt <- x$model$kt
  b0x <- x$model$b0x
  gc <- x$model$gc
  N <- x$model$model$N
  fan.col <- colorRampPalette(c(colour, rgb(1, 1, 1)))
  fan.col.y <- colorRampPalette(c(colour.y, rgb(1, 1, 1)))
  n.fan <- length(probs) / 2
  type <- "l"

  #Calculate number of plots and rows
  nPlots <- 2 * N + (!is.null(ax)) + 2 * (!is.null(gc))
  is.nonparametric <- function(bx) {
    is.character(bx) && bx == "NP"
  }
  if (parametricbx == FALSE) {  #substract the parametric plots
    nParametricbx <- ifelse(is.null(x$model$model$periodAgeFun), 0,
                            sum(sapply(x$model$model$periodAgeFun,
                                       FUN = function(x) !is.nonparametric(x)))) +
      (!is.null(x$model$model$cohortAgeFun) && !is.nonparametric(x$model$model$cohortAgeFun))

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
    axBoot <- sapply(X = x$bootParameters, FUN = function(x) x$ax)
    plotParameterFan(x = ages, y = ax, yBoot = axBoot,
                     main = expression(paste(alpha[x], " vs. x", "")),
                     xlab = "age", ylab = "", probs = probs,
                     fan.col = fan.col, n.fan = n.fan, ...)
    axBoot <- sapply(X = y$bootParameters, FUN = function(x) x$ax)
    plotParameterFan(x = ages, y = ax, yBoot = axBoot,
                     main = expression(paste(alpha[x], " vs. x", "")),
                     xlab = "age", ylab = "", probs = probs,
                     fan.col = fan.col.y, n.fan = n.fan, new = FALSE, ...)
  }

  if (!is.null(ax) && nCol == 2 && parametricbx == TRUE) {
    frame()  # add and empty plot to the rigth of ax
  }

  # bx, kt
  if (N > 0) {
    for (i in 1:N) {
      #bx
      if (parametricbx == TRUE || is.nonparametric(x$model$model$periodAgeFun[[i]])) {
        bxBoot <- sapply(X = x$bootParameters, FUN = function(x) x$bx[, i])
        plotParameterFan(x = ages, y = bx[, i], yBoot = bxBoot,
                         main = substitute(paste(beta[x]^{(i)}, " vs. x", ""),
                                           list(i = i)),
                         xlab = "age", ylab = "", probs = probs,
                         fan.col = fan.col, n.fan = n.fan, ...)
        bxBoot <- sapply(X = y$bootParameters, FUN = function(x) x$bx[, i])
        plotParameterFan(x = ages, y = bx[, i], yBoot = bxBoot,
                         main = substitute(paste(beta[x]^{(i)}, " vs. x", ""),
                                           list(i = i)),
                         xlab = "age", ylab = "", probs = probs,
                         fan.col = fan.col.y, n.fan = n.fan, new = FALSE, ...)

      }
      #kt
      ktBoot <- sapply(X = x$bootParameters, FUN = function(x) x$kt[i, ])
      plotParameterFan(x = years, y = kt[i, ], yBoot = ktBoot,
                       main = substitute(paste(kappa[t]^{(i)}, " vs. t", ""),
                                         list(i = i)),
                       xlab = "year", ylab = "", probs = probs,
                       fan.col = fan.col, n.fan = n.fan, ...)
      ktBoot <- sapply(X = y$bootParameters, FUN = function(x) x$kt[i, ])
      plotParameterFan(x = years, y = kt[i, ], yBoot = ktBoot,
                       main = substitute(paste(kappa[t]^{(i)}, " vs. t", ""),
                                         list(i = i)),
                       xlab = "year", ylab = "", probs = probs,
                       fan.col = fan.col.y, n.fan = n.fan, new = FALSE, ...)
    }
  }


  if (!is.null(gc) == TRUE) {
    #bx0
    if (parametricbx == TRUE || is.nonparametric(x$model$model$cohortAgeFun)) {
      b0xBoot <- sapply(X = x$bootParameters, FUN = function(x) x$b0x)
      plotParameterFan(x = ages, y = b0x, yBoot = b0xBoot,
                       main = substitute(paste(beta[x]^{(i)}, " vs. x", ""), list(i = 0)),
                       xlab = "age", ylab = "", probs = probs,
                       fan.col = fan.col, n.fan = n.fan, ...)
      b0xBoot <- sapply(X = y$bootParameters, FUN = function(x) x$b0x)
      plotParameterFan(x = ages, y = b0x, yBoot = b0xBoot,
                       main = substitute(paste(beta[x]^{(i)}, " vs. x", ""), list(i = 0)),
                       xlab = "age", ylab = "", probs = probs,
                       fan.col = fan.col.y, n.fan = n.fan, new = FALSE, ...)
    }
    #gc
    gcBoot <- sapply(X = x$bootParameters, FUN = function(x) x$gc)
    plotParameterFan(x = cohorts, y = gc, yBoot = gcBoot,
                     main = expression(paste(gamma[t-x], " vs. t-x", "")),
                     xlab = "cohort", ylab = "", probs = probs,
                     fan.col = fan.col, n.fan = n.fan, ...)
    gcBoot <- sapply(X = y$bootParameters, FUN = function(x) x$gc)
    plotParameterFan(x = cohorts, y = gc, yBoot = gcBoot,
                     main = expression(paste(gamma[t-x], " vs. t-x", "")),
                     xlab = "cohort", ylab = "", probs = probs,
                     fan.col = fan.col.y, n.fan = n.fan, new = FALSE, ...)

  }
  par(oldpar)
}

#' Plot fanchart of the parameters
#'
#' @keywords internal
plotParameterFan <- function(x, y, yBoot, main, xlab, ylab, probs, fan.col,
                             n.fan, new = TRUE, ...) {
  if(new) plot(NULL, ylab = ylab, xlab = xlab, main = main, xlim = range(x),
       ylim = range(yBoot, finite = TRUE))

  #Compute quantiles
  pp <- apply(t(yBoot), 2, quantile, probs = probs / 100, na.rm = TRUE)

  #find intervals of consecutive not missing values
  n <- length(y)
  i_end <- 0
  while (i_end < n) {
    #start of the interval
    i_start <- i_end + 1
    while (i_start < n && is.na(y[i_start])) {
      i_start <- i_start + 1
    }
    #end of the interval
    i_end <- min(i_start + 1, n)
    while (i_end < n && !is.na(y[i_end + 1])) {
      i_end <- i_end + 1
    }
    if (i_start != i_end) {
      fanplot::fan(pp[, i_start:i_end], data.type = "values",
                   start = x[i_start], probs = probs, fan.col = fan.col,
                   n.fan = n.fan, ln = NULL)

    }
  }
}
