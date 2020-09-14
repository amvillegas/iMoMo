#' Create a Lee-Carter model for improvment rates
#'
#' Utility function to initialise a \code{iMoMo} object representing a
#' Lee-Carter model on improvment rates.
#'
#' The created model  has predictor structure
#' \deqn{\log(\mu_{xt}/\mu_{x,t-1}) = \alpha_x + \beta_x\kappa_t.}
#' To ensure identifiability one of the  following constraints is imposed
#' \deqn{\sum_t\kappa_t = 0,\,\kappa_1 = 0,\, \kappa_n = 0}
#' depending on the value of \code{const}, and
#' \deqn{\sum_x\beta_x = X.}
#' where \eqn{X}  is the number of ages in the data
#'
#' @inheritParams iMoMo
#' @param const defines the constraint to impose to the period index of the
#'  model to ensure identifiability. The alternatives are
#'  \code{"sum"}(default),  \code{"last"} and \code{"first"} which apply
#'  constraints \eqn{\sum_t\kappa_t = 0}, \eqn{\kappa_n = 0} and
#'  \eqn{\kappa_1 = 0}, respectively.
#'
#' @return An object of class \code{"StMoMo"}.
#'
#' @seealso \code{\link{iMoMo}}
#'
#' @export
lci <- function(type = c("indirect", "direct"), const = c("sum", "last", "first")) {
  type <- match.arg(type)
  const <- match.arg(const)
  constLC <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    N <- length(ages)
    c1 <- switch(const, sum = mean(kt[1, ], na.rm = TRUE),
                 first = kt[1, 1], last = tail(kt[1, ], 1))
    ax <- ax + c1 * bx[, 1]
    kt[1, ] <- kt[1, ] - c1
    c2 <- sum(bx[, 1], na.rm = TRUE)
    bx[, 1] <- bx[, 1] / c2 * N
    kt[1, ] <- kt[1, ] * c2 / N
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  iMoMo(staticAgeFun = TRUE, periodAgeFun = 'NP', type = type,
         constFun = constLC)
}


#' Create a Cairns-Blake-Dowd model for improvement rates
#'
#' @inheritParams iMoMo
#'
#' @export
cbdi <- function(type = c("indirect", "direct")) {
  type <- match.arg(type)
  f1 <- function(x,ages) x - mean(ages)
  iMoMo(type = type, staticAgeFun = FALSE, periodAgeFun=c("1", f1))
}


#' Create an Age-Period-Cohort model for improvement rates
#'
#' @inheritParams iMoMo
#'
#' @export
apci <- function(type = c("indirect", "direct")) {
  type <- match.arg(type)
  constAPC <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    nYears <- dim(kt)[2]
    x <- ages
    t <- 1:nYears
    c <- (1 - tail(ages, 1)):(nYears - ages[1])
    #\sum g(c)=0  and  \sum cg(c)=0
    phiReg <- lm(gc ~ 1 + c, na.action = na.omit)
    phi <- coef(phiReg)
    gc <- gc - phi[1] - phi[2] * c
    ax <- ax + phi[1] - phi[2] *x
    kt <- kt + phi[2] * t
    #\sum k(t)=0
    c1 <- mean(kt, na.rm = TRUE)
    kt <- kt - c1
    ax <- ax + c1
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  iMoMo(type = type, staticAgeFun = TRUE, periodAgeFun = "1",
         cohortAgeFun = "1", constFun = constAPC)
}
