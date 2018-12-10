#' Transform fitted parameters of direct improvement
#' rated model into the appropriate output format
#' @keywords internal
transParamiMoMoD <- function(out){
  if(!is.null(out$kt)) out$kt <- out$kt * -1
  if(!is.null(out$gc)) out$gc <- out$gc * -1
  if(!is.null(out$ax)) out$ax <- out$ax * -1
  out
}



#' Transform fitted parameters of observed improvement
#' rated model into the appropriate output format
#' @keywords internal
transParamiMoMoO <- function(out){
  if(!is.null(out$kt)) out$kt <- out$kt * -1
  if(!is.null(out$gc)) out$gc <- out$gc * -1
  if(!is.null(out$ax)) out$ax <- out$ax * -1
  out
}



#' Transform fitted parameters of fitted improvement
#' rated model into the appropriate output format
#' @keywords internal
transParamiMoMoF <- function(out){
  out$years <- out$years[2:length(out$years)]
  out$cohorts <- out$cohorts[2:length(out$cohorts)]
  if (out$model$N > 0) {
    for (i in 1:out$model$N) {
      ci <- out$kt[i, 1]
      out$kt[i, ] <- out$kt[i, ] - ci
      out$kt[i, ] <- c(0, diff(out$kt[i, ]))
      out$ax <- out$ax + ci * out$bx[, i]

    }
    out$kt <- matrix(out$kt[,-1], nrow = out$model$N, ncol = length(out$years),
                     dimnames = list(1:out$model$N, out$years))
  }

  if (!is.null(out$gc)) {
    c0 <- out$gc[1]
    out$gc <- out$gc - c0
    out$ax <- out$ax + c0 * out$b0x
    out$gc <- diff(out$gc)
  }
  out$Ax <-  out$ax
  out$ax <- NULL
  out <- c(out, list(ax=NULL))
  if(!is.null(out$kt)) out$kt <- out$kt * -1
  if(!is.null(out$gc)) out$gc <- out$gc * -1
  out
}


#' Transform fitted parameters of fitted improvement
#' rated model into the appropriate output format
#' @keywords internal
transParamiMoMoCIF <- function(out){
  #Get the model parameters
  ages <- out$ages
  nAges <- length(ages)
  kt <- out$kt
  bx <- out$bx
  b0x <- out$b0x
  gc <- out$gc
  #Get static age function and constant improvements
  axTemp <- out$ax[grep(pattern = "^[-]?[[:digit:]]+$",
                        names(out$ax))]
  dxTemp <- out$ax[grep(pattern = ":t",
                        names(out$ax))]
  names(dxTemp) <- sub(pattern = ":t", replacement = "" ,
                       x = names(dxTemp))
  ax <- rep(0,nAges)
  names(ax) <- ages
  ax[names(axTemp)] <- axTemp
  ax[is.na(ax)] <- 0

  dx <- rep(0,nAges)
  names(dx) <- ages
  dx[names(dxTemp)] <- dxTemp
  dx[is.na(dx)] <- 0

  #Apply trasnformations
  t0 <- out$years[1]
  ax <- ax +  t0 * dx

  constPar <- out$model$constFun2(ax = ax, bx = bx, kt = kt,
                               b0x = b0x, gc = gc, wxt = out$wxt,
                               ages = ages)
  constPar <- out$model$constFun3(ax = constPar$ax, dx = dx, bx = constPar$bx,
                                  kt = constPar$kt, b0x = constPar$b0x,
                                  gc = constPar$gc, wxt = out$wxt,
                                  ages = ages)
  out$ax <- constPar$ax
  dx <- constPar$dx
  out$bx <- constPar$bx
  out$kt <- constPar$kt
  out$b0x <- constPar$b0x
  out$gc <- constPar$gc
  if (out$model$N > 0){
    constPar <- constRemoveTrends(ax = constPar$ax, bx = constPar$bx,
                                  kt = constPar$kt, dx = dx)

    out$ax <- constPar$ax
    out$bx <- constPar$bx
    out$kt <- constPar$kt
    dx <- constPar$dx
  }
  #Transform to improvement rate setting
  out <- transParamiMoMoF(out)
  out$ax <- -dx

#   #Apply initial model constraints
#   constPar <- out$model$constFun2(ax = out$ax, bx = out$bx, kt = out$kt,
#                                   b0x = out$b0x, gc = out$gc, wxt = out$wxt[,-1],
#                                   ages = ages)
#   out$ax <- constPar$ax
#   out$bx <- constPar$bx
#   out$kt <- constPar$kt
#   out$b0x <- constPar$b0x
#   out$gc <- constPar$gc


  if(is.null(out$kt)) out <- c(out, list(kt=NULL))
  if(is.null(out$bx)) out <- c(out, list(bx=NULL))
  if(is.null(out$b0x)) out <- c(out, list(b0x=NULL))
  if(is.null(out$gc)) out <- c(out, list(gc=NULL))

  out
}

#' Remove trends from kt
#' @keywords internal
constRemoveTrends <- function(ax, bx, kt, dx){
  n <- dim(kt)[2]
  K <- dim(kt)[1]
  tt <- 1:n
  for (i in 1:K){
    phi <- coef(lm(kt[i, ] ~ tt))
    kt[i, ] <- kt[i, ] - phi[1] - phi[2] * tt
    ax <- ax + phi[1] * bx[, i]
    dx <- dx + phi[2] * bx[, i]
  }
  list(ax = ax, bx = bx, kt = kt, dx = dx)
}

