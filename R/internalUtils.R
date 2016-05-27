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
