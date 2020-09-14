#' Simple VAR estimation
#'
#'Function to fit VAR mole of the form:
#' \deqn{\Delta^d y_t = C+Dt+\sum_{i=1}^p A_i \Delta^d y_{t-i} + \epsilon_t}
#' where \eqn{C} and  \eqn{D} are \eqn{K}-dimensional vectors for parameters and
#' \eqn{A_1,...,A_p} are \eqn{K\times K} matrices of autoregressive parampeters.
#'
#' @param y a matrix containing the multivariate timeseries
#' @param p integer for the lag order (default is \code{p=1})
#' @param d degree of differencing (default is \code{d=0})
#' @param type string indicating whether a constand a trend are included. If
#' \code{"const"} only \eqn{C} is included. If \code{"trend"} only \eqn{D} is
#' included.  If \code{"both"} both \eqn{C} and \eqn{D} are included.
#'
#' @return A list of class \code{"simpleVAR2"} with components:
#'
#' \item{K}{Dimension of the vector \eqn{y}. \code{K=1} if univariate}
#' \item{A}{A 3-dimensional array with the vector autregressive matrices, \eqn{A_1,...,A_p}.}
#' \item{C}{Vector of constant parameters.}
#' \item{D}{Vector of trend parametesr.}
#' \item{p}{Integer indicating the lag order.}
#' \item{d}{Degree of differencing.}
#' \item{n}{number of time periods used in the fitting.}
#' \item{Sigma}{Variance-covariance matrix.}
#' \item{resid}{Residuals of the model.}
#' \item{fittingModel}{Output of the model call underlying the fitting.}
#' \item{y}{Data used in the fitting}
#' \item{dy}{Difference data used in the fitting}
#' \item{datamat}{Array including all hte lags used in fitting the model.}
#'
#' @export
simpleVAR2 <- function(y,p=1,d=0, type = c("const", "trend", "both", "none")){

  type <- match.arg(type)
  if(is.array(y)){
    n<-dim(y)[1]
    K<-dim(y)[2]
  }else{
    n <- length(y)
    K<-1
    y <- array(y,c(n,K),dimnames=list(names(y)))
  }

  # if(K<2){
  #   stop("The matrix 'y' should contain at least two variables. For univariate analysis consider ar() and arima() in package stats.\n")
  # }

  if(p!=0){
    #difference the data
    dy <- y
    if(d>0) dy <- diff(y,differences=d)
    datamat <- NULL
    dy_i <- y
    for (i in 1:(d+1)){
      datamat[[i]]<- dy_i
      dy_i <- diff(dy_i)
    }

    #fit the time series model
    if (K < 2){
      xconst <- rep(1, n-d)
      xtrend <- (d+1):n
      yreg <- dy[(p+1):(n-d), 1]
      xreg <- NULL
      for (i in 1:p) {
        xreg <- cbind(xreg, dy[(p+1-i):(n - i - d), 1])
      }
      if ((type == "const")  | (type == "both")) {
        xreg <- cbind(xreg, xconst[(p+1):(n-d)])
      }
      if ((type == "trend")  | (type == "both")) {
        xreg <- cbind(xreg, xtrend[(p+1):(n-d)])
      }
      fittingModel <- lm(yreg ~ -1 + xreg)
      A <- NULL
      if (p > 0){
        for(i in 1:p){
          A[[i]] <- array(fittingModel$coef[i], dim = c(1,1))
        }
      }
      if(type %in% c("const", "both")){
        C<-fittingModel$coef[p + 1]        #Intercept
      }else{
        C<-rep(0,K)
      }
      if(type %in% c("both", "trend")){
        D<-fittingModel$coef[p + 1 + (type == "both")]         #trend
      }else{
        D<-rep(0,K)
      }
      Sigma <- array(var(resid(fittingModel)), dim = c(1,1))

    } else {
      fittingModel <- VAR(dy,p=p,type=type)
      #prepare the output
      A <- vars::Acoef(fittingModel) #Autoregressive matrices
      B <- vars::Bcoef(fittingModel)
      if(type %in% c("const", "both")){
        C<-B[,"const"]         #Intercept
      }else{
        C<-rep(0,K)
      }
      if(type %in% c("both", "trend")){
        D<-B[,"trend"]         #Intercept
      }else{
        D<-rep(0,K)
      }
      Sigma <- cov(resid(fittingModel))
    }

    out<-list(K=K,A=A,C=C,D=D,p=p,d=d,n=n,Sigma=Sigma,
              resid=resid(fittingModel),fittingModel=fittingModel,
              y=y,dy=dy,datamat=datamat)
  }else{
    if(type %in% c("const")){ #Random walk with drift
      out<-simpleVAR(y=y,p=p,d=d)
      out$D <- rep(0,K)
      out$n <- n
      if(K > 1){
        out$Sigma <- cov(out$resid)
      } else {
        out$Sigma <- var(out$resid)
      }
    }
    else{
      stop("Ramdom walks and random walks with trends not handled at the moment\n")
    }
  }
  class(out)<-"simpleVAR2"
  out
}


#' Forecast a VAR from class simpleVAR2
#'
#' Function for obtaininng forcast of a VAR model
#'
#' @param object An object of class \code{"simpleVAR2."}.
#' @param h Number of periods for the forecasted series
#' @param colnames Name of the columns (time periods)
#'
#' @return A matrix with the forecast series
#'
#' @export
forecast.simpleVAR2 <- function(object, h = 10, colnames=NULL){
  n.ahead <- h

  A <- object$A
  C <- object$C
  D <- object$D
  p <- object$p
  K <- object$K
  d <- object$d
  n <- object$n

  #forecast the differenced series
  dyt_for <- array(0,dim=c(n.ahead,K))
  dy_t <-  tail(object$dy,object$p)
  for (t in 1:n.ahead){
    dyt_for[t,] <- C + D*(n+t)
    i <- 1
    while(i<=p){
      dyt_for[t,] <- dyt_for[t,] + t(A[[i]]%*%array(dy_t[p+1-i,],dim=c(K,1)))
      i<-i+1
    }
    i <- 1
    while(i<p){
      dy_t[i,] <- dy_t[i+1,]
      i<-i+1
    }
    dy_t[p,]<-dyt_for[t,]

  }
  #undifference the data
  j <- d
  yt_for<-dyt_for
  while(j>0){
    yt_t_1 <- tail(object$datamat[[j]],1)
    for(t in 1:n.ahead){
      yt_for[t,] <- yt_for[t,] + yt_t_1
      yt_t_1 <- yt_for[t,]
    }
    j<-j-1
  }
  dimnames(yt_for) <- list(colnames)
  yt_for
}


#' Simulate a VAR model
#'
#' Returns one simulated path of the VAR model in \code{object}.
#'
#' @param object: An object of class \code{"simpleVAR2"}.
#' @param nsim number of periods for the simulated series.
#' @param seed either \code{NULL} or an integer that will be used in a
#' call to \code{\link{set.seed}} before simulating the time series.
#' The default, \code{NULL} will not change the random generator state.
#' @param colnames Name of the columns (time periods)
#' @export
simulate.simpleVAR2 <- function(object, nsim = 10, seed = NULL, colnames=NULL){

  if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)
  if (is.null(seed))
    RNGstate <- .Random.seed
  else {
    R.seed <- .Random.seed
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  n.path <- 1
  n.ahead <- nsim
  A <- object$A
  C <- object$C
  D <- object$D
  p <- object$p
  K <- object$K
  d <- object$d
  n <- object$n

  #generate innovations
  u <- mvrnorm(n = n.ahead, mu = rep(0, object$K), Sigma=object$Sigma)

  #generate the differenced series
  dyt_sim <- array(0,dim=dim(u))
  dy_t <-  tail(object$dy,object$p)
  for (t in 1:n.ahead){
    dyt_sim[t,] <- C + D*(n+t) + u[t,]
    i <- 1
    while(i<=p){
      dyt_sim[t,] <- dyt_sim[t,] + t(A[[i]]%*%array(dy_t[p+1-i,],dim=c(K,1)))
      i<-i+1
    }
    i <- 1
    while(i<p){
      dy_t[i,] <- dy_t[i+1,]
      i<-i+1
    }
    dy_t[p,]<-dyt_sim[t,]

  }

  #undifference the data
  j <- d
  yt_sim<-dyt_sim
  while(j>0){
      yt_t_1 <- tail(object$datamat[[j]],1)
      for(t in 1:n.ahead){
        yt_sim[t,] <- yt_sim[t,] + yt_t_1
        yt_t_1 <- yt_sim[t,]
      }
      j<-j-1
  }
  dimnames(yt_sim) <- list(colnames)
  yt_sim

}
