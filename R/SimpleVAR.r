#' Simple VAR estimation
#'
#' Function to fit VAR mole of the form:
#' \deqn{\Delta^d y_t = C + \sum_{i=1}^p A_i \Delta^d y_{t-i} + \epsilon_t}
#' where \eqn{C} is a \eqn{K}-dimensional vector of parameters and
#' \eqn{A_1,...,A_p} are \eqn{K\times K} matrices of autoregressive parampeters.
#'
#' @param y a matrix containing the multivariate timeseries
#' @param p Integer for the lag order (default is p=1)
#' @param d Degree of differencing (default is d=0)
#' @export
simpleVAR <- function(y,p=1,d=0){
  if(is.array(y)){
    n<-dim(y)[1]
    K<-dim(y)[2]
  }else{
    n <- length(y)
    K<-1
    y <- array(y,c(n,K),dimnames=list(names(y)))
  }

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
  fittingModel <- ar.ols(dy,aic=F,order.max=p,demean=FALSE,intercept=TRUE)

  #prepare the output
  A <- NULL
  i<-1
  while(i<=p){
    A[[i]]<-fittingModel$ar[i,,]
    i<-i+1
  }
  out<-list(K=K,A=A,C=fittingModel$x.intercept,p=p,d=d,Sigma=fittingModel$var.pred,
            resid=fittingModel$resid,fittingModel=fittingModel,y=y,dy=dy,datamat=datamat)
  class(out)<-"simpleVAR"
  out

}

#' Forecast a VAR from class simpleVAR
#'
#' Forecast a VAR from class simpleVAR
#'
#' @param object An object of class \code{"simpleVAR"}.
#' @param h Number of periods for the simulated series
#' @param colnames Name of the columns (time periods)
#'
#' @return A matrix with the forecast series
#'
#' @export
forecast.simpleVAR <- function(object, h = 10, colnames=NULL){
  n.ahead <- h
  A <- object$A
  C <- object$C
  p <- object$p
  K <- object$K
  d <- object$d

  #forecast the differenced series
  dyt_for <- array(0,dim=c(n.ahead,K))
  dy_t <-  tail(object$dy,object$p)
  for (t in 1:n.ahead){
    dyt_for[t,] <- C
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
#' @param object: An object of class \code{"simpleVAR1"}.
#' @param nsim number of periods for the simulated series.
#' @param seed either \code{NULL} or an integer that will be used in a
#' call to \code{\link{set.seed}} before simulating the time series.
#' The default, \code{NULL} will not change the random generator state.
#' @param colnames Name of the columns (time periods)
#' @export
simulate.simpleVAR <- function(object, nsim = 10, seed = NULL, colnames=NULL){
  n.path <- 1
  n.ahead <- nsim
  A <- object$A
  C <- object$C
  p <- object$p
  K <- object$K
  d <- object$d

  #generate innovations
  u <- mvrnorm(n=n.ahead,mu=rep(0,object$K),Sigma=object$Sigma)

  #generate the differenced series
  dyt_sim <- array(0,dim=dim(u))
  dy_t <-  tail(object$dy,object$p)
  for (t in 1:n.ahead){
    dyt_sim[t,] <- C + u[t,]
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
