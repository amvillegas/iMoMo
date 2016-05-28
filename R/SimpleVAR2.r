#******************************************************************************************
# Filename: SimpleVAR2.r
# Programmer: Andres Villegas
#             andresmauriciovillegas@gmail.com
#             Andres.Villegas.1@cass.city.ac.uk
#             Andres.Villegas@hymans.co.uk
# Date created: 13/03/2014
# Modification list:
#            - 30/07/2014: Code updated to allow for linear trends
#
# Description: Class for estimating, projecting and simulating Vector Auto Regressive models
#
#              Delta^d y_t = C+Dt+\sum_{i=1}^p A_i Delta^d y_{t-i} + e_t
#******************************************************************************************

#------------------------------------------------------------------------------------------
#                                 ESTIMATION
#------------------------------------------------------------------------------------------
#' Simple VAR estimation
#'
#' y: a matrix containing the multivariate timeseries
#' p: Integer for the lag order (default is p=1)
#' d: Degree of differencing (default is d=0)
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

      # include.constant<- (type == "const")  | (type == "both")
      # include.drift<- (type == "trend")  | (type == "both")
      # fittingModel <- Arima(dy[,1], order = c(p, 0, 0),
      #                       include.drift = include.drift,
      #                       include.constant = include.constant, method = "CSS")
      # A <- NULL
      # if (p > 0){
      #   for(i in 1:p){
      #     A[[i]] <- array(fittingModel$coef[i], dim = c(1,1))
      #   }
      # }
      # if(type %in% c("const", "both")){
      #   C<-fittingModel$coef["intercept"]        #Intercept
      # }else{
      #   C<-rep(0,K)
      # }
      # if(type %in% c("both", "trend")){
      #   D<-fittingModel$coef["drift"]         #trend
      # }else{
      #   D<-rep(0,K)
      # }
      # Sigma <- array(var(resid(fittingModel)), dim = c(1,1))

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


#------------------------------------------------------------------------------------------
#                                 FORECAST
#------------------------------------------------------------------------------------------
#' Forecast a VAR from class simpleVAR
#'
#'
#'  ARGUMENT
#'
#' object: An object of class "ets", "Arima" or "ar".
#' n.ahead: Number of periods for the simulated series
#' colnames: Name of the columns (time periods)
#' @export
forecast.simpleVAR2 <- function(object, h = 10, colnames=NULL){
  n.ahead <- h
  # Forecast a VAR from class simpleVAR
  #
  #
  #  ARGUMENT
  #
  # object: An object of class "ets", "Arima" or "ar".
  # n.ahead: Number of periods for the simulated series
  # colnames: Name of the columns (time periods)

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

#------------------------------------------------------------------------------------------
#                                 SIMULATE
#------------------------------------------------------------------------------------------

#' Simulate n.path of VAR from class simpleVAR
#'
#'
#'  ARGUMENT
#'
#' object: An object of class "ets", "Arima" or "ar".
#' n.ahead: Number of periods for the simulated series
#' n.path: Number of paths to be simulated
#' colnames: Name of the columns (time periods)
#' @export
simulate.simpleVAR2 <- function(object, nsim = 10, seed = NULL, colnames=NULL){
  n.path <- 1
  n.ahead <- nsim
  # Simulate n.path of VAR from class simpleVAR
  #
  #
  #  ARGUMENT
  #
  # object: An object of class "ets", "Arima" or "ar".
  # n.ahead: Number of periods for the simulated series
  # n.path: Number of paths to be simulated
  # colnames: Name of the columns (time periods)

  A <- object$A
  C <- object$C
  D <- object$D
  p <- object$p
  K <- object$K
  d <- object$d
  n <- object$n
  onePathSim<-function(k){
    #generate innovations
    u <- mvrnorm(n=n.ahead,mu=rep(0,object$K),Sigma=object$Sigma)

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
  #run for all the simulations
  onePathSim()

}
