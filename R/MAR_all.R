#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom stats rnorm
NULL

#' A funtion to generate a multivariate autoregressive process (MAR) in time
#' series
#'
#' The function \code{MAR} is used for generating MAR model(s) for examples
#' of the functions \code{SNSeg_Uni}, \code{SNSeg_Multi}, and \code{SNSeg_HD}.
#'
#' @param n the size (length) of time series to be generated
#' @param reptime the number of time series to be generated
#' @param rho value of autocorrelation
#'
#' @return Returns a matrix of the simulated MAR processes. The number of columns
#' of this matrix is equivalent to the value of input argument \code{reptime}, and
#' the number of rows is the value of input argument \code{n}.
#'
#' @examples
#' MAR(n = 1000, reptime = 2, rho = -0.7)
#'
#' @export MAR
MAR <- function(n, reptime, rho){
  inter <- matrix(0, n+30, reptime)
  epsilon <- matrix(rnorm((n+30)*reptime,0,1),(n+30),reptime)
  for (j in 1:(n+29)){
    inter[j+1,1:reptime] <- epsilon[j+1,1:reptime]+rho*inter[j,1:reptime]
  }
  return(inter[31:(n+30),1:reptime])
}

#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom stats rnorm
NULL

#' A funtion to generate a multivariate autoregressive process (MAR) model in
#' time series for testing change points based on variance and
#' autocovariance
#'
#' The function \code{MAR_Variance} is used for generating MAR model(s) for
#' examples  of the functions \code{SNSeg_Uni}, \code{SNSeg_Multi}, and \code{SNSeg_HD}.
#'
#' @param reptime The number of time series to be generated
#' @param type The type of time series for simulation, which includes V1, V2, V3
#' , A1, A2 and A3. The V-beginnings are for testing the variance, and the
#' A-beginnings are for testing the autocorrelation. The simulated time series
#' come from supplement of Zhao et al. (2022) <doi:10.1111/rssb.12552>.
#' Default \code{type} is \code{V3}.
#'
#' The time length and "true change-points locations" (cps) for each \code{type} are as follows:
#' \code{V1}: cps at 400 and 750 with a time length of 1024.
#' \code{V2}: cps at 125, 532 and 704 with a time length of 1024.
#' \code{V3}: cps at 512 and 768 with a time length of 1024.
#' \code{A1}: cps at 400 and 750 with a time length of 1024.
#' \code{A2}: cps at 50 with a time length of 1024.
#' \code{A3}: cps at 512 and 768 with a time length of 1024.
#'
#' @return Returns a matrix of the simulated MAR processes. The number of columns
#' of this matrix is equivalent to the value of input argument \code{reptime}.
#'
#' @examples
#' MAR_Variance(reptime = 2, type = "V1")
#'
#' @export MAR_Variance
MAR_Variance <- function(reptime, type='V3'){
  burnin <- 50
  if(type=='V1'){
    n <- 1024
    inter <- matrix(0, n+burnin, reptime)
    epsilon <- matrix(rnorm((n+burnin)*reptime,0,1),(n+burnin),reptime)
    cp_sets <- c(0,400,750,1024)
    no_seg <- length(cp_sets)-1
    rho_sets <- list(c(0.5),c(0.5),c(0.5))
    sd_sets <- list(c(1),c(2),c(1))
    for(j in 1:(burnin-1)){
      inter[j+1,1:reptime] <- sd_sets[[1]]*epsilon[j+1,1:reptime]+rho_sets[[1]]*inter[j,1:reptime]
    }
    for(index in 1:no_seg){ # Mean shift
      tau1 <- cp_sets[index]+burnin
      tau2 <- cp_sets[index+1]+burnin-1
      for(j in tau1:tau2){
        inter[j+1,1:reptime] <- sd_sets[[index]]*epsilon[j+1,1:reptime]+rho_sets[[index]]%*%inter[j:(j-length(rho_sets[[index]])+1),1:reptime]
      }
    }
  }
  if(type=='V2'){
    n <- 1024
    inter <- matrix(0, n+burnin, reptime)
    epsilon <- matrix(rnorm((n+burnin)*reptime,0,1),(n+burnin),reptime)
    cp_sets <- c(0,125,532,704,1024)
    no_seg <- length(cp_sets)-1
    rho_sets <- list(c(0.7),c(0.3),c(0.9),c(0.1))
    sd_sets <- list(c(1,0.6),c(1,0.3),c(1),c(1,-0.5))
    for(j in 2:(burnin-1)){
      inter[j+1,1:reptime] <- sd_sets[[1]]%*%epsilon[(j+1):(j-length(sd_sets[[1]])+2),1:reptime]+rho_sets[[1]]*inter[j,1:reptime]
    }
    for(index in 1:no_seg){ # Mean shift
      tau1 <- cp_sets[index]+burnin
      tau2 <- cp_sets[index+1]+burnin-1
      for(j in tau1:tau2){
        inter[j+1,1:reptime] <- sd_sets[[index]]%*%epsilon[(j+1):(j-length(sd_sets[[index]])+2),1:reptime]+rho_sets[[index]]*inter[j,1:reptime]
      }
    }
  }
  if(type=='V3'){
    n <- 1024
    inter <- matrix(0, n+burnin, reptime)
    epsilon <- matrix(rnorm((n+burnin)*reptime,0,1),(n+burnin),reptime)
    cp_sets <- c(0,512,768,1024)
    no_seg <- length(cp_sets)-1
    rho_sets <- list(c(0.9), c(1.69,-0.81), c(1.32,-0.81))
    for(j in 1:(burnin-1)){
      inter[j+1,1:reptime] <- epsilon[j+1,1:reptime]+rho_sets[[1]]*inter[j,1:reptime]
    }
    for(index in 1:no_seg){ # Mean shift
      tau1 <- cp_sets[index]+burnin
      tau2 <- cp_sets[index+1]+burnin-1
      for(j in tau1:tau2){
        inter[j+1,1:reptime] <- epsilon[j+1,1:reptime]+rho_sets[[index]]%*%inter[j:(j-length(rho_sets[[index]])+1),1:reptime]
      }
    }
  }
  if(type=='A1'){
    n <- 1024
    inter <- matrix(0, n+burnin, reptime)
    epsilon <- matrix(rnorm((n+burnin)*reptime,0,1),(n+burnin),reptime)
    cp_sets <- c(0,400,750,1024)
    no_seg <- length(cp_sets)-1
    rho_sets <- list(c(0.5),c(0.9),c(0.3))
    for(j in 1:(burnin-1)){
      inter[j+1,1:reptime] <- epsilon[j+1,1:reptime]+rho_sets[[1]]*inter[j,1:reptime]
    }
    for(index in 1:no_seg){ # Mean shift
      tau1 <- cp_sets[index]+burnin
      tau2 <- cp_sets[index+1]+burnin-1
      for(j in tau1:tau2){
        inter[j+1,1:reptime] <- epsilon[j+1,1:reptime]+rho_sets[[index]]%*%inter[j:(j-length(rho_sets[[index]])+1),1:reptime]
      }
    }
  }
  if(type=='A2'){
    n <- 1024
    inter <- matrix(0, n+burnin, reptime)
    epsilon <- matrix(rnorm((n+burnin)*reptime,0,1),(n+burnin),reptime)
    cp_sets <- c(0,50,1024)
    no_seg <- length(cp_sets)-1
    rho_sets <- list(c(0.75),c(-0.5))
    for(j in 1:(burnin-1)){
      inter[j+1,1:reptime] <- epsilon[j+1,1:reptime]+rho_sets[[1]]*inter[j,1:reptime]
    }
    for(index in 1:no_seg){ # Mean shift
      tau1 <- cp_sets[index]+burnin
      tau2 <- cp_sets[index+1]+burnin-1
      for(j in tau1:tau2){
        inter[j+1,1:reptime] <- epsilon[j+1,1:reptime]+rho_sets[[index]]%*%inter[j:(j-length(rho_sets[[index]])+1),1:reptime]
      }
    }
  }
  if(type=='A3'){
    n <- 1024
    inter <- matrix(0, n+burnin, reptime)
    epsilon <- matrix(rnorm((n+burnin)*reptime,0,1),(n+burnin),reptime)
    cp_sets <- c(0,512,768,1024)
    no_seg <- length(cp_sets)-1
    rho_sets <- list(c(-0.9),c(0.9),c(0))
    sd_sets <- list(c(1,0.7),c(1),c(1,-0.7))
    for(j in 2:(burnin-1)){
      inter[j+1,1:reptime] <- sd_sets[[1]]%*%epsilon[(j+1):(j-length(sd_sets[[1]])+2),1:reptime]+rho_sets[[1]]*inter[j,1:reptime]
    }
    for(index in 1:no_seg){ # Mean shift
      tau1 <- cp_sets[index]+burnin
      tau2 <- cp_sets[index+1]+burnin-1
      for(j in tau1:tau2){
        inter[j+1,1:reptime] <- sd_sets[[index]]%*%epsilon[(j+1):(j-length(sd_sets[[index]])+2),1:reptime]+rho_sets[[index]]*inter[j,1:reptime]
      }
    }
  }
  return(inter[(burnin+1):(n+burnin),1:reptime])
}

#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom mvtnorm rmvnorm
NULL

#' A Funtion to generate a multivariate autoregressive process (MAR) model in
#' time series. It is used for testing change-points based on the change in multivariate
#' means or multivariate covariance for multivariate time series. It also works
#' for the change in correlations between two univariate time series.
#'
#' The function \code{MAR_MTS_Covariance} is used to generate MAR model(s) for
#' examples of the functions \code{SNSeg_Uni}, \code{SNSeg_Multi}, and \code{SNSeg_HD}.
#'
#' @param n the size of time series to be generated.
#' @param reptime the number of time series to be generated.
#' @param rho_sets autocorrelations for each univariate time series.
#' @param cp_sets numeric values of the true change-point locations (0, change-point
#' locations and the end point).
#' @param sigma_cross a list of matrices to generate the multivariate covariance
#' matrices.
#'
#' @returns Returns a list of matrices where each matrix is a MAR process. The
#' number of columns for each sub-matrix is equivalent to the value of input
#' argument \code{reptime}.
#'
#' @examples
#' n <- 1000
#' reptime <- 2
#' sigma_cross <- list(4*matrix(c(1,0.8,0.8,1), nrow=2),
#'                       matrix(c(1,0.2,0.2,1), nrow=2),
#'                       matrix(c(1,0.8,0.8,1), nrow=2))
#' cp_sets <- round(c(0,n/3,2*n/3,n))
#' noCP <- length(cp_sets)-2
#' rho_sets <- rep(0.5, noCP+1)
#' MAR_MTS_Covariance(n, reptime, rho_sets, cp_sets, sigma_cross)
#'
#' @export MAR_MTS_Covariance
MAR_MTS_Covariance <- function(n, reptime, rho_sets, cp_sets, sigma_cross){
  no_ts <- dim(sigma_cross[[1]])[1]
  ts_sim <- list()
  burnin <- 100
  no_seg <- length(cp_sets)-1
  for(rep_index in 1:reptime){
    epsilon <- rmvnorm(burnin, mean=rep(0,no_ts), sigma=sigma_cross[[1]])
    inter <- matrix(0, n+burnin, no_ts)
    for (j in 1:(burnin-1)){
      inter[j+1,] <- epsilon[j+1,]+rho_sets[1]*inter[j,]
    }
    for(index in 1:no_seg){ # Mean shift
      tau1 <- cp_sets[index]+burnin
      tau2 <- cp_sets[index+1]+burnin-1
      epsilon <- rmvnorm(tau2-tau1+1, mean=rep(0,no_ts), sigma=sigma_cross[[index]])
      for(j in tau1:tau2){
        inter[j+1,] <- epsilon[(j+1-tau1),]+rho_sets[index]*inter[j,]
      }
    }
    ts_sim[[rep_index]] <- t(inter[(burnin+1):(n+burnin),])
  }
  return(ts_sim)
}


