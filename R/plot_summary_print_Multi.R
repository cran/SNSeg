#' Plotting the output for multivariate time series with dimension no greater
#' than 10
#'
#' Plotting method for S3 objects of class \code{SNSeg_Multi}
#' @method plot SNSeg_Multi
#' @param x a \code{SNSeg_Multi} object
#' @param cpts.col a specification for the color of the vertical lines at
#' the change point estimators, see \link[graphics]{par}
#' @param ... additional graphical arguments, see \link[graphics]{plot}
#' and \link[graphics]{abline}
#' @details
#' The location of each change point estimator is plotted as a vertical line
#' against the input time series.
#' @examples
#' \donttest{
#' # Please run this function before simulation
#' exchange_cor_matrix <- function(d, rho){
#'   tmp <- matrix(rho, d, d)
#'   diag(tmp) <- 1
#'   return(tmp)
#' }
#'
#' # simulation of multivariate time series
#' library(mvtnorm)
#' set.seed(10)
#' d <- 5
#' n <- 600
#' nocp <- 5
#' cp_sets <- round(seq(0, nocp+1 ,1)/(nocp+1)*n)
#' mean_shift <- rep(c(0,2),100)[1:(length(cp_sets)-1)]/sqrt(d)
#' rho_sets <- 0.2
#' sigma_cross <- list(exchange_cor_matrix(d,0))
#' ts <- MAR_MTS_Covariance(n, 2, rho_sets, cp_sets = c(0,n), sigma_cross)
#' ts <- ts[1][[1]]
#'
#' # Test for the change in multivariate means
#' # grid_size defined
#' result <- SNSeg_Multi(ts, paras_to_test = "mean", confidence = 0.99,
#'                       grid_size_scale = 0.05, grid_size = 45)
#' # plot method
#' plot(result)
#'
#' }
#'
#' @importFrom graphics abline plot
#' @export
plot.SNSeg_Multi <- function(x, cpts.col='red', ...) {
  ts <- x$ts
  if(dim(ts)[1]>dim(ts)[2]){
    ts <- t(ts)
  }
  for(i in 1:dim(ts)[1]){
    plot(ts[i,], type='l', xlab = "Time", ylab = "Value",
         main=paste("SN Segmentation plot for Time Series", i))
    abline(h=x$est_cp, col=cpts.col)
  }
}

#' Summary of SN-based change-point estimates for multivariate time series with
#' dimension no greater than 10
#'
#' Summary method for objects of class \code{SNSeg_Multi}
#' @method summary SNSeg_Multi
#' @param object a \code{SNSeg_Multi} object
#' @param ... not in use
#' @details Provide information about estimated change-point locations, the
#' parameter tested by SN-based procedures, the confidence level, the \code{grid_size},
#' and the critical value of the SN-based test.
#' @examples
#' \donttest{
#' # Please run this function before simulation
#' exchange_cor_matrix <- function(d, rho){
#'   tmp <- matrix(rho, d, d)
#'   diag(tmp) <- 1
#'   return(tmp)
#' }
#'
#' # simulation of multivariate time series
#' library(mvtnorm)
#' set.seed(10)
#' d <- 5
#' n <- 600
#' nocp <- 5
#' cp_sets <- round(seq(0, nocp+1 ,1)/(nocp+1)*n)
#' mean_shift <- rep(c(0,2),100)[1:(length(cp_sets)-1)]/sqrt(d)
#' rho_sets <- 0.2
#' sigma_cross <- list(exchange_cor_matrix(d,0))
#' ts <- MAR_MTS_Covariance(n, 2, rho_sets, cp_sets = c(0,n), sigma_cross)
#' ts <- ts[1][[1]]
#'
#' # Test for the change in multivariate means
#' # grid_size defined
#' result <- SNSeg_Multi(ts, paras_to_test = "mean", confidence = 0.99,
#'                       grid_size_scale = 0.05, grid_size = 45)
#' # summary method
#' summary(result)
#' }
#'
#' @export
summary.SNSeg_Multi <- function(object, ...) {
  cpts <- object$est_cp
  para <- object$paras_to_test
  if(length(cpts) > 0){
    cat(paste0('There are ', length(cpts), ' change-points detected at ', 100*object$confidence,
               'th confidence level based on the change in multivariate ', cat(para)))
    cat('\n')
    cat('\n')
    cpts <- paste(cpts, collapse = ',')
    cat(paste('The detected change-point location(s) are', cpts, 'with a grid_size of', object$grid_size))
    cat('\n')
    cat('\n')
    cat(paste('The critical value of SN-based test is', object$critical_value))
  }else{
    cat(paste('No change-point was detected based on the change in multivariate ', cat(para)))
    cat('\n')
    cat('\n')
    cat(paste('The critical value of SN-based test is', object$critical_value))
  }
}

#' Print SN-based change-point estimates for multivariate time series with
#' dimension no greater than 10
#'
#' Print method for objects of class \code{SNSeg_Multi}
#' @method print SNSeg_Multi
#' @param x a \code{SNSeg_Multi} object
#' @param ... not in use
#' @examples
#' \donttest{
#' # Please run this function before simulation
#' exchange_cor_matrix <- function(d, rho){
#'   tmp <- matrix(rho, d, d)
#'   diag(tmp) <- 1
#'   return(tmp)
#' }
#'
#' # simulation of multivariate time series
#' library(mvtnorm)
#' set.seed(10)
#' d <- 5
#' n <- 600
#' nocp <- 5
#' cp_sets <- round(seq(0, nocp+1 ,1)/(nocp+1)*n)
#' mean_shift <- rep(c(0,2),100)[1:(length(cp_sets)-1)]/sqrt(d)
#' rho_sets <- 0.2
#' sigma_cross <- list(exchange_cor_matrix(d,0))
#' ts <- MAR_MTS_Covariance(n, 2, rho_sets, cp_sets = c(0,n), sigma_cross)
#' ts <- ts[1][[1]]
#'
#' # Test for the change in multivariate means
#' # grid_size defined
#' result <- SNSeg_Multi(ts, paras_to_test = "mean", confidence = 0.99,
#'                       grid_size_scale = 0.05, grid_size = 45)
#' # print method
#' print(result)
#' }
#'
#' @export
print.SNSeg_Multi <- function(x, ...) {
  cpts <- x$est_cp
  if(length(cpts)>0){
    cpts <- paste(cpts, collapse = ',')
    cat(paste('The detected change-point location(s) are', cpts))
  }else{
    cat(paste('No change-point was found'))
  }
}
