#' Plotting the output for high-dimensional time series with dimension greater
#' than 10
#'
#' Plotting method for S3 objects of class \code{SNSeg_HD}
#' @method plot SNSeg_HD
#' @param x a \code{SNSeg_HD} object
#' @param cpts.col a specification for the color of the vertical lines at
#' the change point estimators, see \link[graphics]{par}
#' @param ts_index The index number(s) of the univariate time series to be plotted.
#' Users should enter a positive integer or a vector of positive integers that are
#' no greater than the dimension of the input time series. The default is the
#' first 5 time series, i.e., \code{ts_index = c(1:5)}.
#' @param ... additional graphical arguments, see \link[graphics]{plot}
#' and \link[graphics]{abline}
#' @details
#' The location of each change point estimator is plotted as a vertical line
#' against the input time series.
#' @examples
#' \donttest{
#' n <- 500
#' p <- 50
#' nocp <- 5
#' cp_sets <- round(seq(0,nocp+1,1)/(nocp+1)*n)
#' num_entry <- 5
#' kappa <- sqrt(4/5)
#' mean_shift <- rep(c(0,kappa),100)[1:(length(cp_sets)-1)]
#' set.seed(1)
#' ts <- matrix(rnorm(n*p,0,1),n,p)
#' no_seg <- length(cp_sets)-1
#' for(index in 1:no_seg){
#'   tau1 <- cp_sets[index]+1
#'   tau2 <- cp_sets[index+1]
#'   ts[tau1:tau2,1:num_entry] <- ts[tau1:tau2,1:num_entry] +
#'     mean_shift[index]
#' }
#'
#' # grid_size defined
#' result <- SNSeg_HD(ts, confidence = 0.9, grid_size_scale  = 0.05,
#'                    grid_size = 40)
#' # plot the 1st, 3rd and 5th time series
#' plot(result, cpts.col = 'red', ts_index = c(1,3,5))
#' }
#'
#' @importFrom graphics abline plot
#' @export
plot.SNSeg_HD <- function(x, cpts.col='red', ts_index = c(1:5), ...) {
  ts <- x$ts
  if(dim(ts)[1]>dim(ts)[2]){
    ts <- t(ts)
  }
  for(i in 1:ts_index){
    plot(ts[i,], type='l', xlab = "Time", ylab = "Value",
         main=paste("SN Segmentation plot for Time Series", i))
    abline(h=x$est_cp, col=cpts.col)
  }
}

#' Summary of SN-based change-point estimates for high-dimensional time series
#' with dimension greater than 10
#'
#' Summary method for objects of class \code{SNSeg_HD}
#' @method summary SNSeg_HD
#' @param object a \code{SNSeg_HD} object
#' @param ... not in use
#' @details Provide information about estimated change-point locations, the
#' parameter tested by SN-based procedures, the confidence level, the \code{grid_size},
#' and the critical value of the SN-based test.
#' @examples
#' \donttest{
#' n <- 500
#' p <- 50
#' nocp <- 5
#' cp_sets <- round(seq(0,nocp+1,1)/(nocp+1)*n)
#' num_entry <- 5
#' kappa <- sqrt(4/5)
#' mean_shift <- rep(c(0,kappa),100)[1:(length(cp_sets)-1)]
#' set.seed(1)
#' ts <- matrix(rnorm(n*p,0,1),n,p)
#' no_seg <- length(cp_sets)-1
#' for(index in 1:no_seg){
#'   tau1 <- cp_sets[index]+1
#'   tau2 <- cp_sets[index+1]
#'   ts[tau1:tau2,1:num_entry] <- ts[tau1:tau2,1:num_entry] +
#'     mean_shift[index]
#' }
#'
#' # grid_size defined
#' result <- SNSeg_HD(ts, confidence = 0.9, grid_size_scale  = 0.05,
#'                    grid_size = 40)
#' # summary method
#' summary(result)
#' }
#'
#' @export
summary.SNSeg_HD <- function(object, ...) {
  cpts <- object$est_cp
  if(length(cpts) > 0){
    cat(paste0('There are ', length(cpts), ' change-points detected at ', 100*object$confidence,
               'th confidence level based on the change in high-dimensional means'))
    cat('\n')
    cat('\n')
    cpts <- paste(cpts, collapse = ',')
    cat(paste('The detected change-point location(s) are', cpts, 'with a grid_size of', object$grid_size))
    cat('\n')
    cat('\n')
    cat(paste('The critical value of SN-based test is', object$critical_value))
  }else{
    cat(paste('No change-point was detected based on the change in high dimensional means'))
    cat('\n')
    cat('\n')
    cat(paste('The critical value of SN-based test is', object$critical_value))
  }
}

#' Print SN-based change-point estimates for high-dimensional time series with
#' dimension greater than 10
#'
#' Print method for objects of class \code{SNSeg_HD}
#' @method print SNSeg_HD
#' @param x a \code{SNSeg_HD} object
#' @param ... not in use
#' @examples
#' \donttest{
#' n <- 500
#' p <- 50
#' nocp <- 5
#' cp_sets <- round(seq(0,nocp+1,1)/(nocp+1)*n)
#' num_entry <- 5
#' kappa <- sqrt(4/5)
#' mean_shift <- rep(c(0,kappa),100)[1:(length(cp_sets)-1)]
#' set.seed(1)
#' ts <- matrix(rnorm(n*p,0,1),n,p)
#' no_seg <- length(cp_sets)-1
#' for(index in 1:no_seg){
#'   tau1 <- cp_sets[index]+1
#'   tau2 <- cp_sets[index+1]
#'   ts[tau1:tau2,1:num_entry] <- ts[tau1:tau2,1:num_entry] +
#'     mean_shift[index]
#' }
#'
#' # grid_size defined
#' result <- SNSeg_HD(ts, confidence = 0.9, grid_size_scale  = 0.05,
#'                    grid_size = 40)
#' # print method
#' print(result)
#' }
#'
#' @export
print.SNSeg_HD <- function(x, ...) {
  cpts <- x$est_cp
  if(length(cpts)>0){
    cpts <- paste(cpts, collapse = ',')
    cat(paste('The detected change-point location(s) are', cpts))
  }else{
    cat(paste('No change-point was found'))
  }
}
