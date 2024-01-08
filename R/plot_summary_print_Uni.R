#' Plotting the output for univariate or bivariate time series (testing
#' the change in correlation between bivariate time series)
#'
#' Plotting method for S3 objects of class \code{SNSeg_Uni}
#' @method plot SNSeg_Uni
#' @param x a \code{SNSeg_Uni} object
#' @param cpts.col a specification for the color of the vertical lines at
#' the change point estimators, see \link[graphics]{par}
#' @param ... additional graphical arguments, see \link[graphics]{plot}
#' and \link[graphics]{abline}. Users are allowed to enter their own title for
#' the univariate time series plot. The bivariate time series does not contain
#' this option.
#' @details
#' The location of each change point estimator is plotted as a vertical line
#' against the input time series.
#' @examples
#' \donttest{
#' set.seed(7)
#' ts <- MAR_Variance(2, "V1")
#' ts <- ts[,2]
#' # test the change in a single parameter (variance)
#' # grid_size defined
#' result <- SNSeg_Uni(ts, paras_to_test = "variance", confidence = 0.9,
#'                     grid_size_scale = 0.05, grid_size = 67,
#'                     plot_SN = FALSE, est_cp_loc = TRUE)
#' plot(result, cpts.col='red')
#' }
#'
#' @importFrom graphics abline plot par
#' @export
plot.SNSeg_Uni <- function(x, cpts.col='red', ...) {
  ts <- x$ts
  if(min(dim(as.matrix(ts))) == 1){ # univariate time series
    plot(ts, type = 'l', xlab = 'Time', ylab = "Value", ...)
    abline(h=x$est_cp, col=cpts.col)
  }else{ # bivariate time series for change in bivariate correlation
    if(dim(ts)[2] == 2) {ts <- t(ts)}
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(mfrow=c(1,2))
    plot(ts[1,], xlab = "Time", ylab = "Value", type = 'l',
         main="SN Segmentation plot for Time Series 1")
    abline(h=x$est_cp, col=cpts.col)
    plot(ts[2,], xlab = "Time", ylab = "Value", type = 'l',
         main="SN Segmentation plot for Time Series 2")
    abline(h=x$est_cp, col=cpts.col)
  }
}

#' Summary of SN-based change-point estimates for univariate or bivariate time
#' series (testing the change in correlation between bivariate time series)
#'
#' Summary method for objects of class \code{SNSeg_Uni}
#' @method summary SNSeg_Uni
#' @param object a \code{SNSeg_Uni} object
#' @param ... not in use
#' @details Provide information about estimated change-point locations, the
#' parameter tested by SN-based procedures, the confidence level, the \code{grid_size},
#' and the critical value of the SN-based test.
#' @examples
#' \donttest{
#' set.seed(7)
#' ts <- MAR_Variance(2, "V1")
#' ts <- ts[,2]
#' # test the change in a single parameter (variance)
#' # grid_size defined
#' result <- SNSeg_Uni(ts, paras_to_test = "variance", confidence = 0.9,
#'                     grid_size_scale = 0.05, grid_size = 67,
#'                     plot_SN = FALSE, est_cp_loc = TRUE)
#' summary(result)
#' }
#'
#' @export
summary.SNSeg_Uni <- function(object, ...) {
  cpts <- object$est_cp
  para <- object$paras_to_test
  if(length(cpts) > 0){
    cat(paste0('There are ', length(cpts), ' change-point(s) detected at ', 100*object$confidence, 'th confidence level based on the change in '))
    if(length(para) == 1){
      if(inherits(para, 'function')){
        cat(paste0('the user-defined functional.\n'))
      }else if(inherits(para, 'character')){
        cat(paste0('the single ', para, ' parameter.\n'))
      }else{
        cat(paste0('the ', 100*para, 'th quantile.\n'))
      }
    }else{
      cat(paste0('multiple parameters.\n'))
      cat('\n')
      cat(paste0('The parameters being tested are '))
      for(p in para){
        if(is.na(suppressWarnings(as.numeric(p)))){
          cat(paste0(p), ', ')
        }else{
          cat(paste0(100*as.numeric(p), 'th quantile, '))
        }
      }
      cat('\n')
    }
    cat('\n')
    cat(paste('The critical value of SN-based test is', object$critical_value))
    cat('\n')
    cat('\n')
    cpts <- paste(cpts, collapse = ',')
    cat(paste('The detected change-point location(s) are', cpts, 'with a grid_size of', object$grid_size))
  }else{
    cat(paste('No change-point was detected based on the change in '))
    if(length(para) == 1){
      if(inherits(para, 'function')){
        cat(paste0('the user-defined functional.\n'))
      }else if(inherits(para, 'character')){
        cat(paste0('the single ', para, ' parameter.\n'))
      }else{
        cat(paste0('the ', 100*para, 'th quantile.\n'))
      }
    }else{
      cat(paste0('multiple parameters.\n'))
      cat('\n')
      cat(paste0('The parameters being tested are '))
      for(p in para){
        if(inherits(p, 'character')){
          cat(paste0(p), ', ')
        }else{
          cat(paste0(100*p, 'th quantile, '))
        }
      }
      cat('\n')
    }
    cat('\n')
    cat(paste('The critical value of SN-based test is', object$critical_value))
  }
}

#' Print SN-based change-point estimates for univariate or bivariate time
#' series (testing the change in correlation between bivariate time series)
#'
#' Print method for objects of class \code{SNSeg_Uni}
#' @method print SNSeg_Uni
#' @param x a \code{SNSeg_Uni} object
#' @param ... not in use
#' @examples
#' \donttest{
#' set.seed(7)
#' ts <- MAR_Variance(2, "V1")
#' ts <- ts[,2]
#' # test the change in a single parameter (variance)
#' # grid_size defined
#' result <- SNSeg_Uni(ts, paras_to_test = "variance", confidence = 0.9,
#'                     grid_size_scale = 0.05, grid_size = 67,
#'                     plot_SN = FALSE, est_cp_loc = TRUE)
#' print(result)
#' }
#'
#' @export
print.SNSeg_Uni <- function(x, ...) {
  cpts <- x$est_cp
  if(length(cpts)>0){
    cpts <- paste(cpts, collapse = ',')
    cat(paste('The detected change-point location(s) are', cpts))
  }else{
    cat(paste('No change-point was found'))
  }
}
