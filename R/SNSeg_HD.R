#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Self-normalization (SN) based change points estimation for high dimensional
#' time series for changes in high-dimensional means (SNHD).
#'
#' The function \code{SNSeg_HD} is a SNHD change point
#' estimation procedure.
#'
#' @param ts A high-dimensional time series represented as a matrix with p
#' columns, where each column is a univariate time series. The dimension p for
#' ts should be at least 10.
#' @param confidence Confidence level of SN tests as a numeric value. Available
#' choices of confidence levels contain 0.9, 0.95, 0.99, 0.995 and 0.999. The
#' default is set to 0.9.
#' @param grid_size_scale numeric value of the trimming parameter and only in
#' use if grid_size = NULL. Users are allowed to choose any grid_size_scale
#' between 0.05 and 0.5. A warning will be given if it is outside the
#' range.
#' @param grid_size Local window size h to compute the critical value for SN
#' test. Since grid_size = n*grid_size_scale, where n is the length of time
#' series, this function will compute the grid_size_scale by diving n from
#' grid_size when it is not NULL.
#' @param plot_SN  Boolean value to plot the time series or not.
#' The default setting is FALSE.
#' @param est_cp_loc Boolean value to plot a red solid vertical line for
#' estimated change-point locations if est_cp_loc = TRUE.
#' @param ts_index The index number(s) of the univariate time series to be plotted.
#' Users should enter a positive integer or a vector of positive integers that are
#' no greater than the dimension of the input time series. The default is the
#' first 5 time series, i.e., \code{ts_index = c(1:5)}.
#' @return SNSeg_HD returns an S3 object of class "SNSeg_HD" including the time
#' series, the local window size to cover a change point, the estimated change-point
#' locations, the confidence level and the critical value of the SN test. It also
#' generates time series segmentation plot when \code{plot_SN = TRUE}.
#' \describe{
#'   \item{\code{ts}}{A numeric matrix of the input time series.}
#'   \item{\code{grid_size}}{A numeric value of the window size.}
#'   \item{\code{SN_sweep_result}}{A list of n matrices where each matrix
#'   consists of four columns: (1) SN-based test statistic for each change-point
#'   location (2) Change-point location  (3) Lower bound of the window h and
#'   (4) Upper bound of the window h.}
#'   \item{\code{est_cp}}{A vector containing the locations of the estimated
#'   change-points.}
#'   \item{\code{confidence}}{Confidence level of SN test as a numeric value.}
#'   \item{\code{critical_value}}{Critical value of the SN-based test statistic.}
#' }
#'
#' Users can apply the functions \code{summary.SN} to compute the parameter estimate
#' of each segment separated by the detected change-points. An additional function
#' \code{plot.SN} can be used to plot the time series with estimated change-points.
#' Users can set the option \code{plot_SN = TRUE} or use the function \code{plot.SN}
#' to plot the time series.
#'
#' It deserves to note that some change-points could be missing due to the constraint
#' on \code{grid_size_scale} or related \code{grid_size} that \code{grid_size_scale}
#' has a minimum value of 0.05. Therefore, SNCP claims no change-points within the
#' first n*\code{grid_size_scale} or the last n*\code{grid_size_scale} time points.
#' This is a limitation of the function \code{SNSeg_HD}.
#'
#' For more examples of \code{SNSeg_HD} see the help vignette:
#' \code{vignette("SNSeg", package = "SNSeg")}
#'
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
#' # Estimated change-point locations
#' result$est_cp
#'
#' # For more examples, please run the following command:
#' # vignette("SNSeg", package = "SNSeg")
#' }
#'
#' @export
SNSeg_HD <- function(ts, confidence = 0.9, grid_size_scale = 0.05,
                     grid_size = NULL, plot_SN = FALSE, est_cp_loc = TRUE,
                     ts_index = c(1:5)){
  if(is.null(ts))
  {stop("Input of ts is missing!")}
  if(dim(ts)[1]<dim(ts)[2]){
    ts <- t(ts)
  }

  if(!(confidence %in% c(0.9,0.95,0.99,0.995,0.999))){
    stop("Confidence must be one of 0.9, 0.95, 0.99, 0.995 and 0.999!")
  }
  n <- dim(ts)[1]
  critical_values_HD <- SNSeg::critical_values_HD
  # interpolation for grid_size_scale
  if((is.null(grid_size))&(is.null(grid_size_scale))){
    grid_size_scale <- 0.05
    posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
    critical_value <- critical_values_HD[1,posi_confidence]
    grid_size <- floor(grid_size_scale*n)
    warning("Undefined value detected for both grid_size and grid_size_scale! The system would use 0.05 as the default for grid_size_scale.")
  } else if((is.null(grid_size))&(!is.null(grid_size_scale))){
    if(grid_size_scale<0.05){
      grid_size_scale <- 0.05
      posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
      critical_value <- critical_values_HD[1,posi_confidence]
      grid_size <- floor(grid_size_scale*n)
      warning("Detected the grid_size_scale is less than 0.05. The system would use 0.05 for grid_size_scale.")
    }
    if(grid_size_scale>0.5){
      grid_size_scale <- 0.5
      posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
      critical_value <- critical_values_HD[10,posi_confidence]
      grid_size <- floor(grid_size_scale*n)
      warning("Detected the grid_size_scale is greater than 0.5. The system would use 0.5 for grid_size_scale.")
    }
    if(grid_size_scale>=0.05 & grid_size_scale<=0.5){
      if((grid_size_scale %in% critical_values_HD[,1])){
        posi_cri <- match(as.character(grid_size_scale),critical_values_HD[,1])
        posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
        critical_value <- critical_values_HD[posi_cri,posi_confidence]
        grid_size <- floor(grid_size_scale*n)
      } else if(!(grid_size_scale %in% critical_values_HD[,1])){
        grid_size <- floor(grid_size_scale*n)
        posi_epsilon <- match(grid_size_scale,sort(c(critical_values_HD[,1],grid_size_scale)))
        posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
        critical_vector <- c(critical_values_HD[(posi_epsilon-1):posi_epsilon,posi_confidence])
        epsilon_vector <- c(critical_values_HD[(posi_epsilon-1):posi_epsilon,1])
        critical_value <- approx(epsilon_vector,critical_vector,xout = grid_size_scale)$y
      }
    }
  } else if(!is.null(grid_size)){
    grid_size_scale <- grid_size/n
    if(grid_size_scale<0.05){
      grid_size_scale <- 0.05
      grid_size <- round(grid_size_scale*n)
      posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
      critical_value <- critical_values_HD[1,posi_confidence]
      warning("Detected the grid_size_scale is less than 0.05 from the current grid_size. The system would use 0.05 for grid_size_scale.")
    }
    if(grid_size_scale>0.5){
      grid_size_scale <- 0.5
      grid_size <- round(grid_size_scale*n)
      posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
      critical_value <- critical_values_HD[10,posi_confidence]
      warning("Detected the grid_size_scale is greater than 0.5 from the current grid_size. The system would use 0.5 for grid_size_scale.")
    }
    if(grid_size_scale>=0.05 & grid_size_scale<=0.5){
      if((grid_size_scale %in% critical_values_HD[,1])){
        posi_cri <- match(as.character(grid_size_scale),critical_values_HD[,1])
        posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
        critical_value <- critical_values_HD[posi_cri,posi_confidence]
      } else if(!(grid_size_scale %in% critical_values_HD[,1])){
        posi_epsilon <- match(grid_size_scale,sort(c(critical_values_HD[,1],grid_size_scale)))
        posi_confidence <- match(as.character(confidence),colnames(critical_values_HD))
        critical_vector <- c(critical_values_HD[(posi_epsilon-1):posi_epsilon,posi_confidence])
        epsilon_vector <- c(critical_values_HD[(posi_epsilon-1):posi_epsilon,1])
        critical_value <- approx(epsilon_vector,critical_vector,xout = grid_size_scale)$y
      }
    }
  }

  SN_sweep_result <- SN_sweep_mean_HD(ts, grid_size)
  SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)

  if(plot_SN == TRUE){
    for(i in ts_index){
      plot(ts[,i], xlab = "Time", ylab = "Value",
           main=paste0("SN Segmentation plot for Time Series ",i,""))
      if(est_cp_loc){
        abline(v = SN_result, col = 'red')
      }
    }
  }

  final_result <- structure(
    list(
      ts = ts, grid_size = grid_size,
      SN_sweep_result = SN_sweep_result, est_cp = SN_result,
      confidence = confidence, critical_value = critical_value
    ), class = 'SNSeg_HD'
  )
  final_result

}
