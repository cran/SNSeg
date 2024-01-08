#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @importFrom utils data
#' @importFrom stats approx
#' @importFrom stats cov
#' @importFrom mvtnorm rmvnorm
NULL

#' Self-normalization (SN) based change points estimation for multivariate time
#' series
#'
#' The function \code{SNSeg_Multi} is a SN-based change-points estimation
#' procedure for a multivariate time series based on changes in the multivariate
#' means or covariance matrix.
#'
#' @param ts A multivariate time series represented as a matrix with p columns,
#' where each column is a univariate time series. The dimension p for ts should
#' be at least 2.
#' @param paras_to_test Type of the parameter as a string for which SN
#' algorithms test. Available choices include \code{mean} and \code{covariance}.
#' @param confidence Confidence level of SN tests as a numeric value. Available
#' choices of confidence levels contain 0.9, 0.95, 0.99, 0.995 and 0.999. The
#' default is set to 0.9.
#' @param grid_size_scale  numeric value of the trimming parameter and only in
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
#' estimated change-point locations if est_cp_loc = TRUE
#' @return \code{SNSeg_Multi} returns an S3 object of class "SNSeg_Multi" including
#' the time series, the type of parameter to be tested, the local window size to
#' cover a change point, the estimated change-point locations, the confidence level
#' and the critical value of the SN test. It also generates time series segmentation
#' plot when \code{plot_SN = TRUE}.
#' \describe{
#'   \item{\code{ts}}{A numeric matrix of the input time series.}
#'   \item{\code{paras_to_test}}{the parameter used for the SN test as character.}
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
#' This is a limitation of the function \code{SNSeg_Multi}.
#'
#' For more examples of \code{SNSeg_Multi} see the help vignette:
#' \code{vignette("SNSeg", package = "SNSeg")}
#'
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
#' # Estimated change-point locations
#' result$est_cp
#'
#' # For more examples, please run the following command:
#' # vignette("SNSeg", package = "SNSeg")
#' }
#'
#' @export
SNSeg_Multi <- function(ts, paras_to_test = "mean", confidence = 0.9,
                        grid_size_scale = 0.05, grid_size = NULL,
                        plot_SN = FALSE, est_cp_loc = TRUE){

  if(is.null(ts))
  {stop("Input of ts is missing!")}
  if(!(paras_to_test %in% c("mean","covariance"))){
    stop("paras_to_test must be either mean or covariance!")}
  if(!(confidence %in% c(0.9,0.95,0.99,0.995,0.999))){
    stop("Confidence must be one of 0.9, 0.95, 0.99, 0.995 and 0.999!")
  }

  if(!inherits(ts, 'matrix')){
    stop("ts must be a matrix with dimension at least 2!")
  }
  if(dim(ts)[1] > dim(ts)[2]){ts <- t(ts)}

  n <- dim(ts)[2]
  d <- dim(ts)[1]
  if(paras_to_test=="covariance"){
    d <- d*(d+1)/2
  }
  if(d > 10){
    stop("The dimension of parameters must be smaller than or equal to 10!")
  }

  critical_values_multi <- SNSeg::critical_values_multi
  # interpolation for grid_size_scale
  if((is.null(grid_size))&(is.null(grid_size_scale))){
    grid_size_scale <- 0.05
    SN_critical <- critical_values_multi[critical_values_multi[,2]==d & critical_values_multi[,1]==grid_size_scale,]
    posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
    critical_value <- SN_critical[1,posi_confidence]
    grid_size <- round(grid_size_scale*n)
    warning("Undefined value detected for both grid_size and grid_size_scale! The system would use 0.05 as the default for grid_size_scale.")
  } else if((is.null(grid_size))&(!is.null(grid_size_scale))){
    if(grid_size_scale<0.05){
      grid_size_scale <- 0.05
      SN_critical <- critical_values_multi[critical_values_multi[,2]==d & critical_values_multi[,1]==grid_size_scale,]
      posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
      critical_value <- SN_critical[1,posi_confidence]
      grid_size <- round(grid_size_scale*n)
      warning("Detected the grid_size_scale is less than 0.05. The system would use 0.05 for grid_size_scale.")
    } else if(grid_size_scale>0.5){
      grid_size_scale <- 0.5
      SN_critical <- critical_values_multi[critical_values_multi[,2]==d & critical_values_multi[,1]==grid_size_scale,]
      posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
      critical_value <- SN_critical[1,posi_confidence]
      grid_size <- round(grid_size_scale*n)
      warning("Detected the grid_size_scale is greater than 0.5. The system would use 0.5 for grid_size_scale.")
    } else if(grid_size_scale>=0.05 & grid_size_scale<=0.5){
      if((grid_size_scale %in% critical_values_multi[,1])){
        SN_critical <- critical_values_multi[critical_values_multi[,2]==d & critical_values_multi[,1]==grid_size_scale,]
        posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
        critical_value <- SN_critical[1,posi_confidence]
        grid_size <- round(grid_size_scale*n)
      } else if(!(grid_size_scale %in% critical_values_multi[,1])){
        grid_size <- round(grid_size_scale*n)
        SN_critical <- critical_values_multi[critical_values_multi[,2]==d,]
        posi_epsilon <- match(grid_size_scale,sort(c(SN_critical[,1],grid_size_scale)))
        posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
        critical_vector <- c(SN_critical[(posi_epsilon-1):posi_epsilon,posi_confidence])
        epsilon_vector <- c(SN_critical[(posi_epsilon-1):posi_epsilon,1])
        critical_value <- approx(epsilon_vector,critical_vector,xout = grid_size_scale)$y
      }
    }
  } else if(!is.null(grid_size)){
    grid_size_scale <- grid_size/n
    if(grid_size_scale<0.05){
      grid_size_scale <- 0.05
      grid_size <- round(grid_size_scale*n)
      SN_critical <- critical_values_multi[critical_values_multi[,2]==d & critical_values_multi[,1]==grid_size_scale,]
      posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
      critical_value <- SN_critical[1,posi_confidence]
      warning("Detected the grid_size_scale is less than 0.05 from the current grid_size. The system would use 0.05 for grid_size_scale.")
    }
    if(grid_size_scale>0.5){
      grid_size_scale <- 0.5
      grid_size <- round(grid_size_scale*n)
      SN_critical <- critical_values_multi[critical_values_multi[,2]==d & critical_values_multi[,1]==grid_size_scale,]
      posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
      critical_value <- SN_critical[1,posi_confidence]
      warning("Detected the grid_size_scale is greater than 0.5 from the current grid_size. The system would use 0.5 for grid_size_scale.")
    }
    if(grid_size_scale>=0.05 & grid_size_scale<=0.5){
      if((grid_size_scale %in% critical_values_multi[,1])){
        SN_critical <- critical_values_multi[critical_values_multi[,2]==d & critical_values_multi[,1]==grid_size_scale,]
        posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
        critical_value <- SN_critical[1,posi_confidence]
      } else if(!(grid_size_scale %in% critical_values_multi[,1])){
        SN_critical <- critical_values_multi[critical_values_multi[,2]==d,]
        posi_epsilon <- match(grid_size_scale,sort(c(SN_critical[,1],grid_size_scale)))
        posi_confidence <- match(as.character(confidence),colnames(critical_values_multi))
        critical_vector <- c(SN_critical[(posi_epsilon-1):posi_epsilon,posi_confidence])
        epsilon_vector <- c(SN_critical[(posi_epsilon-1):posi_epsilon,1])
        critical_value <- approx(epsilon_vector,critical_vector,xout = grid_size_scale)$y
      }
    }
  }

  if(paras_to_test == "mean"){
    SN_sweep_result <- SN_sweep_mean_MTS(ts, grid_size)
    SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)
  }

  if(paras_to_test == "covariance"){
    SN_sweep_result <- SN_sweep_covmatrix_MTS(ts, grid_size)
    SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)
  }

  if(plot_SN == TRUE){
    for(i in 1:dim(ts)[1]){
      plot(ts[i,], xlab = "Time", ylab = "Value",
           main=paste0("SN Segmentation plot for Time Series ",i,""))
      if(est_cp_loc){
        abline(v = SN_result, col = 'red')
      }
    }
  }

  final_result <- structure(
    list(
      ts = ts, paras_to_test = paras_to_test, grid_size = grid_size,
      SN_sweep_result = SN_sweep_result, est_cp = SN_result,
      confidence = confidence, critical_value = critical_value
    ), class = 'SNSeg_Multi'
  )
  final_result
}
