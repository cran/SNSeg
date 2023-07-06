#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @importFrom utils data
#' @importFrom stats approx
NULL

#' SN-based test statistic segmentation plot for univariate, mulitivariate
#' and high-dimensional time series
#'
#' The function \code{max_SNsweep} allows users to compute and plot the SN-based
#' test statistics along with the identified change-points from functions
#' SNSeg_Uni, SNSeg_Multi, or SNSeg_HD.
#'
#' @param SN_result The output of functions SNSeg_Uni, SNSeg_Multi or SNSeg_HD.
#' @param plot_SN A boolean value to return an SN-based segmentation plot if
#' plot_SN = TRUE.
#' @param est_cp_loc A boolean value to plot a red solid vertical line for
#' estimated change-point locations if est_cp_loc = TRUE.
#' @param critical_loc A boolean value to plot a blue dashed horizontal line for
#' the critical value if critical_loc = TRUE
#'
#' @return Returns a vector of numeric values of calculated SN-based statistics
#' for each time point. It also generates a SN-based test statistics segmentation
#' plot with the estimated change-points.
#'
#' For more examples of \code{max_SNsweep} please see the SNSeg vignette:
#' \code{vignette("SNSeg", package = "SNSeg")}
#'
#' @examples
#' \donttest{
#' set.seed(7)
#' n <- 2000
#' reptime <- 2
#' cp_sets <- round(n*c(0,cumsum(c(0.5,0.25)),1))
#' mean_shift <- c(0.4,0,0.4)
#' rho <- -0.7
#' ts <- MAR(n, reptime, rho)
#' no_seg <- length(cp_sets)-1
#' for(index in 1:no_seg){
#'   tau1 <- cp_sets[index]+1
#'   tau2 <- cp_sets[index+1]
#'   ts[tau1:tau2,] <- ts[tau1:tau2,] + mean_shift[index]
#' }
#' ts <- ts[,2]
#' result <- SNSeg_Uni(ts, paras_to_test = "mean", confidence = 0.9,
#'                     grid_size_scale = 0.05, grid_size = 116,
#'                     plot_SN = FALSE, est_cp_loc = FALSE)
#'
#' # Generate SN-based test statistic segmentation plot
#' # To get the computed SN-based statistics, please run the command "test_stat"
#' test_stat <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
#'                          critical_loc = TRUE)
#'
#' # For more examples of \code{max_SNsweep} see the help vignette:
#' # \code{vignette("SNSeg", package = "SNSeg")}
#' }
#'
#' @export max_SNsweep

max_SNsweep <- function(SN_result, plot_SN = TRUE,
                        est_cp_loc = TRUE, critical_loc = TRUE){
  # SN-based test statistic for each time point within the window set
  max_matrix <- function(matrix_data){
    if(is.null(dim(matrix_data))){
      return(0)
    }else{
      return(max(matrix_data[,1]))
    }
  }

  SN_test_stat <- sapply(SN_result$SN_sweep_result, max_matrix)

  if(plot_SN){
    plot(SN_test_stat, xlab = "Time", ylab = "Value",
          main="SN Test Statistic Segmentation Plot")
    if(est_cp_loc){abline(v = SN_result$est_cp, col = 'red')}
    if(critical_loc){abline(h = SN_result$critical_value, col = 'blue', lty = 2)}
  }

  return("SN_test_statistic" = SN_test_stat)
}
