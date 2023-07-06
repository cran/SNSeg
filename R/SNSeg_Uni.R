#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @importFrom utils data
#' @importFrom stats approx
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom stats pnorm
#' @importFrom mvtnorm rmvnorm
NULL

#' Self-normalization (SN) based change point estimates for univariate time
#' series
#'
#' The function \code{SNSeg_Uni} is a SN change point estimation procedure for a
#' univariate time series based on the change in a single or multiple parameters
#' . It also detect changes in correlation between two univariate time series.
#'
#' @param ts A univariate time series expressed as a numeric vector. when the
#' argument paras_to_test is specified as "bivcor", the correlation between
#' bivariate time series, the input ts must be an n by 2 matrix
#' @param paras_to_test The parameters that SN algorithm aim to examine, which
#' are presented as a string, a number, or a combination of both. Available
#' choices of paras_to_test include "mean", "variance", "acf", "bivcor"
#' and a numeric value of quantile between 0 and 1. In the scenario where the
#' input ts is a univariate time series, users are allowed to enter a
#' combination of parameters for paras_to_test except "bivcor".
#' @param confidence Confidence level of SN tests as a numeric value. Available
#' choices of confidence levels contain 0.9, 0.95, 0.99, 0.995 and 0.999. The
#' default is set to 0.9.
#' @param grid_size_scale A numeric value of the trimming parameter and only in
#' use if grid_size = NULL. Users are allowed to choose any grid_size_scale
#' between 0.05 and 0.5. A warning will be given if it is outside the
#' range.
#' @param grid_size Local window size h to compute the critical value for SN
#' test. Since grid_size = n*grid_size_scale, where n is the length of time
#' series, this function will compute the grid_size_scale by dividing n from
#' grid_size when it is not NULL.
#' @param plot_SN  Boolean value to plot the time series or not.
#' The default setting is FALSE.
#' @param est_cp_loc Boolean value to plot a red solid vertical line for
#' estimated change-point locations if est_cp_loc = TRUE
#'
#' @return \code{SNSeg_Uni} returns a list of objects, including the type of
#' parameter to be tested, the local window size to cover a change point, the
#' estimated change-point locations, the confidence level and the critical value
#' of the SN test. It also generates a time series segmentation plot when
#' \code{plot_SN = TRUE}.
#' \describe{
#'   \item{\code{paras_to_test}}{a character, numeric value or vector of the
#'   parameter(s) used for the SN test.}
#'   \item{\code{grid_size}}{A numeric value of the window size.}
#'   \item{\code{SN_sweep_result}}{A list of matrices where each matrix
#'   consists of four columns: (1) SN-based test statistic for each change-point
#'   location (2) Change-point location  (3) Lower bound of the local window and
#'   (4) Upper bound of the local window.}
#'   \item{\code{est_cp}}{A vector containing the locations of the estimated
#'   change-points.}
#'   \item{\code{confidence}}{Confidence level of SN test as a numeric value.}
#'   \item{\code{critical_value}}{Critical value of the SN-based test statistic.}
#' }
#'
#' For more examples of \code{SNSeg_Uni} see the help vignette:
#' \code{vignette("SNSeg", package = "SNSeg")}
#'
#' @examples
#' \donttest{
#' # code to simulate a univariate time series
#' set.seed(7)
#' ts <- MAR_Variance(2, "V1")
#' ts <- ts[,2]
#' # test the change in a single parameter (variance)
#' # grid_size defined
#' result <- SNSeg_Uni(ts, paras_to_test = "variance", confidence = 0.9,
#'                     grid_size_scale = 0.05, grid_size = 67,
#'                     plot_SN = TRUE, est_cp_loc = TRUE)
#' # estimated change-point locations
#' result$est_cp
#' # For more examples of change in a single or multiple parameters, please run
#' # the command: vignette("SNSeg", package = "SNSeg")
#' }
#'
#' @export SNSeg_Uni

SNSeg_Uni <- function(ts, paras_to_test, confidence = 0.9,
                      grid_size_scale = 0.05, grid_size = NULL,
                      plot_SN = TRUE, est_cp_loc = TRUE){
  # combine two SNSeg_Uni functions
  # check if number of parameters > 1
  if(length(paras_to_test) == 1){ # SNSeg_Uni_single_para
    SNSeg_Uni_single_para(ts, type = paras_to_test, confidence = confidence,
                          grid_size_scale = grid_size_scale,
                          grid_size = grid_size, plot_SN = plot_SN,
                          est_cp_loc = est_cp_loc)
  } else{ # SNSeg_Uni_multi_para
    SNSeg_Uni_multi_para(ts, paras_to_test = paras_to_test,
                         confidence = confidence,
                         grid_size_scale = grid_size_scale,
                         grid_size = grid_size, plot_SN = plot_SN,
                         est_cp_loc = est_cp_loc)
  }
}
