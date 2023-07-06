#' SNSeg: An R Package for Time Series Segmentation via Self-Normalization (SN)
#'
#' The SNSeg package provides three functions for multiple change point
#' estimation using SN-based algorithms: \code{SNSeg_Uni}, \code{SNSeg_Multi} and \code{SNSeg_HD}.
#' Three critical value tables (\code{critical_values_single},
#' \code{critical_values_multi} and \code{critical_values_HD}) were attached.
#' Functions \code{MAR}, \code{MAR_Variance} and \code{MAR_MTS_Covariance} can be utilized
#' to generate time series data that are used for the functions \code{SNSeg_Uni}, \code{SNSeg_Multi} and \code{SNSeg_HD}.
#' The function \code{max_SNsweep} enables users to compute the SN test
#' statistic and make the segmentation plot for these statistics.
#'
#' @section SNSeg_Uni:
#' \code{SNSeg_Uni} provides SN-based change point estimates for a univariate
#' time series based on changes in a single parameter or multiple parameters.
#'
#' For the parameters of the SN test, the function
#' \code{SNSeg_Uni} offers mean, variance, acf, bivariate
#' correlation and numeric quantiles as available options. To visualize the
#' estimated change points, users can set "plot_SN = TRUE" and "est_cp_loc = TRUE"
#' to generate the time series segmentation plot. The output comprises of the
#' parameter(s), the window size, and the estimated change point locations.
#'
#' @section SNSeg_Multi:
#' \code{SNSeg_Multi} provides SN-based change point estimates for multivariate
#' time series based on changes in multivariate means or covariance matrix.
#' Different from the function \code{SNSeg_Uni}, \code{SNSeg_Multi} does not
#' contain the option to generate the time series segmentation plot.
#'
#' @section SNSeg_HD:
#' \code{SNSeg_HD} provides SN-based change point estimates for a
#' high-dimensional time series based on changes in high-dimensional means.
#'
#' @section critical values table:
#' The package \code{SNSeg} provides three critical values table.
#'
#' Table \code{critical_values_single} tabulates critical values of SN-based
#' change point estimates based on the change in a single parameter.
#'
#' Table \code{critical_values_multi} tabulates critical values of SN-based
#' change point estimates based on changes in multiple parameters.
#'
#' Table \code{critical_values_HD} tabulates critical values of of SN-based
#' change point estimates based on changes in high-dimensional means.
#'
#' @docType package
#' @name SNSeg
NULL
#> NULL
