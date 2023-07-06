#' Critical Values of Self-Normalization (SN) based test statistic for the
#' change in a single parameter (SNCP)
#'
#' A dataset containing the critical value for SN-based change point estimates
#' based on the change in a single parameter.
#'
#' @name critical_values_single
#' @docType data
#' @format A data frame with 6 variables:
#' \describe{
#' \item{\code{epsilon}}{value used to compute grid_size_scale and SN-based
#' test statistic}
#' \item{\code{0.9}}{critical value at confidence level \code{0.9}}
#' \item{\code{0.95}}{critical value at confidence level \code{0.95}}
#' \item{\code{0.99}}{critical value at confidence level \code{0.99}}
#' \item{\code{0.995}}{critical value at confidence level \code{0.995}}
#' \item{\code{0.999}}{critical value at confidence level \code{0.999}}
#' }
#'
"critical_values_single"

#' Critical Values of Self-Normalization (SN) based test statistic for changes
#' in multiple parameters (SNCP)
#'
#' A dataset containing the critical value of SN-based change point estimates
#' based on simultaneous changes in multiple parameters.
#'
#' @name critical_values_multi
#' @docType data
#' @format A data frame with 7 variables:
#' \describe{
#' \item{\code{epsilon}}{value used to compute grid_size_scale and SN-based
#' test statistic}
#' \item{\code{p}}{dimension of the multi-parameters}
#' \item{\code{0.9}}{critical value at confidence level \code{0.9}}
#' \item{\code{0.95}}{critical value at confidence level \code{0.95}}
#' \item{\code{0.99}}{critical value at confidence level \code{0.99}}
#' \item{\code{0.995}}{critical value at confidence level \code{0.995}}
#' \item{\code{0.999}}{critical value at confidence level \code{0.999}}
#' }
#'
"critical_values_multi"

#' Critical Values of Self-Normalization (SN) based test statistic for changes
#' in high-dimensional means (SNHD)
#'
#' A dataset containing the critical value of SN-based change point estimates
#' based on changes in high-dimensional means.
#'
#' @name critical_values_HD
#' @docType data
#' @format A data frame with 6 variables:
#' \describe{
#' \item{\code{epsilon}}{value used to compute grid_size_scale and SN-based
#' test statistic}
#' \item{\code{0.9}}{critical value at confidence level \code{0.9}}
#' \item{\code{0.95}}{critical value at confidence level \code{0.95}}
#' \item{\code{0.99}}{critical value at confidence level \code{0.99}}
#' \item{\code{0.995}}{critical value at confidence level \code{0.995}}
#' \item{\code{0.999}}{critical value at confidence level \code{0.999}}
#' }
#'
"critical_values_HD"
