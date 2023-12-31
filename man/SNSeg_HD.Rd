% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNSeg_HD.R
\name{SNSeg_HD}
\alias{SNSeg_HD}
\title{Self-normalization (SN) based change points estimation for high dimensional
time series for changes in high-dimensional means (SNHD).}
\usage{
SNSeg_HD(ts, confidence = 0.9, grid_size_scale = 0.05, grid_size = NULL)
}
\arguments{
\item{ts}{A high-dimensional time series represented as a matrix with p
columns, where each column is a univariate time series. The dimension p for
ts should be at least 10.}

\item{confidence}{Confidence level of SN tests as a numeric value. Available
choices of confidence levels contain 0.9, 0.95, 0.99, 0.995 and 0.999. The
default is set to 0.9.}

\item{grid_size_scale}{numeric value of the trimming parameter and only in
use if grid_size = NULL. Users are allowed to choose any grid_size_scale
between 0.05 and 0.5. A warning will be given if it is outside the
range.}

\item{grid_size}{Local window size h to compute the critical value for SN
test. Since grid_size = n*grid_size_scale, where n is the length of time
series, this function will compute the grid_size_scale by diving n from
grid_size when it is not NULL.}
}
\value{
SNSeg_HD returns a list of numeric objects, including the local window
size to cover a change point, the estimated change-point locations, the
confidence level and the critical value of the SN test.
\describe{
\item{\code{grid_size}}{A numeric value of the window size.}
\item{\code{SN_sweep_result}}{A list of n matrices where each matrix
consists of four columns: (1) SN-based test statistic for each change-point
location (2) Change-point location  (3) Lower bound of the window h and
(4) Upper bound of the window h.}
\item{\code{est_cp}}{A vector containing the locations of the estimated
change-points.}
\item{\code{confidence}}{Confidence level of SN test as a numeric value.}
\item{\code{critical_value}}{Critical value of the SN-based test statistic.}
}

For more examples of \code{SNSeg_HD} see the help vignette:
\code{vignette("SNSeg", package = "SNSeg")}
}
\description{
The function \code{SNSeg_HD} is a SNHD change point
estimation procedure.
}
\examples{
\donttest{
n <- 500
p <- 50
nocp <- 5
cp_sets <- round(seq(0,nocp+1,1)/(nocp+1)*n)
num_entry <- 5
kappa <- sqrt(4/5)
mean_shift <- rep(c(0,kappa),100)[1:(length(cp_sets)-1)]
set.seed(1)
ts <- matrix(rnorm(n*p,0,1),n,p)
no_seg <- length(cp_sets)-1
for(index in 1:no_seg){
  tau1 <- cp_sets[index]+1
  tau2 <- cp_sets[index+1]
  ts[tau1:tau2,1:num_entry] <- ts[tau1:tau2,1:num_entry] +
    mean_shift[index]
}

# grid_size defined
result <- SNSeg_HD(ts, confidence = 0.9, grid_size_scale  = 0.05,
                   grid_size = 40)
# Estimated change-point locations
result$est_cp

# For more examples, please run the following command:
# vignette("SNSeg", package = "SNSeg")
}

}
