% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAR_all.R
\name{MAR_Variance}
\alias{MAR_Variance}
\title{A funtion to generate a multivariate autoregressive process (MAR) model in
time series for testing change points based on variance and
autocovariance}
\usage{
MAR_Variance(reptime, type = "V3")
}
\arguments{
\item{reptime}{The number of time series to be generated}

\item{type}{The type of time series for simulation, which includes V1, V2, V3
, A1, A2 and A3. The V-beginnings are for testing the variance, and the
A-beginnings are for testing the autocorrelation. The simulated time series
come from supplement of Zhao et al. (2022) \url{doi:10.1111/rssb.12552}.
Default \code{type} is \code{V3}.

The time length and "true change-points locations" (cps) for each \code{type} are as follows:
\code{V1}: cps at 400 and 750 with a time length of 1024.
\code{V2}: cps at 125, 532 and 704 with a time length of 1024.
\code{V3}: cps at 512 and 768 with a time length of 1024.
\code{A1}: cps at 400 and 750 with a time length of 1024.
\code{A2}: cps at 50 with a time length of 1024.
\code{A3}: cps at 512 and 768 with a time length of 1024.}
}
\value{
Returns a matrix of the simulated MAR processes. The number of columns
of this matrix is equivalent to the value of input argument \code{reptime}.
}
\description{
The function \code{MAR_Variance} is used for generating MAR model(s) for
examples  of the functions \code{SNSeg_Uni}, \code{SNSeg_Multi}, and \code{SNSeg_HD}.
}
\examples{
MAR_Variance(reptime = 2, type = "V1")

}
