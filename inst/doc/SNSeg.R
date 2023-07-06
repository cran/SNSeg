## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(SNSeg)

## -----------------------------------------------------------------------------
# Please run the following function before running examples:
mix_GauGPD <- function(u,p,trunc_r,gpd_scale,gpd_shape){
  indicator <- u<p
  rv <- rep(0, length(u))
  rv[indicator>0] <- qtruncnorm(u[indicator>0]/p,a=-Inf,b=trunc_r)
  rv[indicator<=0] <- qgpd((u[indicator<=0]-p)/(1-p), loc=trunc_r, scale=gpd_scale,shape=gpd_shape)
  return(rv)
}

## -----------------------------------------------------------------------------
set.seed(7)
n <- 2000
reptime <- 2
cp_sets <- round(n*c(0,cumsum(c(0.5,0.25)),1))
mean_shift <- c(0.4,0,0.4)
rho <- -0.7
ts <- MAR(n, reptime, rho)
no_seg <- length(cp_sets)-1
for(index in 1:no_seg){ # Mean shift
  tau1 <- cp_sets[index]+1   
  tau2 <- cp_sets[index+1]
  ts[tau1:tau2,] <- ts[tau1:tau2,] + mean_shift[index]
}
ts <- ts[,2]
# grid_size undefined
result <- SNSeg_Uni(ts, paras_to_test = "mean", confidence = 0.9,
                    grid_size_scale = 0.05, grid_size = NULL, 
                    plot_SN = FALSE, est_cp_loc = FALSE)
# grid_size defined & generate time series segmentation plot
result <- SNSeg_Uni(ts, paras_to_test = "mean", confidence = 0.9,
                    grid_size_scale = 0.05, grid_size = 116, 
                    plot_SN = TRUE, est_cp_loc = TRUE)
# Estimated change-point locations
result$est_cp
# plot the SN-based test statistic
SN_stat <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

## -----------------------------------------------------------------------------
set.seed(7)
ts <- MAR_Variance(2, "V1")
ts <- ts[,2]
# grid_size defined
result <- SNSeg_Uni(ts, paras_to_test = "variance", confidence = 0.9,
                    grid_size_scale = 0.05, grid_size = 67, 
                    plot_SN = TRUE, est_cp_loc = TRUE)
# Estimated change-point locations
result$est_cp
# plot the SN-based test statistic
SN_stat <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

## -----------------------------------------------------------------------------
set.seed(7)
ts <- MAR_Variance(2, "A3")
ts <- ts[,2]
# grid_size defined
result <- SNSeg_Uni(ts, paras_to_test = "acf", confidence = 0.9,
          grid_size_scale = 0.05, grid_size = 92, plot_SN = TRUE,
          est_cp_loc = TRUE)
# Estimated change-point locations
result$est_cp
# plot the SN-based test statistic
SN_stat <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

## -----------------------------------------------------------------------------
library(mvtnorm)
set.seed(7)
n <- 1000
sigma_cross <- list(4*matrix(c(1,0.8,0.8,1), nrow=2),
                    matrix(c(1,0.2,0.2,1), nrow=2),
                    matrix(c(1,0.8,0.8,1), nrow=2))
cp_sets <- round(c(0,n/3,2*n/3,n))
noCP <- length(cp_sets)-2
rho_sets <- rep(0.5, noCP+1)
ts <- MAR_MTS_Covariance(n, 2, rho_sets, cp_sets, sigma_cross)
ts <- ts[1][[1]]
# grid_size defined
result <- SNSeg_Uni(ts, paras_to_test = "bivcor", confidence = 0.9,
                    grid_size_scale = 0.05, grid_size = 77, 
                    plot_SN = TRUE, est_cp_loc = TRUE)
# Estimated change-point locations
result$est_cp
# plot the SN-based test statistic
SN_stat <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

## -----------------------------------------------------------------------------
library(truncnorm)
library(evd)
set.seed(7)
n <- 1000
cp_sets <- c(0,n/2,n)
noCP <- length(cp_sets)-2
reptime <- 2
rho <- 0.2
# AR time series with no change-point (mean, var)=(0,1)
ts <- MAR(n, reptime, rho)*sqrt(1-rho^2)
trunc_r <- 0
p <- pnorm(trunc_r)
gpd_scale <- 2
gpd_shape <- 0.125
for(ts_index in 1:reptime){
  ts[(cp_sets[2]+1):n, ts_index] <- mix_GauGPD(pnorm(ts[(cp_sets[2]+1):n, ts_index]),
p,trunc_r,gpd_scale,gpd_shape)
}
ts <- ts[,2]
# grid_size undefined
# test in 90% quantile
result <- SNSeg_Uni(ts, paras_to_test = 0.9, confidence = 0.9,
                    grid_size_scale = 0.066, grid_size = NULL,
                    plot_SN = TRUE, est_cp_loc = TRUE)
# Estimated change-point locations
result$est_cp
# plot the SN-based test statistic
SN_stat <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

## -----------------------------------------------------------------------------
set.seed(7)
n <- 1000
cp_sets <- c(0,333,667,1000)
no_seg <- length(cp_sets)-1
rho <- 0
# AR time series with no change-point (mean, var)=(0,1)
ts <- MAR(n, 2, rho)*sqrt(1-rho^2)
no_seg <- length(cp_sets)-1
sd_shift <- c(1,1.6,1)
for(index in 1:no_seg){ # Mean shift
  tau1 <- cp_sets[index]+1
  tau2 <- cp_sets[index+1]
  ts[tau1:tau2,] <- ts[tau1:tau2,]*sd_shift[index]
}
d <- 2
ts <- ts[,2]

# Test in 90th and 95th quantile with grid_size undefined
result <- SNSeg_Uni(ts, paras_to_test = c(0.9, 0.95), confidence = 0.9,
                    grid_size_scale = 0.05, grid_size = NULL, 
                    plot_SN = FALSE, est_cp_loc = FALSE)

# Test in 90th quantile and the variance with grid_size undefined
result <- SNSeg_Uni(ts, paras_to_test = c(0.9, 'variance'),
                    confidence = 0.95, grid_size_scale = 0.078,
                    grid_size = NULL, plot_SN = FALSE, 
                    est_cp_loc = FALSE)

# Test in 90th quantile, variance and acf with grid_size undefined
result <- SNSeg_Uni(ts, paras_to_test = c(0.9,'variance', "acf"),
                    confidence = 0.9, grid_size_scale = 0.064,
                    grid_size = NULL, plot_SN = TRUE, 
                    est_cp_loc = TRUE)
# Estimated change-point locations
result$est_cp
# Test in 60th quantile, mean, variance and acf with grid_size defined
result <- SNSeg_Uni(ts, paras_to_test = c(0.6, 'mean', "variance",   
                    "acf"), confidence = 0.9, grid_size_scale = 0.05,
                    grid_size = 83, plot_SN = TRUE, est_cp_loc = TRUE)
# Estimated change-point locations
result$est_cp
# SN-based test statistic segmentation plot
SN_stst <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

## -----------------------------------------------------------------------------
# Please run this function before running examples:
exchange_cor_matrix <- function(d, rho){
  tmp <- matrix(rho, d, d)
  diag(tmp) <- 1
  return(tmp)
}

## -----------------------------------------------------------------------------
library(mvtnorm)
set.seed(10)
d <- 5
n <- 1000
cp_sets <- round(n*c(0,cumsum(c(0.075,0.3,0.05,0.1,0.05)),1))
mean_shift <- c(-3,0,3,0,-3,0)/sqrt(d)
mean_shift <- sign(mean_shift)*ceiling(abs(mean_shift)*10)/10
rho_sets <- 0.5
sigma_cross <- list(exchange_cor_matrix(d,0))
ts <- MAR_MTS_Covariance(n, 2, rho_sets, cp_sets=c(0,n), sigma_cross)
noCP <- length(cp_sets)-2
no_seg <- length(cp_sets)-1
for(rep_index in 1:2){
  for(index in 1:no_seg){ # Mean shift
    tau1 <- cp_sets[index]+1
    tau2 <- cp_sets[index+1]
    ts[[rep_index]][,tau1:tau2] <- ts[[rep_index]][,tau1:tau2] + mean_shift[index]
    }
}
ts <- ts[1][[1]]

# grid_size undefined
result <- SNSeg_Multi(ts, paras_to_test = "mean", confidence = 0.95,
                      grid_size_scale = 0.079, grid_size = NULL)
# grid_size defined
result <- SNSeg_Multi(ts, paras_to_test = "mean", confidence = 0.99,
                      grid_size_scale = 0.05, grid_size = 65)
# Estimated change-point locations
result$est_cp
# SN-based test statistic segmentation plot
SN_stst <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

## -----------------------------------------------------------------------------
library(mvtnorm)
set.seed(10)
reptime <- 2
d <- 4
n <- 1000
sigma_cross <- list(exchange_cor_matrix(d,0.2),
                    2*exchange_cor_matrix(d,0.5),
                    4*exchange_cor_matrix(d,0.5))
rho_sets <- c(0.3,0.3,0.3)
mean_shift <- c(0,0,0) # with mean change
cp_sets <- round(c(0,n/3,2*n/3,n))
ts <- MAR_MTS_Covariance(n, reptime, rho_sets, cp_sets, sigma_cross)
noCP <- length(cp_sets)-2
no_seg <- length(cp_sets)-1
for(rep_index in 1:reptime){
  for(index in 1:no_seg){ # Mean shift
    tau1 <- cp_sets[index]+1
    tau2 <- cp_sets[index+1]
    ts[[rep_index]][,tau1:tau2] <- ts[[rep_index]][,tau1:tau2] +
      mean_shift[index]
    }
}
ts <- ts[[1]]

# grid_size undefined
result <- SNSeg_Multi(ts, paras_to_test = "covariance", 
                      confidence = 0.9, grid_size_scale = 0.05, 
                      grid_size = NULL)
# grid_size defined
result <- SNSeg_Multi(ts, paras_to_test = "covariance", 
                      confidence = 0.9, grid_size_scale = 0.05, 
                      grid_size = 81)
# Estimated change-point locations
result$est_cp
# SN-based test statistic segmentation plot
SN_stst <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

## -----------------------------------------------------------------------------
n <- 600
p <- 100
nocp <- 5
cp_sets <- round(seq(0,nocp+1,1)/(nocp+1)*n)
num_entry <- 5
kappa <- sqrt(4/5) # Wang et al(2020)
mean_shift <- rep(c(0,kappa),100)[1:(length(cp_sets)-1)]
set.seed(1)
ts <- matrix(rnorm(n*p,0,1),n,p)
no_seg <- length(cp_sets)-1
for(index in 1:no_seg){ # Mean shift
  tau1 <- cp_sets[index]+1
  tau2 <- cp_sets[index+1]
  ts[tau1:tau2,1:num_entry] <- ts[tau1:tau2,1:num_entry] +
    mean_shift[index] # sparse change
}

# grid_size undefined
result <- SNSeg_HD(ts, confidence = 0.9, grid_size_scale = 0.05,
                   grid_size = NULL)
# grid_size defined
result <- SNSeg_HD(ts, confidence = 0.9, grid_size_scale  = 0.05,
                   grid_size = 52)
# Estimated change-point locations
result$est_cp
# SN-based test statistic segmentation plot
SN_stst <- max_SNsweep(result, plot_SN = TRUE, est_cp_loc = TRUE,
                       critical_loc = TRUE)

