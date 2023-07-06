#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @importFrom utils data
#' @importFrom stats approx
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom mvtnorm rmvnorm
NULL

SNSeg_Uni_single_para <- function(ts, type = "mean", confidence = 0.9,
                      grid_size_scale = 0.05, grid_size = NULL,
                      plot_SN = TRUE, est_cp_loc = TRUE){
  if(is.null(ts))
    {stop("Input of ts is missing!")}
  if(!((type %in% c("mean","variance","acf","bivcor")) | ((type>=0)&(type<=1)) )){
    stop("paras_to_test must come from one of the categories: mean, variance, acf,
         bivcor, or a percentage between 0 and 1!")}
  if(!(confidence %in% c(0.9,0.95,0.99,0.995,0.999))){
    stop("Confidence must be one of 0.9, 0.95, 0.99, 0.995 and 0.999!")
  }

  if((type %in% c("mean","variance","acf")) | (inherits(type,'numeric'))){
    if(!inherits(ts, 'numeric')){
      stop("ts must be numeric!")
    }
    n <- length(ts)
  }
  if(type == "bivcor"){
    if(!inherits(ts,'matrix')){
      stop("ts must be a matrix!")
    }
    if(!((dim(ts)[1] == 2) | (dim(ts)[2] == 2))){
      stop("ts must be 2-dimensional!")
    }

    if(dim(ts)[2] == 2) {ts <- t(ts)}
    n <- dim(ts)[2]
  }

  critical_values_single <- SNSeg::critical_values_single
  # interpolation for grid_size_scale
  if((is.null(grid_size))&(is.null(grid_size_scale))){
    grid_size_scale <- 0.05
    posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
    critical_value <- critical_values_single[1,posi_confidence]
    grid_size <- floor(grid_size_scale*n)
    warning("Undefined value detected for both grid_size and grid_size_scale! The system would use 0.05 as the default for grid_size_scale.")
  } else if((is.null(grid_size))&(!is.null(grid_size_scale))){
    if(grid_size_scale<0.05){
      grid_size_scale <- 0.05
      posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
      critical_value <- critical_values_single[1,posi_confidence]
      grid_size <- floor(grid_size_scale*n)
      warning("Detected the grid_size_scale is less than 0.05. The system would use 0.05 for grid_size_scale.")
    }
    if(grid_size_scale>0.5){
      grid_size_scale <- 0.5
      posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
      critical_value <- critical_values_single[18,posi_confidence]
      grid_size <- floor(grid_size_scale*n)
      warning("Detected the grid_size_scale is greater than 0.5. The system would use 0.5 for grid_size_scale.")
    }
    if(grid_size_scale>=0.05 & grid_size_scale<=0.5){
      if((grid_size_scale %in% critical_values_single[,1])){
        posi_cri <- match(as.character(grid_size_scale),critical_values_single[,1])
        posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
        critical_value <- critical_values_single[posi_cri,posi_confidence]
        grid_size <- floor(grid_size_scale*n)
      } else if(!(grid_size_scale %in% critical_values_single[,1])){
        grid_size <- floor(grid_size_scale*n)
        posi_epsilon <- match(grid_size_scale,sort(c(critical_values_single[,1],grid_size_scale)))
        posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
        critical_vector <- c(critical_values_single[(posi_epsilon-1):posi_epsilon,posi_confidence])
        epsilon_vector <- c(critical_values_single[(posi_epsilon-1):posi_epsilon,1])
        critical_value <- approx(epsilon_vector,critical_vector,xout = grid_size_scale)$y
      }
    }
  } else if(!is.null(grid_size)){
    grid_size_scale <- grid_size/n
    if(grid_size_scale<0.05){
      grid_size_scale <- 0.05
      grid_size <- round(grid_size_scale*n)
      posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
      critical_value <- critical_values_single[1,posi_confidence]
      warning("Detected the grid_size_scale is less than 0.05 from the current grid_size. The system would use 0.05 for grid_size_scale.")
    }
    if(grid_size_scale>0.5){
      grid_size_scale <- 0.5
      grid_size <- round(grid_size_scale*n)
      posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
      critical_value <- critical_values_single[18,posi_confidence]
      warning("Detected the grid_size_scale is greater than 0.5 from the current grid_size. The system would use 0.5 for grid_size_scale.")
    }
    if(grid_size_scale>=0.05 & grid_size_scale<=0.5){
      if((grid_size_scale %in% critical_values_single[,1])){
        posi_cri <- match(as.character(grid_size_scale),critical_values_single[,1])
        posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
        critical_value <- critical_values_single[posi_cri,posi_confidence]
      } else if(!(grid_size_scale %in% critical_values_single[,1])){
        posi_epsilon <- match(grid_size_scale,sort(c(critical_values_single[,1],grid_size_scale)))
        posi_confidence <- match(as.character(confidence),colnames(critical_values_single))
        critical_vector <- c(critical_values_single[(posi_epsilon-1):posi_epsilon,posi_confidence])
        epsilon_vector <- c(critical_values_single[(posi_epsilon-1):posi_epsilon,1])
        critical_value <- approx(epsilon_vector,critical_vector,xout = grid_size_scale)$y
      }
    }
  }

  # Mean
  if(type == "mean"){

    # SN change points estimate
    SN_sweep_result <- SN_sweep_mean(ts, grid_size)
    SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)

    if(plot_SN){
      plot(ts, xlab = "Time", ylab = "Value",
           main="SN Segmentation plot for Univariate Mean")
      if(est_cp_loc){
        abline(v = SN_result, col = 'red')
      }
    }
  }

  # Variance
  else if (type == "variance"){

    # SN change points estimate
    SN_sweep_result <- SN_sweep_variance(ts, grid_size)
    SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)

    if(plot_SN){
      plot(ts, xlab = "Time", ylab = "Value",
           main="SN Segmentation plot for Univariate Variance")
      if(est_cp_loc){
        abline(v = SN_result, col = 'red')
      }
    }
  }

  # ACF
  else if(type == "acf"){

    # SN change points estimate
    SN_sweep_result <- SN_sweep_acf(ts, grid_size)
    SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)

    if(plot_SN){
      plot(ts, xlab = "Time", ylab = "Value",
           main="SN Segmentation plot for Univariate ACF")
      if(est_cp_loc){
        abline(v = SN_result, col = 'red')
      }
    }
  }

  # bivariate correlation
  else if(type == "bivcor"){

    # SN change points estimate
    SN_sweep_result <- SN_sweep_bivcor(ts, grid_size)
    SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)

    if(plot_SN){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      par(mfrow=c(1,2))
      plot(ts[1,], xlab = "Time", ylab = "Value",
           main="SN Segmentation plot for Bivariate Correlation")
      if(est_cp_loc){
        abline(v = SN_result, col = 'red')
      }
      plot(ts[2,], xlab = "Time", ylab = "Value",
           main="SN Segmentation plot for Bivariate Correlation")
      if(est_cp_loc){
        abline(v = SN_result, col = 'red')
      }
    }
  }

  # Quantile
  else if(inherits(type, 'numeric')){
    quantile_level <- type
    type <- paste(quantile_level*100, '% quantile', sep = '')

    # SN change points estimate
    SN_sweep_result <- SN_sweep_quantile(ts, grid_size, quantile_level)
    SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)

    if(plot_SN){
      plot(ts, xlab = "Time", ylab = "Value",
           main=paste0("SN Segmentation plot for ",quantile_level*100,"th Quantile"))
      if(est_cp_loc){
        abline(v = SN_result, col = 'red')
      }
    }
  }
  return(list("paras_to_test" = type, "grid_size" = grid_size,
              "SN_sweep_result" = SN_sweep_result, "est_cp" = SN_result,
              "confidence" = confidence, "critical_value" = critical_value))
}

