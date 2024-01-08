#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @importFrom utils data
#' @importFrom stats approx
#' @importFrom stats acf
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom stats pnorm
NULL

SNSeg_Uni_multi_para <- function(ts, paras_to_test = c(0.9,0.95),
                                confidence = 0.9, grid_size_scale = 0.05,
                                grid_size = NULL, plot_SN = TRUE,
                                est_cp_loc = TRUE){
  if(is.null(ts))
  {stop("Input of ts is missing!")}
  if(!inherits(ts, 'numeric')){
    stop("ts must be numeric!")
  }
  if(!(confidence %in% c(0.9,0.95,0.99,0.995,0.999))){
    stop("Confidence must be one of 0.9, 0.95, 0.99, 0.995 and 0.999!")
  }

  d <- length(paras_to_test)
  paras_to_test_list <- lapply(paras_to_test, function(col) {
    if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
      as.numeric(as.character(col))
    } else {
      col
    }
  })
  for(para in paras_to_test_list){
    if(!is.numeric(para)){
      if(!(para %in% c("mean","variance","acf")))
        stop("Categorical inputs of paras_to_test
             must be within mean, variance or acf!")
    }
    if(is.numeric(para)){
      if((para<=0)|(para>=1))
        stop("Numeric inputs of paras_to_test must be between 0 and 1!")
    }
  }
  n <- length(ts)
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

  SN_sweep_result <- suppressWarnings(SN_sweep_multiparameter(ts, grid_size, paras_to_test))
  SN_result <- SN_divisive_path(start=1, end=n, grid_size, SN_sweep_result, critical_value=critical_value)

  if(plot_SN){
    plot(ts, xlab = "Time", ylab = "Value",
         main="SN Segmentation plot for Multi-Parameters")
    if(est_cp_loc){
      abline(v = SN_result, col = 'red')
    }
  }
  final_result <- structure(
    list(
      ts = ts, paras_to_test = paras_to_test, grid_size = grid_size,
      SN_sweep_result = SN_sweep_result, est_cp = SN_result,
      confidence = confidence, critical_value = critical_value
    ), class = 'SNSeg_Uni'
  )
  final_result
}
