#' @useDynLib SNSeg, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @importFrom stats acf
NULL

#' Parameter estimates of each segment separated by Self-Normalization (SN) based
#' change-point estimates
#'
#' The function \code{SNSeg_estimate} computes parameter estimates of each segment
#' that are separated by the SN-based change-point estimates.
#'
#' @param SN_result An S3 object served as the output of the functions \code{SNSeg_Uni},
#' \code{SNSeg_Multi}, or \code{SNSeg_HD}.
#'
#' @return \code{SNSeg_estimate} returns an S3 object of class "SNSeg_estimate" including
#' the parameter estimates of each segment separated by the SN-based change-point
#' estimates.
#'
#' \enumerate{
#'   \item If the time series is univariate, for a single parameter change, the output
#'   contains parameter estimates for one of the followings: \code{mean}, \code{variance},
#'   \code{acf}, \code{quantile}, or \code{general}, which can be referred to the change
#'   in a single mean, variance, autocorrelation, a given quantile level, or a general
#'   functional. For multi-parameter changes, the output can be a combination of
#'   \code{mean}, \code{variance}, \code{acf}, and a dataframe with each quantile level
#'   depending on the type of parameters (argument \code{paras_to_test} of \code{SNSeg_Uni},
#'   \code{SNSeg_Multi}, or \code{SNSeg_HD}) that users select.
#'
#'   \item If the time series is multivariate with a dimension no greater than 10, the output
#'   contains parameter estimates for one of the followings: \code{bivcor}, \code{multi_mean},
#'   or \code{covariance}, which can be referred to the change in correlation between
#'   bivariate time series and the change in multivariate means or covariance between
#'   multivariate time series.
#'
#'   \item If the time series is high-dimensional with a dimension greater than 10, the output
#'   contains the parameter estimate \code{HD_mean} to represent the change in high-dimensional
#'   means.
#' }
#'
#' For more examples of \code{SNSeg_estimate} see the help vignette:
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
#' # variance estimates of the separated segments
#' SNSeg_estimate(SN_result = result)
#'
#' # For more examples of SNSeg_estimate, please run
#' # the command: vignette("SNSeg", package = "SNSeg")
#' }
#'
#' @export
SNSeg_estimate <- function(SN_result){
  # SN_result must be a S3 class oblect
  ts <- SN_result$ts
  cpt_loc <- c(1,SN_result$est_cp,max(dim(as.matrix(ts))))
  para <- SN_result$paras_to_test
  if(inherits(ts, 'matrix')){
    if(dim(ts)[1]>dim(ts)[2]){ts <- t(ts)}
    # check dimension
    if(!is.null(para)){
      if(para == 'bivcor'){ # bivariate correlation
        SN_estimate <- list()
        SN_estimate$bivcor <- c()
        for(i in 1:(length(cpt_loc)-1)){
          ts_segment1 <- ts[1,cpt_loc[i]:cpt_loc[i+1]]
          ts_segment2 <- ts[2,cpt_loc[i]:cpt_loc[i+1]]
          SN_estimate$bivcor <- c(SN_estimate$bivcor, cor(ts_segment1,ts_segment2))
        }
      }else if(para == 'covariance'){
        SN_estimate <- list()
        matrix_cov <- list()
        for(i in 1:(length(cpt_loc)-1)){
          ts_segment <- ts[,cpt_loc[i]:cpt_loc[i+1]]
          matrix_cov[[i]] <- cov(t(ts_segment))
        }
        SN_estimate$covariance <- matrix_cov
      }else{ # multi-mean
        SN_estimate <- list()
        matrix_mean <- matrix(0,nrow=length(cpt_loc)-1,ncol=dim(ts)[1])
        for(i in 1:(length(cpt_loc)-1)){
          ts_segment <- ts[,cpt_loc[i]:cpt_loc[i+1]]
          matrix_mean[i,] <- rowMeans(ts_segment)
        }
        matrix_mean <- data.frame(matrix_mean)
        colnames(matrix_mean) <- paste0('ts',1:dim(ts)[1])
        SN_estimate$multi_mean <- matrix_mean
      }
    }else{ # high-dimensional mean
      SN_estimate <- list()
      matrix_mean <- matrix(0,nrow=length(cpt_loc)-1,ncol=dim(ts)[2])
      for(i in 1:(length(cpt_loc)-1)){
        ts_segment <- ts[,cpt_loc[i]:cpt_loc[i+1]]
        matrix_mean[i,] <- rowMeans(ts_segment)
      }
      matrix_mean <- data.frame(matrix_mean)
      colnames(matrix_mean) <- paste0('ts',1:dim(ts)[2])
      SN_estimate$HD_mean <- matrix_mean
    }
  } else{
    SN_estimate <- list()

    if(length(para) == 1){
      if(!inherits(para, 'function')){
        if(para == 'mean'){
          SN_estimate$mean <- c()
          for(i in 1:(length(cpt_loc)-1)){
            ts_segment <- ts[cpt_loc[i]:cpt_loc[i+1]]
            SN_estimate$mean <- c(SN_estimate$mean, mean(ts_segment))
          }
        } else if(para == 'variance'){
          SN_estimate$variance <- c()
          for(i in 1:(length(cpt_loc)-1)){
            ts_segment <- ts[cpt_loc[i]:cpt_loc[i+1]]
            SN_estimate$variance <- c(SN_estimate$variance, var(ts_segment))
          }
        } else if(para == 'acf'){
          SN_estimate$acf <- c()
          for(i in 1:(length(cpt_loc)-1)){
            ts_segment <- ts[cpt_loc[i]:cpt_loc[i+1]]
            SN_estimate$acf <- c(SN_estimate$acf, as.vector(acf(ts_segment,lg.max=1,plot=F)$acf)[2])
          }
        } else if(inherits(para, 'numeric')){ # quantile
          SN_estimate$quantile <- c()
          for(i in 1:(length(cpt_loc)-1)){
            ts_segment <- ts[cpt_loc[i]:cpt_loc[i+1]] # return the one that user is interested in
            SN_estimate$quantile <- c(SN_estimate$quantile, as.numeric(quantile(ts_segment,para)))
          }
        }
      }else{
        # general functional
        SN_estimate$general <- c()
        for(i in 1:(length(cpt_loc)-1)){
          ts_segment <- ts[cpt_loc[i]:cpt_loc[i+1]] # return the one that user is interested in
          SN_estimate$general <- c(SN_estimate$general, as.numeric(para(ts_segment)))
        }
      }

    }else{
      para_list <- suppressWarnings(as.numeric(para))
      para_list <- para_list[!is.na(para_list)]
      if(length(para_list != 0)){
        qt_matrix <- matrix(0,nrow=length(cpt_loc)-1,ncol=length(para_list))
        for(i in 1:(length(cpt_loc)-1)){
          ts_segment <- ts[cpt_loc[i]:cpt_loc[i+1]]
          if('mean' %in% para){SN_estimate$mean <- c(SN_estimate$mean, mean(ts_segment))}
          if('variance' %in% para){SN_estimate$variance <- c(SN_estimate$variance, var(ts_segment))}
          if('acf' %in% para){SN_estimate$acf <- c(SN_estimate$acf, as.vector(acf(ts_segment,lg.max=1,plot=F)$acf)[2])}
          qt_matrix[i,] <- as.numeric(quantile(ts_segment,para_list))
        }
        qt_matrix <- data.frame(qt_matrix)
        colnames(qt_matrix) <- para_list
        SN_estimate$quantile <- qt_matrix
      } else{
        for(i in 1:(length(cpt_loc)-1)){
          ts_segment <- ts[cpt_loc[i]:cpt_loc[i+1]]
          if('mean' %in% para){SN_estimate$mean <- c(SN_estimate$mean, mean(ts_segment))}
          if('variance' %in% para){SN_estimate$variance <- c(SN_estimate$variance, var(ts_segment))}
          if('acf' %in% para){SN_estimate$acf <- c(SN_estimate$acf, as.vector(acf(ts_segment,lg.max=1,plot=F)$acf)[2])}
        }
      }
    }
  }
  final_result <- structure(
    SN_estimate, class = 'SNSeg_estimate'
  )
  final_result
}

