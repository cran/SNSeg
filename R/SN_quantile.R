SN_sweep_quantile <- function(data, grid_size, q){
  n <- length(data)
  substat <- list()
  substat[1:(grid_size-1)] <- NA
  for (k in grid_size:(n-grid_size)){
    pre_grid_no <- floor(k/grid_size)
    post_grid_no <- floor((n-k)/grid_size)
    pre_grid_sets <- k-(pre_grid_no:1)*grid_size+1
    post_grid_sets <- k+(1:post_grid_no)*grid_size
    sn_grid_stat <- c()
    for(pre_grid_position in pre_grid_sets){
      for(post_grid_position in post_grid_sets){
        sn_grid_stat <- rbind(sn_grid_stat,
                              c(SN_test_quantile(ts=data[pre_grid_position:post_grid_position], k=k-pre_grid_position+1, q=q),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}

#' @importFrom stats quantile
NULL

SN_test_quantile <- function(ts, k, q){
  n <- length(ts)
  quantile1 <- quantile(ts[1:k],q)
  quantile2 <- quantile(ts[(k+1):n],q)
  D <- k*(n-k)/n^1.5*(quantile1-quantile2)
  inter1 <- (1:k)*((k-1):0)*cumsum_quantile_constrast_Cpp(ts[1:k],'L',q)
  inter2 <- (0:(n-k-1))*((n-k):1)*cumsum_quantile_constrast_Cpp(ts[(k+1):n],'R',q)
  M1 <- sum(inter1^2)/n^2/k^2
  M2 <- sum(inter2^2)/n^2/(n-k)^2
  test_SN <- D^2/(M1+M2)
  return(test_SN)
}
