SN_sweep_variance <- function(data, grid_size){
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
                              c(SN_test_variance(ts=data[pre_grid_position:post_grid_position], k=k-pre_grid_position+1),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}

#' @importFrom stats var
NULL

SN_test_variance <- function(ts, k){
  n <- length(ts)
  var1 <- var(ts[1:k])*(k-1)/k
  var2 <- var(ts[(k+1):n])*(n-k-1)/(n-k)
  D <- k*(n-k)/n^1.5*(var1-var2)
  inter1 <- (1:k)*((k-1):0)*cumsum_variance_constrast(ts[1:k],'L')
  inter2 <- (0:(n-k-1))*((n-k):1)*cumsum_variance_constrast(ts[(k+1):n],'R')
  M1 <- sum(inter1^2)/n^2/k^2
  M2 <- sum(inter2^2)/n^2/(n-k)^2
  test_SN <- D^2/(M1+M2)
  return(test_SN)
}
