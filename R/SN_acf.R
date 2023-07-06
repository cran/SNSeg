SN_sweep_acf <- function(data, grid_size){
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
                              c(SN_test_acf(ts=data[pre_grid_position:post_grid_position], k=k-pre_grid_position+1),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}

SN_test_acf <- function(ts, k){
  n <- length(ts)
  acf1 <- cpp_acf(ts[1:k])
  acf2 <- cpp_acf(ts[(k+1):n])
  D <- k*(n-k)/n^1.5*(acf1-acf2)
  # inter1 <- (1:k)*((k-1):0)*cumsum_acf_constrast_Cpp(ts[1:k],'L')
  # inter2 <- (0:(n-k-1))*((n-k):1)*cumsum_acf_constrast_Cpp(ts[(k+1):n],'R')
  inter1 <- (1:k)*((k-1):0)*cumsum_acf_constrast(ts[1:k],'L')
  inter1[is.na(inter1)] <- 0 # prevent the rare case where x1=x2 at the begining of the series
  inter2 <- (0:(n-k-1))*((n-k):1)*cumsum_acf_constrast(ts[(k+1):n],'R')
  inter2[is.na(inter2)] <- 0 # prevent the rare case where x1=x2 at the begining of the series
  M1 <- sum(inter1^2)/n^2/k^2
  M2 <- sum(inter2^2)/n^2/(n-k)^2
  test_SN <- D^2/(M1+M2)
  return(test_SN)
}
