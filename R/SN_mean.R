SN_sweep_mean <- function(data, grid_size){
  n <- length(data)
  substat <- list()
  substat[1:(grid_size-1)] <- NA
  for (k in grid_size:(n-grid_size)) {
    pre_grid_no <- floor(k/grid_size)
    post_grid_no <- floor((n-k)/grid_size)
    pre_grid_sets <- k-(pre_grid_no:1)*grid_size+1
    post_grid_sets <- k+(1:post_grid_no)*grid_size
    sn_grid_stat <- c()
    for(pre_grid_position in pre_grid_sets){
      for(post_grid_position in post_grid_sets){
        sn_grid_stat <- rbind(sn_grid_stat,
                              c(SN_test(ts=data[pre_grid_position:post_grid_position], k=k-pre_grid_position+1),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}

SN_test <- function(ts, k){
n <- length(ts)
mean1 <- mean(ts[1:k])
mean2 <- mean(ts[(k+1):n])
inter1 <- cumsum(ts[1:k])-(1:k)*mean1
inter2 <- cumsum(ts[n:(k+1)])-(1:(n-k))*mean2
M1 <- sum(inter1^2)/n^2
M2 <- sum(inter2^2)/n^2
test_SN <- n*(((n-k)*k/n^2)*(mean1-mean2))^2/(M1+M2)
return(test_SN)
}
