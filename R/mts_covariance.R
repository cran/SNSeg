SN_test_covmatrix_MTS <- function(ts, k){
  # Assuming the mean does not change (and =0), this essentially amounts to mean change detection in (EX^2, EXY, EY^2)
  n <- dim(ts)[2]
  d <- dim(ts)[1]
  covmatrix1 <- covmatrix2 <- c()
  inter1 <- inter2 <- c()
  for(ts_index1 in 1:d){
    for(ts_index2 in ts_index1:d){
      covmatrix1 <- c(covmatrix1, mean(ts[ts_index1,1:k]*ts[ts_index2,1:k]))
      covmatrix2 <- c(covmatrix2, mean(ts[ts_index1,(k+1):n]*ts[ts_index2,(k+1):n]))
      inter1 <- cbind(inter1, cumsum(ts[ts_index1,1:k]*ts[ts_index2,1:k]))
      inter2 <- cbind(inter2, cumsum(ts[ts_index1,(k+1):n]*ts[ts_index2,(k+1):n]))
    }
  }
  D <- k*(n-k)/n^1.5*(covmatrix1-covmatrix2)
  inter1 <- inter1-(1:k)%*%t(covmatrix1)
  inter2 <- inter2-(1:(n-k))%*%t(covmatrix2)

  no_ts <- d*(d+1)/2
  M1 <- M2 <- matrix(0, no_ts, no_ts)
  for(index1 in 1:(no_ts-1)){
    for(index2 in (index1+1):no_ts){
      M1[index1, index2] <- M1[index2, index1] <- sum(inter1[,index1]*inter1[,index2])
      M2[index1, index2] <- M2[index2, index1] <- sum(inter2[,index1]*inter2[,index2])
    }
  }
  for(index1 in 1:no_ts){
    M1[index1, index1] <- sum(inter1[,index1]^2)
    M2[index1, index1] <- sum(inter2[,index1]^2)
  }
  test_SN <- t(D)%*%solve(M1+M2)%*%D*n^2
  return(test_SN)
}

SN_sweep_covmatrix_MTS <- function(data, grid_size){
  n <- dim(data)[2]
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
                              c(SN_test_covmatrix_MTS(ts=data[,pre_grid_position:post_grid_position], k=k-pre_grid_position+1),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}

