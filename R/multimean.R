SN_test_mean_MTS <- function(ts, k){
  no_ts <- dim(ts)[1]
  n <- dim(ts)[2]
  mean1 <- apply(ts[,1:k],1,mean)
  mean2 <- apply(ts[,(k+1):n],1,mean)
  D <- k*(n-k)/n^1.5*(mean1-mean2)
  inter1 <- apply(ts[,1:k],1,cumsum)-(1:k)%*%t(mean1)
  inter2 <- apply(ts[,n:(k+1)],1,cumsum)-(1:(n-k))%*%t(mean2)
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

SN_sweep_mean_MTS <- function(data, grid_size){
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
                              c(SN_test_mean_MTS(ts=data[,pre_grid_position:post_grid_position], k=k-pre_grid_position+1),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}

