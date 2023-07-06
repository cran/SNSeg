D <- function(S,k,s,e){
  #S: matrix containing partial sums
  #k: potential change point location
  #s: starting point of the segment
  #e: ending point of the segment
  if(e-k <= 1 || k-s+1 <= 1){ # should be k-s+1 <= 1 not k-s <= 1
    return(0)
  }
  else{
    result <- (e-k)*(e-k-1)*S[s,k] + (k-s+1)*(k-s)*S[k+1,e] - (k-s)*(e-k-1)*(S[s,e]-S[s,k]-S[k+1,e])
    return(result)
  }
}

# New divisive algorithm
SN_sweep_mean_HD <- function(data, grid_size, equal_scale=T){
  n <- dim(data)[1]
  p <- dim(data)[2]
  L1 <- matrix(0,n,p)
  L2 <- numeric(n)
  L1[1,] <- data[1,]
  L2[1] <- sum(data[1,]^2)
  for(i in 2:n){
    L1[i,] <- L1[i-1,]+data[i,]
    L2[i] <- L2[i-1]+sum(data[i,]^2)
  }
  S <- matrix(0,n,n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(i==1){
        S[i,j] <- sum(L1[j,]^2)-L2[j]
      }
      else{
        S[i,j] <- sum((L1[j,]-L1[i-1,])^2)-(L2[j]-L2[i-1])
      }
    }
  }

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
        # calculate test statistic for T(t1,k,t2)
        num <- (D(S,k,pre_grid_position,post_grid_position))^2
        denom_L <- denom_R <- 0
        for(t in (pre_grid_position+1):(k-2)){ # Left V
          denom_L <- denom_L + (D(S,t,pre_grid_position,k))^2
        }
        for(t in (k+2):(post_grid_position-2)){ # Right V
          denom_R <- denom_R + (D(S,t,k+1,post_grid_position))^2
        }
        if(equal_scale){
          denom <- (denom_L+denom_R)/(post_grid_position-pre_grid_position+1)
        }else{
          denom <- denom_L/(k-pre_grid_position+1) + denom_R/(post_grid_position-k)
        }
        SN_test_HD_tmp <- num/denom
        sn_grid_stat <- rbind(sn_grid_stat, c(SN_test_HD_tmp, pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}
