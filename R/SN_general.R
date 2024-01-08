SN_test_general <- function(ts,k,functional){
  n <- length(ts)
  g1 <- functional(ts[1:k])
  g2 <- functional(ts[(k+1):n])
  D <- k*(n-k)/n^1.5*(g1-g2)
  M1 <- M2 <- 0
  for(i in 1:k){
    inter <- (i*(k-i)*(functional(ts[1:i])-functional(ts[(i+1):k]))/(n*k))^2
    M1 <- M1+ifelse(is.na(inter), 0, inter)
  }
  for(i in (k+1):n){
    inter <- ((n-i+1)*(i-k-1)*(functional(ts[i:n])-functional(ts[(k+1):(i-1)]))/(n*(n-k)))^2
    M2 <- M2+ifelse(is.na(inter), 0, inter)
  }
  test_SN <- D^2/(M1+M2)
  return(test_SN)
}

SN_sweep_general <- function(data, grid_size, functional){
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
                              c(SN_test_general(ts=data[pre_grid_position:post_grid_position],
                                                k=k-pre_grid_position+1, functional=functional),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}
