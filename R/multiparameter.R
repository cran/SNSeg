SN_test_multiparameter <- function(ts, k, paras_to_test){
  n <- length(ts)
  no_para <- length(paras_to_test)
  D <- inter1 <- inter2 <- c()
  for(para in paras_to_test){ # multiparameter test statistics
    if(para=='mean'){
      mean1 <- mean(ts[1:k])
      mean2 <- mean(ts[(k+1):n])
      tmp_D <- k*(n-k)/n^1.5*(mean1-mean2)
      tmp_inter1 <- cumsum_mean_constrast(ts[1:k],'L')
      tmp_inter2 <- cumsum_mean_constrast(ts[(k+1):n],'R')
    }
    if(para=='variance'){
      var1 <- var(ts[1:k])*(k-1)/k
      var2 <- var(ts[(k+1):n])*(n-k-1)/(n-k)
      tmp_D <- k*(n-k)/n^1.5*(var1-var2)
      tmp_inter1 <- cumsum_variance_constrast(ts[1:k],'L')
      tmp_inter2 <- cumsum_variance_constrast(ts[(k+1):n],'R')
    }
    if(para=='acf'){
      acf1 <- cpp_acf(ts[1:k])
      acf2 <- cpp_acf(ts[(k+1):n])
      tmp_D <- k*(n-k)/n^1.5*(acf1-acf2)
      tmp_inter1 <- cumsum_acf_constrast(ts[1:k],'L')
      tmp_inter1[is.na(tmp_inter1)] <- 0 # prevent the rare case where x1=x2 at the begining of the series
      tmp_inter2 <- cumsum_acf_constrast(ts[(k+1):n],'R')
      tmp_inter2[is.na(tmp_inter2)] <- 0 # prevent the rare case where x1=x2 at the begining of the series
    }
    # if(is.numeric(para)){ # quantile
    if(!is.na(as.numeric(para))){
      para <- as.numeric(para)
      q <- para
      quantile1 <- quantile(ts[1:k],q)
      quantile2 <- quantile(ts[(k+1):n],q)
      tmp_D <- k*(n-k)/n^1.5*(quantile1-quantile2)
      tmp_inter1 <- cumsum_quantile_constrast_Cpp(ts[1:k],'L',q)
      tmp_inter2 <- cumsum_quantile_constrast_Cpp(ts[(k+1):n],'R',q)
    }
    D <- c(D, tmp_D)
    inter1 <- cbind(inter1, tmp_inter1)
    inter2 <- cbind(inter2, tmp_inter2)
  }

  multiplier1 <- ((1:k)*((k-1):0))^2/n^2/k^2
  multiplier2 <- ((0:(n-k-1))*((n-k):1))^2/n^2/(n-k)^2
  M1 <- M2 <- matrix(0, no_para, no_para)
  for(index1 in 1:(no_para-1)){
    for(index2 in (index1+1):no_para){
      M1[index1, index2] <- M1[index2, index1] <- sum(inter1[,index1]*inter1[,index2]*multiplier1)
      M2[index1, index2] <- M2[index2, index1] <- sum(inter2[,index1]*inter2[,index2]*multiplier2)
    }
  }
  for(index1 in 1:no_para){
    M1[index1, index1] <- sum(inter1[,index1]^2*multiplier1)
    M2[index1, index1] <- sum(inter2[,index1]^2*multiplier2)
  }
  test_SN <- t(D)%*%solve(M1+M2)%*%D
  return(test_SN)
}

SN_sweep_multiparameter <- function(data, grid_size, paras_to_test){
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
                              c(SN_test_multiparameter(ts=data[pre_grid_position:post_grid_position], k=k-pre_grid_position+1, paras_to_test=paras_to_test),
                                pre_grid_position, k, post_grid_position))
      }
    }
    substat[[k]] <- sn_grid_stat
  }
  substat[(n-grid_size+1):n] <- NA
  return(substat)
}
