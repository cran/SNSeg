SN_divisive_path <- function(start, end, grid_size, SN_sweep_result, critical_value){
  n <- end-start+1
  if(n<(2*grid_size)){
    return(NULL)
  }
  substat <- rep(0,n)
  for (k in grid_size:(n-grid_size)){
    actual_position <- k+start-1
    sn_grid_stat <- SN_sweep_result[[actual_position]]
    sn_grid_stat <- sn_grid_stat[sn_grid_stat[,2]>=start,]
    if(is.null(dim(sn_grid_stat))){
      sn_grid_stat <- sn_grid_stat[sn_grid_stat[4]<=end]
    }else{
      sn_grid_stat <- sn_grid_stat[sn_grid_stat[,4]<=end,]
    }
    if(is.null(dim(sn_grid_stat))){
      substat[k] <- sn_grid_stat[1]
    }else{
      substat[k] <- max(sn_grid_stat[,1])
    }
  }
  if(sum(substat>critical_value)==0){
    return(NULL)
  }else{
    est_cp <- which.max(substat)+start-1
    return(c(SN_divisive_path(start,est_cp,grid_size,SN_sweep_result,critical_value),
             est_cp,
             SN_divisive_path(est_cp+1,end,grid_size,SN_sweep_result,critical_value)))
  }
}
