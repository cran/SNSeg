cumsum_mean_constrast <- function(ts, type){
  n <- length(ts)
  result <- rep(0,n)
  if(type=='L'){
    # first sum
    ydm <- ts
    y_bar_f <- cumsum(ydm)/(1:n)
    # second sum
    ydm <- rev(ts)
    y_bar_b <- rev(cumsum(ydm)/(1:n))
    result[1:(n-1)] <- y_bar_f[1:(n-1)]-y_bar_b[2:n]
  }else{
    # first sum
    ydm <- ts
    y_bar_f <- cumsum(ydm)/(1:n)
    # second sum
    ydm <- rev(ts)
    y_bar_b <- rev(cumsum(ydm)/(1:n))
    result[2:n] <- y_bar_b[2:n]-y_bar_f[1:(n-1)]
  }
  return(result)
}

cumsum_acf_constrast <- function(ts, type){
  n <- length(ts)
  result <- rep(0,n)
  if(type=='L'){
    # first sum
    ydm <- ts
    y_bar <- cumsum(ydm)/(1:n)
    acf0_cumsum <- cumsum(ydm^2) - y_bar^2*(1:n)
    acf1_cumsum <- c(0,cumsum(ydm[1:(n-1)]*ydm[2:n]))-c(0,y_bar[-1]*(cumsum(ydm[-n])+cumsum(ydm[-1])-(1:(n-1))*y_bar[-1]))
    acf1_f <- acf1_cumsum/acf0_cumsum
    # second sum
    ydm <- rev(ts)
    y_bar <- cumsum(ydm)/(1:n)
    acf0_cumsum <- cumsum(ydm^2) - y_bar^2*(1:n)
    acf1_cumsum <- c(0,cumsum(ydm[1:(n-1)]*ydm[2:n]))-c(0,y_bar[-1]*(cumsum(ydm[-n])+cumsum(ydm[-1])-(1:(n-1))*y_bar[-1]))
    acf1_b <- rev(acf1_cumsum/acf0_cumsum)
    result[2:(n-2)] <- acf1_f[2:(n-2)]-acf1_b[3:(n-1)]
  }else{
    # first sum
    ydm <- ts
    y_bar <- cumsum(ydm)/(1:n)
    acf0_cumsum <- cumsum(ydm^2) - y_bar^2*(1:n)
    acf1_cumsum <- c(0,cumsum(ydm[1:(n-1)]*ydm[2:n]))-c(0,y_bar[-1]*(cumsum(ydm[-n])+cumsum(ydm[-1])-(1:(n-1))*y_bar[-1]))
    acf1_f <- acf1_cumsum/acf0_cumsum
    # second sum
    ydm <- rev(ts)
    y_bar <- cumsum(ydm)/(1:n)
    acf0_cumsum <- cumsum(ydm^2) - y_bar^2*(1:n)
    acf1_cumsum <- c(0,cumsum(ydm[1:(n-1)]*ydm[2:n]))-c(0,y_bar[-1]*(cumsum(ydm[-n])+cumsum(ydm[-1])-(1:(n-1))*y_bar[-1]))
    acf1_b <- rev(acf1_cumsum/acf0_cumsum)
    result[3:(n-1)] <- acf1_b[3:(n-1)]-acf1_f[2:(n-2)]
  }
  return(result)
}

cumsum_bivcor_constrast <- function(ts, type){
  n <- dim(ts)[2]
  result <- rep(0,n)
  xdm <- ts[1,]; ydm <- ts[2,]
  if(type=='L'){
    # first sum
    x_bar <- cumsum(xdm)/(1:n)
    y_bar <- cumsum(ydm)/(1:n)
    xy <- cumsum(xdm*ydm)
    x2 <- cumsum(xdm^2)
    y2 <- cumsum(ydm^2)
    bivcor_f <- (xy-(1:n)*x_bar*y_bar)/sqrt(x2-(1:n)*x_bar^2)/sqrt(y2-(1:n)*y_bar^2)

    # second sum
    xdm <- rev(xdm)
    ydm <- rev(ydm)
    x_bar <- cumsum(xdm)/(1:n)
    y_bar <- cumsum(ydm)/(1:n)
    xy <- cumsum(xdm*ydm)
    x2 <- cumsum(xdm^2)
    y2 <- cumsum(ydm^2)
    bivcor_b <- rev((xy-(1:n)*x_bar*y_bar)/sqrt(x2-(1:n)*x_bar^2)/sqrt(y2-(1:n)*y_bar^2))

    result[2:(n-2)] <- bivcor_f[2:(n-2)]-bivcor_b[3:(n-1)]
  }else{
    # first sum
    x_bar <- cumsum(xdm)/(1:n)
    y_bar <- cumsum(ydm)/(1:n)
    xy <- cumsum(xdm*ydm)
    x2 <- cumsum(xdm^2)
    y2 <- cumsum(ydm^2)
    bivcor_f <- (xy-(1:n)*x_bar*y_bar)/sqrt(x2-(1:n)*x_bar^2)/sqrt(y2-(1:n)*y_bar^2)

    # second sum
    xdm <- rev(xdm)
    ydm <- rev(ydm)
    x_bar <- cumsum(xdm)/(1:n)
    y_bar <- cumsum(ydm)/(1:n)
    xy <- cumsum(xdm*ydm)
    x2 <- cumsum(xdm^2)
    y2 <- cumsum(ydm^2)
    bivcor_b <- rev((xy-(1:n)*x_bar*y_bar)/sqrt(x2-(1:n)*x_bar^2)/sqrt(y2-(1:n)*y_bar^2))

    result[3:(n-1)] <- bivcor_b[3:(n-1)]-bivcor_f[2:(n-2)]
  }
  result[is.na(result)] <- 0
  return(result)
}

cumsum_variance_constrast <- function(ts, type){
  n <- length(ts)
  result <- rep(0,n)
  if(type=='L'){
    # first sum
    ydm <- ts
    y_bar <- cumsum(ydm)/(1:n)
    acf0_cumsum_f <- cumsum(ydm^2)/(1:n) - y_bar^2
    # second sum
    ydm <- rev(ts)
    y_bar <- cumsum(ydm)/(1:n)
    acf0_cumsum_b <- rev(cumsum(ydm^2)/(1:n) - y_bar^2)
    result[2:(n-2)] <- acf0_cumsum_f[2:(n-2)]-acf0_cumsum_b[3:(n-1)]
  }else{
    # first sum
    ydm <- ts
    y_bar <- cumsum(ydm)/(1:n)
    acf0_cumsum_f <- cumsum(ydm^2)/(1:n) - y_bar^2
    # second sum
    ydm <- rev(ts)
    y_bar <- cumsum(ydm)/(1:n)
    acf0_cumsum_b <- rev(cumsum(ydm^2)/(1:n) - y_bar^2)
    result[3:(n-1)] <- acf0_cumsum_b[3:(n-1)]-acf0_cumsum_f[2:(n-2)]
  }
  return(result)
}
