cmAverage = function(x, frequency = numeric())  {
  #' @title centralized moving average
  #' @param x a number series.
  #' @param frequency get a \code{frequency} period moving average.
  #' @export cmAverage
  #' @examples 
  #' x = c(2,3,4,5,2,7)
  #' cmAverage(x, 3)
  #' x = c(2,3,4,5,2,7)
  #' cmAverage(x,4)

  
  # 中心化移动平均
  # version:1.0.0
  # author:LJ
  
  # x  --数值向量, 不允许有缺失值
  # frequency  --移动平均的项数
  
  
  # error control 
  if (any(is.na(x)))   stop("'x' couldn't have 'NA'.")
  
  T = length(x)
  
  
  # odd ma subfunction
  oddma = function(x, frequency) {
    
    ma = rep(NA, T)  # store ma
    k = (frequency-1)/2
    
    for (t in (k+1):(T-k)) {
      ma[t] = mean(x[(t-k):(t+k)])
    }
    
    # head fill
    for (t in (1:k ))  {
      ma[t] = ((k+2-t)*x[t] + sum(x[1:(t+k)]) - x[t]) / (2*k+1)
    }
    # tail fill
    for (t in ((T-k+1):T)) {
      ma[t] = ((t+1-(T-k))*x[t] + sum(x[(t-k):T]) - x[t]) / (2*k+1)
      
    }
    
    ma
  }
  
  # even ma subfunction
  evenma = function(x, frequency)  {
    ma1 = rep(NA, T)
    ma2 = ma1 # store
    k = frequency / 2
    # step 1
    for (t in k:(T-k))  {
      ma1[t] = mean(x[(t-k+1):(t+k)])
    }
    ## head fill
    for (t in 1:(k-1)) {
      ma1[t] = ((k-t+1)*x[t] + sum(x[1:(t+k)]) - x[t]) / (2*k) 
    }
    ## tail fill
    for (t in (T-k+1):T) {
      ma1[t] = ((t-(T-k)+1)*x[t] + sum(x[(t-k+1):T]) - x[t]) / (2*k)
    }
    
    # step 2 
    for (t in 2:T) {
      ma2[t] = (ma1[t-1]+ma1[t]) / 2
    }
    ma2[1] = ma1[1]
    
    # return ans
    ma2
    
  }
  
  
  # odd ma
  if ((frequency %% 2) != 0 ) {
    ma = oddma(x=x, frequency = frequency)
  }
  # even ma
  else ma = evenma(x=x, frequency = frequency)
  
  
  ma  
  
  
}



