
SpMAverage = function(x)  {
  #' Spencer moving average 
  #' 
  #' return a 15 period moving average.
  #' @param x numeric vector.
  #' @examples 
  #' x = c(1:20)
  #' SpMAverage(x)
  #' @export SpMAverage
  
  
  wgt = c(-3,-6,-5,3,21,46,67,74,67,46,21,3,-5,-6,-3)   #sum is 320
  
  if (any(is.na(x))) stop("'x' has 'NA'.")
  
  T = length(x)
  if (T<15) stop("length of 'x' should greater than 15.")
  ma = rep(NA,T)
  
  for (t in (8:(T-7)))  {
    ma[t] = sum(x[(t-7):(t+7)] * wgt) / sum(wgt)
  }
  
  # head  //简单化,有几项就对应权重的几项,保证当期前系数为最大的权重74即可.(@LJ)
  for (t in (1:7)) {
    ma[t] = sum(x[1:(t+7)] * wgt[(9-t):length(wgt)]) / sum(wgt[(9-t):length(wgt)])
  }
  # tail
  for (t in (T-6):T) {
    ma[t] = sum(x[(t-7):T] * wgt[1:(T+8-t)]) / sum(wgt[1:(T+8-t)])
  }
  
  ma
}




