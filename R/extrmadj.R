#' adjuste extreme values  
#' @description 
#' Adjusting extreme values of series.
#' @param x numeric vector.
#' @param limit threshold for recognizing extreme value.
#' @param subs_method indicate method for substituting extreme values, \code{spencer} or \code{nearest}, 
#' the former means substituted by \emph{spencer} value, latter means substituted by mean of nearest 
#' none extreme values. 
#' @export extrmadj
#' @examples 
#' 
#' ##
#' usad1 = USAccDeaths
#' usad1[2:5] = mean(USAccDeaths) * 3  # first obs will become a extreme since the trend has been changed by 2:5
#' (ext1 = extrmadj(usad1))
#' plot(ext1$extr_infr$x, type = "l") # original
#' lines(ext1$x.nextr, type = "l", col = "red") # substituted
#' ## 
#' usad2 = USAccDeaths
#' usad2[15] = - mean(USAccDeaths) * 3
#' (ext2 = extrmadj(usad2))
#' plot(ext2$extr_infr$x, type = "l") # original
#' lines(ext2$x.nextr, type = "l", col = "red") # substituted
#' 
#' ## substituted by mean of nearest non-extrmadj
#' (ext3 = extrmadj(usad2, subs_method = "near"))
#' plot(ext3$extr_infr$x, type = "l") # original
#' lines(ext3$x.nextr, type = "l", col = "red") # substituted
#' 
# x:向量,季节调整后的原始值序列
# limit:界限, 默认为3.5,级表示超过3.5倍的标准差为异常值.
# 
#  列表: 包含一个数据框/一个向量/阈值. 其中数据框含有原始序列/是否异常值/调整后序列三列; 向量即为调整后序列; 阈值即表示波动超过它则认为是异常值.
# x.nextr : 极端值修正后的x序列
# x : 原序列
# x.isextreme : 标识原序列中该行是否为极端值, TRUE:是极端值, FALSE:不是极端值.
# x.isextrextr : 极端异常值, 迭代1000次仍无法修正的.
# limits : 判断上下界限(比例的), 最近一次迭代的结果.
# iteration : 迭代次数


extrmadj = function(x, limit = 3.5, subs_method = c("spencer", "nearest"))  {
  
  if (any(is.na(x))) stop("'x' has 'NA'.")
  
  T = length(x)
  
  if (T<4)  stop("'x' too short.")
  
  # subfunction : substituted by Spencer values #
  extreme_sp = function(x, limit = 3.5) {
    x.nextr = x
    sp = SpMAverage(x) #得到Spencer序列
    ratio = x.nextr / sp  #ratio of x and Spencer series.
    sd = sd(ratio)
    mn = mean(ratio)
    threshold_ceil = mn + limit * sd # unit value \ ±3.5个标准差
    threshold_floor = mn - limit * sd 
    
    isextreme = ((ratio > threshold_ceil) | (ratio < threshold_floor))   # boole values express whether it is extreme.
    ext_dx = which(isextreme == TRUE)  # store TRUE 's index in original series.
    l_dx = length(ext_dx)
    
    
    if (l_dx>0L) { # substitute while there are extreme values.
      for (i in (1:l_dx))  {
        t = ext_dx[i] # 定位到原序列中异常值的索引
        x.nextr[t] = sp[t]
        
      }
      
    }
    
    limits = c(threshold_floor, threshold_ceil )  # 上下界限
    extr_infr = data.frame(x=x, x.isextreme=isextreme, x.nextr = x.nextr)
    lst = list(x.nextr = x.nextr, extr_infr = extr_infr, limits = limits)
    lst
    
  }
  
  
  # sub-function : substituted by nearest #
  extreme_nrst = function(x, limit = 3.5)  {
    
    
    # unit value: T,l_dx, k, threshold
    
    
    
    sp = SpMAverage(x) #得到Spencer序列
    # 
    # T = length(x)
    # 
    # if (T<4)  stop("'x' too short.")
    
    x.nextr = x
    ext_num = 999L  # neq 0 is OK.
    
    iter = 1L
    while ((iter <= 1000L) & (ext_num > 0L))  {  # number of iterator reach 1000 or ... 
      # sp = SpMAverage(x.nextr) # 迭代过程中对替换后的序列计算得到新的Spencer序列[这样很容易将谷抬得太高,把峰抹的太平]
      ratio = x.nextr/sp  #ratio of x and Spencer series.
      sd = sd(ratio)
      mn = mean(ratio)
      threshold_ceil = mn + limit * sd # unit value \ ±3.5个标准差
      threshold_floor = mn - limit * sd 
      
      isextreme = ((ratio > threshold_ceil) | (ratio < threshold_floor))   # boole values express whether it is extreme.
      
      if (iter==1L) x.isextreme = isextreme  # one of return values, identify who is extreme value.
      
      ext_dx = which(isextreme == TRUE)  # store TRUE 's index in original series.
      l_dx = length(ext_dx)
      # 
      # if (l_dx<=0L) break   # extreme's number is 0, thus break.
      
      k = 2L  # 表示前后各取两项, 算平均值.
      for (i in (1:l_dx))  {
        x.temp = x.nextr
        t = ext_dx[i] # 定位到原序列中异常值的索引
        left = max((t-k), 1)  # bingo
        right = min((t+k), T)
        x.temp[t] = (sum(x.nextr[left:right], na.rm = TRUE) - x.nextr[t]) / (length(x.nextr[left:right])-1)
        mink = right - left  # 取均值最少需要这么多个非极端值(1:k和T-k:T与中间的元素用到均值项数不一致).
        if (sum(!isextreme[left:right], na.rm = TRUE) >= mink)  isextreme[t] = FALSE  # 取平均的时候有mink个值为非异常值时, 即认为其变为非极端值, 即使数理上偏离Spencer曲线, 所以将对应位置的isextreme标记为FALSE.
        
      }
      ext_num = sum(isextreme, na.rm = TRUE)  # ext_num == 0 时即意味着没有极端值了.
      x.nextr = x.temp
      
      iter = iter + 1L
      
    }
    
    # 设定返回值 #
    limits = c(threshold_floor, threshold_ceil )  # 上下界限
    iteration = iter - 1L  # 迭代了多少次
    if (l_dx > 0 ) {
      warning("'x' maybe has some too extremely extreme values, which coudn't be adjusted by iterating 1000 times.")
      message("recommend use substitute method of 'spencer' when method of 'nearest' takes no effect, i.e. the values of the Spencer curve are substituted for the extreme values in the original series.")
      lst = list(x.nextr = x.nextr, extr_infr = data.frame(x=x, x.isextreme=x.isextreme, x.isextrextr=isextreme, x.nextr = x.nextr), limits = limits, iteration = iteration)
      
      # 列表list
      # x.nextr : 极端值修正后的x序列
      # x : 原序列
      # x.isextreme : 标识原序列中该行是否为极端值, TRUE:是极端值, FALSE:不是极端值.
      # x.isextrextr : 极端异常值, 迭代10000次仍无法修正的.
    }  
    else {
      lst = list(x.nextr = x.nextr, extr_infr = data.frame(x=x, x.isextreme=x.isextreme, x.nextr = x.nextr), limits = limits, iteration = iteration)
    }
    
    lst
    
  }
  
  # select subfunction by "subs_method" #
  if (class(subs_method) == "character") {
    method_n = switch(match.arg(subs_method, c("spencer", "nearest")),
                      spencer = 1L,
                      nearest = 2L
    )
    
  }
  
  else {
    if (length(subs_method) > 1) warning("length of 'subs_method' is greater than 1, the other will be omit if not 1.")
    method_n = subs_method[1]
  }          
  ans = switch (method_n,
                extreme_sp(x = x, limit = limit),
                extreme_nrst(x = x, limit = limit)
  )
  subs_num = sum(ans$extr_infr$x.isextreme, na.rm = TRUE)
  message(paste(subs_num, "extreme values have been substituted."))
  
  ans
  
}




