#' turns of a time series 
#' @description 
#' \code{bbturns} returns turns information of a time series, selected by B-B algorithm.
#' 
#' @details 
#' B-B algorithm proposed by Gerhard Bry and Charlotte Boschan in 1971, see \url{http://papers.nber.org/books/bry_71-1}.
#' 
#' @param x a numeric vector.
#' @param neighbour find local extreme within \eqn{+/- neighbour}. 
#' @param mincycle numeric, one cycle must have at least \code{mincycle} months, if a monthly time series.
#' @param phase numeric, a phase of one cycle must have length of \code{phase} months, phase means peak to peak or trough to trough. 
#' @param extr_limit numeric, some values will be regarded as extreme who out range of \eqn{+/- extr_limit} standard deviation of ratio of irregular and original series.
#' @param end_num peak or trough should be apart form either end of series by at least \code{end_num} months.
#' @param mcd_gap max gap for calculate ratio when getting \emph{MCD} value, default is 8.
#' @aliases bbturns
#' @return 
#' \itemize{
#'  \item \code{isturn} a vector stores turns information of series. 1 means it is a peak turn, -1 means it is trough turn, 0 means neither peak nor trough.
#'  \item \code{isturn_old} a vector stores turns information for last setup.
#'  \item \code{peaktrough_idx} stores index of peak or trough turns of original series.
#'  \item \code{peak_idx} peak index.
#'  \item \code{trough_idx} trough index.
#' }
#' 
#' @examples 
#' data(importexport)
#' bbturns(importexport$ImEx)
#' @export bbturns
#' @author Hejun Liu, Beijing \email{liuhejunlj@@163.com} 
#' 
#' 


bbturns = function(x, neighbour = 5, mincycle = 15, phase = 5, extr_limit = 3.5, end_num = 6, mcd_gap = 8) {
  
# 
#   x : 向量
#   neighbour : 前后多少期内寻找极值
#   end_num : 两端转折点距离端点的最低长度
#   mincycle : 最短周期长度('峰-峰'or'谷-谷')
#   mcd_gap : 计算MCD值得时候计算变换率最大间隔期,默认为8.
#   phase : 分段最低长度
#   extr_limit : 排除异常值所用的阈值,默认超过3.5倍标准差为异常值,以Spencer值替代.

# 子函数 #
# ===========================================================

## 局部极值点 ##

### 在12期移动平均中寻找局部极值 ###

localextr_12 = function(x, neighbour=neighbour, end_num=end_num)   {
  
  if (any(is.na(x))) stop("'x' has 'NA'.")
  
  # err = mean(x) / 1000L  # 在比较大小的时候考虑到精度的影响,允许的误差
  
  srs_len = length(x) 
  isturn = rep(0, srs_len)  
  for (t in end_num:(srs_len-end_num+1)) {
    
    lft= t-neighbour
    rgt = t+neighbour
    mx = max(x[lft:rgt])
    mi = min(x[lft:rgt])
    # 局部极大值:[isturn]标记为1, 极小值标记为-1, 都不是为0.
    if (x[t] >= mx)  isturn[t] = 1
    if (x[t] <= mi)  isturn[t] = -1
  }
  
  peak_idx = which(isturn == 1)  # 峰对应在原序列中的索引
  trough_idx = which(isturn == -1) # 谷对应在原序列中的索引
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  lst = list(isturn = isturn, peaktrough_idx=peaktrough_idx, peak_idx=peak_idx,trough_idx=trough_idx)
  lst
  
}

### 后期在Spencer/MCD/原...中寻找局部极值点的函数 ###

localextr_med = function(x, neighbour=neighbour, isturn, peak_idx, trough_idx)   {
 
  # determine tentative turns by searching local maximum or minimum .
  # 
  # x : 'Spencer curve'/'MCD curve'/'original series'
  # neighbour :  calculate max /min of ±neighbour observations.
  # isturn : a vector stores information of tentative turns determined from '12-month moving average'/'Spencer curve'/'MCD curve'/'original series', peak : 1; trough :-1;otherwise :0. This is previous returns, similarly hereinafter.
  # peak_idx : stores index of peak in 12-month moving average series ...
  # trough_idx : stores index of trough in ...
  
  if (any(is.na(x))) stop("'x' has 'NA'.")
  if (any(is.na(isturn))) stop("'isturn' has 'NA'.")
  if (any(is.na(peak_idx))) stop("'peak_idx' has 'NA'.")
  if (any(is.na(trough_idx))) stop("'trough_idx' has 'NA'.")
  
  isturn_old = isturn
  
  srs_len = length(x) 
  # isturn = rep(0, srs_len)  
  pk_len = length(peak_idx)
  tr_len = length(trough_idx)
  # adjust peak #
  for (i in 1:pk_len) {
    t = peak_idx[i] # locate the 't' of original series 
    lft =max(( t-neighbour), 1)
    rgt = min((t+neighbour), srs_len)
    mx = max(x[lft:rgt])
    
    for (j in (lft:rgt)) {
      if (x[j] >= mx) {
        
        isturn[t] = 0
        isturn[j] = 1  # !!!和上一行的顺序千万不能错, 当'j'遍历到刚好等于't'时,如果顺序反了,即前脚将局部最大值标记为1,后脚又将其标记为了0.顺序对了则能保证位置没变(本曲线和上次曲线的转折点位置)的暂定转折点依然标记为转折点.
        break
      }
      
    }
    
    
  }
  
  # adjust trough 
  for (i in 1:tr_len) {
    t = trough_idx[i] # locate the 't' of original series 
    lft =max(( t-neighbour), 1)
    rgt = min((t+neighbour), srs_len)
    mi = min(x[lft:rgt])
    
    for (j in (lft:rgt)) {
      if (x[j] <= mi) {
        
        isturn[t] = 0
        isturn[j] = -1   
        break
      }
      
    }
  }
  
  
  peak_idx = which(isturn == 1)  # 峰对应在原序列中的索引
  trough_idx = which(isturn == -1) # 谷对应在原序列中的索引
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  lst = list(isturn = isturn,isturn_old=isturn_old, peaktrough_idx=peaktrough_idx, peak_idx=peak_idx,trough_idx=trough_idx)
  lst
  
}

## 峰谷交错 ##

peaktroughcross = function(x, isturn , peaktrough_idx) {
 
  # x : series of original or moving average ...
  # isturn : tentative turns
  # peaktrough_idx : index of peak or trough .

  pt_len = length(peaktrough_idx)
  pt_idx = peaktrough_idx
  isturn_old = isturn
  
  iter = 1L
  
  while (TRUE)   {
    i = 1L
    while (i < pt_len) {
      
      t = pt_idx[i]
      tnext = pt_idx[i+1] # 下一个转折点的索引
      # 峰 #
      if (isturn[t] == 1L) {
        if (isturn[tnext] == 1L) {
          if (x[t] > x[tnext])  {
            isturn[tnext] = 0L
            i = i+1L   # 下一个被排除了, 所以不必再和后面的对比了, 和后面的位移一起实现位移2个. 后面的后面依旧是峰怎么办? 迭代解决!
          }
          else {
            isturn[t] = 0L
          }
        }
      }
      # 谷 #
      if (isturn[t] == -1L) {
        if (isturn[tnext] == -1L) {
          if (x[t] < x[tnext])  {
            isturn[tnext] = 0L
            i = i+1L   # 下一个被排除了, 所以不必再和后面的对比了, 和后面的位移一起实现位移2个. 后面的后面依旧是峰怎么办? 迭代解决!
          }
          else {
            isturn[t] = 0L
          }
        }
      }
      
      i = i + 1L
      
    }
    
    pt_num = isturn[isturn != 0] # 取出标记峰谷的数字(1,-1)
    ref = if (pt_num[1] == 1L)  rep(c(1L,-1L), length = length(pt_num)) else rep(c(-1L,1L), length = length(pt_num))
    if (all(pt_num == ref)) break  # 满足, 则峰谷交错.
    
    iter =+ 1L  # return number of iteration .
  }  
  
  # return values setting 
  
  
  peak_idx = which(isturn == 1)  # 峰对应在原序列中的索引
  trough_idx = which(isturn == -1) # 谷对应在原序列中的索引
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  lst = list(isturn = isturn, isturn_old = isturn_old , peaktrough_idx=peaktrough_idx, peak_idx=peak_idx,trough_idx=trough_idx , iter_num = iter)
  lst
  
}

# test #
# peaktroughcross(x=imex_12, isturn = localextr12$isturn, peaktrough_idx = localextr12$peaktrough_idx)

## 最小周期准则 ##

mincycle = function(x, mincycle = mincycle, isturn, peak_idx , trough_idx)  {
  # 最小周期测试, 默认"峰-峰"或者"谷-谷"的至少跨15个月.
  # 
  # x : 向量
  # isturn : 转折点信息
  # mincycle : 指定最短周期长度
  # peak_idx : '峰'的索引值
  # trough_idx : '谷'的索引值

  
  isturn_old = isturn
  
  pk_len = length(peak_idx)
  tr_len = length(trough_idx)
  # peak-peak #
  i = 0L
  while (i >= pk_len)  {
    i =+ 1L
    t = peak_idx[i]
    tnext = peak_idx[i+1]
    gap = tnext-t+1
    if (gap < mincycle)  {
      if (x[t] > x[tnext]) {
        isturn[tnext] = 0L
        i =+ 1L   # t+1排除后,计步器要多移动一位
      }
      else {
        isturn[t] = 0L
      }
    }
  }
  
  # trough-trough #
  i = 0L
  while (i >= tr_len) {
    i =+ 1L
    t = trough_idx[i]
    tnext = trough_idx[i+1]
    gap = tnext-t+1
    if (gap < mincycle) {
      if (x[t] < x[tnext]) {
        isturn[tnext] = 0L  # !!!如果下一个被排除,那么下一个就不能用来后下下一个谷对比了,因为:当出现t,t+1,t+2连续两个周期不足15月时,如果t+2比t+1更小,那不会产生什么问题,但是如果t+2比t+1大,那么程序会剔除t+2,最后相当于三个谷点钟仅保留了t,t+1和t+2被剔除了,这便有过度约束的嫌疑,因为这种情况保留t和t+2其实就可以了!
        i =+ 1L # 计步器多移一位,跳过t+1,直接比较t+2
        
      }
      else {
        isturn[t] = 0L
      }
      
    }
  }
  
  # 峰谷交错调整 #
  peaktrough_idx = which(isturn != 0L) 
  
  ans = peaktroughcross(x=x, isturn = isturn, peaktrough_idx = peaktrough_idx)
  lst = list(isturn = ans$isturn, isturn_old = isturn_old , peaktrough_idx=ans$peaktrough_idx, peak_idx=ans$peak_idx,trough_idx=ans$trough_idx )
  lst
}




# test #
# debug(mincycle)
# mincycle(x=imex_sp, mincycle = 15, isturn=localextr_sp$isturn, peak_idx=localextr_sp$peak_idx , trough_idx=localextr_sp$trough_idx)


## 两端6月准则 ##

ends_min_gap = function(x, isturn, min_gap = end_num)  {

  # 两端不足6月的暂定转折点要剔除
  # 
  # x : 向量
  # isturn : 转折点信息
  # min_gap : 指定距离端点的最短距离
  # 

  
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  isturn_old = isturn
  x_len = length(x)
  t1 = peaktrough_idx[1]
  tn = peaktrough_idx[length(peaktrough_idx)]
  if (t1 < min_gap) {
    isturn[t] = 0L
  }
  if (tn > (x_len-min_gap+1)) {
    isturn[tn] = 0L
  }
  
  
  peak_idx = which(isturn == 1)  # 峰对应在原序列中的索引
  trough_idx = which(isturn == -1) # 谷对应在原序列中的索引
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  lst = list(isturn = isturn, isturn_old = isturn_old , peaktrough_idx=peaktrough_idx, peak_idx=peak_idx,trough_idx=trough_idx )
  lst
  
}




# test #
# ends_min_gap(x=imex_sp,isturn=mincycle_sp$isturn, min_gap = 6)

## 端点附近峰谷值测试 ##

ends_value = function(x, isturn)  {
  # 两端峰谷值测试
  # 
  # x : 向量
  # isturn : 转折点信息

  
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  isturn_old = isturn
  x_len = length(x)
  t1 = peaktrough_idx[1]
  tn = peaktrough_idx[length(peaktrough_idx)]
  # head turn value test #
  if (isturn[t1]==1L) {
    mx = max(x[1:t1])
    if (x[t1]<mx) {isturn[t1] = 0L}
  }
  if (isturn[t1]==-1L) {
    mi = min(x[1:t1])
    if (x[t1]>mi) {isturn[t1] = 0L}
  }
  
  # tail turn value test #
  if (isturn[tn]==1L) {
    mx = max(x[tn:x_len])
    if (x[tn]<mx) {isturn[tn] = 0L}
  }
  if (isturn[tn]==-1L) {
    mi = min(x[tn:x_len])
    if (x[tn]>mi) {isturn[tn] = 0L}
  }
  
  
  peak_idx = which(isturn == 1)  # 峰对应在原序列中的索引
  trough_idx = which(isturn == -1) # 谷对应在原序列中的索引
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  lst = list(isturn = isturn, isturn_old = isturn_old , peaktrough_idx=peaktrough_idx, peak_idx=peak_idx,trough_idx=trough_idx )
  lst
  
}




# test #
# ends_value(x=imex_nextr$x.nextr,isturn=endsgap_org$isturn)

## 最小分阶段 ##

minphase = function(x, isturn, phase = phase)  {
  # 两端峰谷值测试
  # 
  # x : 向量
  # isturn : 转折点信息
  # phase : 分阶段最低长度,不足该长度的分段则排除
  # 

  
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  isturn_old = isturn
  x_len = length(x)
  pt_len =length(peaktrough_idx)
  
  for (i in (1:(pt_len-1L))) {
    
    t = peaktrough_idx[i]
    tnext = peaktrough_idx[i+1]
    if ((tnext - t + 1) < phase)  {
      isturn[t] = 0L
      isturn[tnext] = 0L
    }
    
  }
  
  
  
  peak_idx = which(isturn == 1)  # 峰对应在原序列中的索引
  trough_idx = which(isturn == -1) # 谷对应在原序列中的索引
  peaktrough_idx = which(isturn != 0) # 峰和谷的索引
  
  lst = list(isturn = isturn, isturn_old = isturn_old , peaktrough_idx=peaktrough_idx, peak_idx=peak_idx,trough_idx=trough_idx )
  lst
  
}




# test #

# minphase(x=imex_nextr$x.nextr,isturn=mincycle_org$isturn, phase=5)









# 12期移动平均曲线中选择暂定转折点 #
# ===================================================================

imex_nextr = extrmadj(x, limit = extr_limit)   # 异常值修正
x_nextr = imex_nextr$x.nextr # 异常值修正后的原序列

x_12ma = cmAverage(x_nextr, 12)
    # x:经过异常值修正后的修正序列,计算得到12期移动平均
localextr12 = localextr_12(x_12ma, neighbour = neighbour, end_num = end_num)

tt_12 = peaktroughcross(x=x_12ma, isturn = localextr12$isturn, peaktrough_idx = localextr12$peaktrough_idx)   # 12MA曲线最终转折点相关结果

# Spencer曲线中...#
# =================================================================

x_sp = SpMAverage(x_nextr)

localextr_sp = localextr_med(x_sp, neighbour = neighbour, isturn = tt_12$isturn, peak_idx = tt_12$peak_idx, trough_idx = tt_12$trough_idx)

mincycle_sp = mincycle(x=x_sp, mincycle = mincycle, isturn=localextr_sp$isturn, peak_idx=localextr_sp$peak_idx , trough_idx=localextr_sp$trough_idx)

tt_sp = ends_min_gap(x=x_sp,isturn=mincycle_sp$isturn, min_gap = end_num) # Spencer曲线最终转折点相关结果


# MCD曲线中... #
# ==================================================================

mcd = MCD(x_nextr, maxgap = mcd_gap)
x_mcd = cmAverage(x_nextr, mcd)

localextr_mcd = localextr_med(x = x_mcd, neighbour = neighbour, isturn = tt_sp$isturn, peak_idx = tt_sp$peak_idx, trough_idx = tt_sp$trough_idx)

mincycle_mcd = mincycle(x=x_mcd, mincycle = mincycle, isturn=localextr_mcd$isturn, peak_idx=localextr_mcd$peak_idx , trough_idx=localextr_mcd$trough_idx)

tt_mcd = ends_min_gap(x_mcd, isturn = localextr_mcd$isturn)

# 原序列(非平滑序列)中... #
# =======================================================================

localextr_org = localextr_med(x = x_nextr, neighbour = max(4, mcd),  isturn = tt_mcd$isturn, peak_idx = tt_mcd$peak_idx, trough_idx = tt_mcd$trough_idx)

endsgap_org = ends_min_gap(x = x_nextr, isturn = localextr_org$isturn, min_gap = end_num)

endsval_org = ends_value(x=x_nextr, isturn=endsgap_org$isturn)

mincycle_org = mincycle(x=x_nextr, mincycle = mincycle, isturn=endsval_org$isturn, peak_idx=endsval_org$peak_idx , trough_idx=endsval_org$trough_idx)

tt_final = minphase(x=x_nextr,isturn=mincycle_org$isturn, phase=phase)


# 返回值 #

tt_final




}




# test #  
# bbturns(bb$ImEx)


