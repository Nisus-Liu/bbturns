MCD = function(x, maxgap = 8) {
  #' @title Capture \emph{MCD} value
  #' @description 
  #' return value of Month of Cyclical Dominance(MCD).
  #' @param maxgap max gap for calculate ratio when getting \emph{MCD} value, default is 8.
  #' @export MCD
  #' 
  #' @examples 
  #' x =seq(80, to = 120, length.out = 100)
  #' MCD(x = x, maxgap = 20)
  #' 
  #' 
  # 计算出一个mcd值,使得循环序列的波动刚好超过不规则序列的波动.
  # x  --数值向量
  # maxgap   --计算变化率的最大间隔
  # Note: 输入的x向量尽量保持在100附近,这样一可避免负值,二可避免小数位太多,导致精度引起的差异.
  # MCD : unite value
  
  
  if (any(is.na(x))) stop("'x' has 'NA'.")
  
  T = length(x)
  
  
  sp = SpMAverage(x)  # 调用SpMAverage()函数
  irr = x/sp *100 #分离出不规则要素序列
  
  for (j in 1:maxgap) {
    irr_r_h = rep(NA,j)
    sp_r_h = rep(NA,j)  # ratios in head.     
    
    # mean of abs of irrigation
    irr_i = irr[1:(T-j)]
    irr_ij = irr[(j+1):T]  # 间隔j的元素组成的序列
    irr_r = abs(irr_ij/irr_i - 1) * 100
    irr_r = c(irr_r_h, irr_r)
    irr_r_m = mean(irr_r, na.rm = T)
    # mean of abs of Spencer
    sp_i = sp[1:(T-j)]
    sp_ij = sp[(j+1):T]
    sp_r = abs(sp_ij/sp_i - 1) * 100
    sp_r = c(sp_r_h, sp_r)
    sp_r_m = mean(sp_r, na.rm = T)
    
    irr_sp =  irr_r_m / sp_r_m * 100 
    if (irr_sp<100 & j>=3) {
      MCD = j
      return(MCD)
    }
  }
  
}


