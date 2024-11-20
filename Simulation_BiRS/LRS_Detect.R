ScanM = function(X, Y, Lmin, Lmax, skip, MB = 1000, alpha = 0.05)
{
  p = ncol(X); n = nrow(X); m = nrow(Y)
  intervals = scan.intervals(p, Lmin, Lmax, skip)
  
  Un = sqrt(n)*(apply(X, 2, mean) - apply(Y, 2, mean))
  EleVar = apply(X, 2, var) + apply(Y, 2, var)*(n/m)
  
  Un = Un/sqrt(EleVar)
  
  thresh = sqrt(2*log(p*Lmax))
  
  SelectStat = rep(0, nrow(intervals))
  for (i in 1:nrow(intervals))
  {
    Uni = Un[intervals[i, 1]:intervals[i, 2]]
    SelectStat[i] = sum(Uni)/sqrt(length(Uni))
  }
  
  index_stop = 1
  startind = c(); endind = c()
  while(index_stop > 0)
  {
    #print(paste0('Greedy Segmentation Step:', greedy_num))
    
    maxstat = max(SelectStat)
    
    if (maxstat < thresh) break;
    
    index_max = which(SelectStat == maxstat)
    #stat = SelectStat[index_max]
    
    st = intervals[index_max, 1]; et = intervals[index_max, 2]
    
    Belong_Func = function(x) return((x[1] > intervals[index_max, 2]) + (x[2] < intervals[index_max, 1]))
    Belong_index = apply(intervals, 1, Belong_Func)
    
    index_out = which(Belong_index == 1)
    intervals = as.matrix(intervals[index_out, ])
    
    SelectStat = SelectStat[index_out]
    #beta_diff = beta_diff[-index_belong]
    
    startind = c(startind, st); endind = c(endind, et)
    index_stop = length(intervals)
  }
  
  SCANM_Res = data.frame(startind, endind)
  
  return(list(SCANM_Res = SCANM_Res))
}

scan.intervals = function(n, Lmin, Lmax, skip)
{
  candidate_length = unique(c(seq(Lmin, Lmax, skip), Lmax))
  intervals = c()
  
  for (i in 1:length(candidate_length))
  {
    lengthi = candidate_length[i]
    starti = 1:(n - lengthi); endi = starti + lengthi
    intervali = cbind(starti, endi)
    
    intervals = rbind(intervals, intervali)
  }
  
  return(intervals)
}
