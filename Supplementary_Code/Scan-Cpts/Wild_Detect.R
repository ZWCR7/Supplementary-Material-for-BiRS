WildDetect = function(X, Y, M, m, MB = 1000, alpha = 0.05)
{
  p = ncol(X)
  intervals = random.intervals(p, M, m)
  
  DCFStat = XFTest(X, Y, MB, alpha)
  Un = DCFStat$UnVec; ThresMat = DCFStat$SneX
  
  EleVar = apply(X, 2, var) + apply(Y, 2, var)
  #ThresV = (apply(ThresMat^2, 1, sum) - sum(EleVar))/sqrt(p) 
  #thresh = quantile(ThresV, 1 - alpha)
  
  ThresInt = matrix(0, MB, nrow(intervals))
  SelectStat = rep(0, nrow(intervals))
  for (i in 1:nrow(intervals))
  {
    Uni = Un[intervals[i, 1]:intervals[i, 2]]
    Thresi = ThresMat[, intervals[i, 1]:intervals[i, 2]]
    
    Vari = EleVar[intervals[i, 1]:intervals[i, 2]]
    
    SelectStat[i] = (sum(Uni^2) - sum(Vari))/sqrt(length(Uni))
    ThresInt[, i] = (apply(Thresi^2, 1, sum) - sum(Vari))/sqrt(length(Uni)) 
  }
  
  ThresV = apply(ThresInt, 1, max)
  thresh = quantile(ThresV, 1 - alpha)
  
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
  
  WDDCF_Res = data.frame(startind, endind)
  
  return(list(WDDCF_Res = WDDCF_Res))
}

random.intervals = function (n, M, m, unique.int = T)
{
  n = as.integer(n)
  M = as.integer(M)
  intervals = matrix(0, nrow = M, ncol = 2)
  intervals[, 1] = floor(runif(M, min = 1, max = n))		# intended starting point
  intervals[, 2] = ceiling(runif(M, min = 1, max = n))	# intended end point
  
  for (i in 1:M)
  {
    tmp = intervals[i, ] 
    while(tmp[1] == tmp[2])
    {	                                # in case the two happen to be equal, take new ones
      tmp[1] = floor(runif(1, min = 1, max = n)) 
      tmp[2] = ceiling(runif(1, min = 1, max = n))
      intervals[i, ] = tmp
    }
    
    if(tmp[1] > tmp[2]) intervals[i, ] <- tmp[2:1]    	# make sure they are ordered
  }
  
  if(unique.int) intervals = unique(intervals)
  
  length.interval = intervals[, 2] - intervals[, 1]
  index_include = which(length.interval >= m)
  
  intervals = intervals[index_include, ]
  return(intervals)
}
