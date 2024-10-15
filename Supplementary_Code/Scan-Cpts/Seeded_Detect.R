SeedDetect = function(X, Y, m, MB = 1000, alpha = 0.05, decay = sqrt(2))
{
  p = ncol(X)
  intervals = seeded.intervals(p, m)
  
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
  
  SDDCF_Res = data.frame(startind, endind)
  
  return(list(SDDCF_Res = SDDCF_Res))
}


seeded.intervals = function(n, m = 10, decay = 2, unique.int = F)
{
  n	= as.integer(n)
  depth	= log(n, base = decay)
  depth	= ceiling(depth)
  M	= sum(2^(1:depth)-1)
  
  boundary_mtx = matrix(NA, ncol = 2)
  colnames(boundary_mtx) = c("st", "end")
  boundary_mtx[1, ] = c(1, n)
  
  depth	<- log(n, base = decay)
  depth	<- ceiling(depth)
  
  
  for(i in 2:depth){
    int_length	<- n * (1/decay)^(i-1)
    
    n_int		<- ceiling(round(n/int_length, 14))*2-1		# sometimes very slight numerical inaccuracies
    
    boundary_mtx	<- rbind(boundary_mtx,
                          cbind(floor(seq(1, n-int_length, length.out = (n_int))), 
                                ceiling(seq(int_length, n, length.out = (n_int)))))
  }
  
  if(unique.int) boundary_mtx = unique(boundary_mtx);
  
  
  length.interval = boundary_mtx[, 2] - boundary_mtx[, 1]
  index_include = which(length.interval >= m)
  
  boundary_mtx = boundary_mtx[index_include, ]
  return(boundary_mtx)
}
