source('SummaryBiRS.R')

WildDetectBiRS = function(X, Y, M, m, trunc, MB = 1000, alpha = 0.05)
{
  p = ncol(X)
  intervals = random.intervals(p, M, m)
  
  totalre = XFTest(X, Y, MB, alpha)
  UnVec = totalre$UnVec
  SneX = totalre$SneX
  
  rm(totalre)
  
  Region_Frame = NULL
  
  Signal_temp = rep(0, p);  Signal = rep(0, p)
  for (i in 1:nrow(intervals))
  {
    foldleni = intervals[i, 2] - intervals[i, 1] + 1
    gapi = intervals[i, 1]:intervals[i, 2]
    
    Uni = UnVec[gapi]; Sni = SneX[, gapi]
    detecti = BiRS_Summary(Uni, Sni, foldleni, trunc, alpha)
    
    Unmaxi = -9999; starti = NULL; endi = NULL
    if (!is.character(detecti))
    {
      startindi = detecti$BSDCF_res$startind + intervals[i, 1] - 1
      endindi = detecti$BSDCF_res$endind + intervals[i, 1] - 1
      
      for (j in 1:length(startindi))
      {
        startindij = startindi[j]; endindij = endindi[j]
        Unmaxij = max(abs(UnVec[startindij:endindij]))
        
        Signal_temp[startindij:endindij] = 1
        Region_Frame = rbind(Region_Frame, c(startindij, endindij, Unmaxij)) 
      }
    }
  }
  
  if (!is.null(Region_Frame))
  {
    SneX_temp = SneX[, which(Signal_temp == 1)]; thres_final = quantile(apply(SneX_temp, 1, max), 1 - alpha)
    
    if (is.null(nrow(Region_Frame)))
    {
      Region_Frame =  matrix(Region_Frame, 1, 3)
    }
    
    for (i in 1:nrow(Region_Frame))
    {
      if (Region_Frame[i, 3] > thres_final)
      {
        Signal[Region_Frame[i, 1]:Region_Frame[i, 2]] = 1
      }
    }
  }
  
  return(Signal)
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
