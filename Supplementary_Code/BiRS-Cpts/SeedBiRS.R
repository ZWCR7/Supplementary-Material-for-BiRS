source('SummaryBiRS.R')

SeedDetectBiRS = function(X, Y, decay = 2, m, trunc, MB = 1000, alpha = 0.05)
{
  p = ncol(X)
  intervals = seeded.intervals(p, m, decay)
  
  totalre = XFTest(X, Y, MB, alpha)
  UnVec = totalre$UnVec
  SneX = totalre$SneX
  
  Signal = rep(0, p)
  #if (!totalre$reject) return(Signal) 
  
  rm(totalre)
  
  Region_Frame = NULL
  
  Signal_temp = rep(0, p)
  
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

seeded.intervals = function(n, m = 10, decay = 2, unique.int = T)
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
