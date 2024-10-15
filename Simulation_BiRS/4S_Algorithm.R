#super scalable short segment (4S) detection algorithm

SSSSDetect = function(X, Y, length.gap, Lmin, MB = 1000, alpha = 0.05)
{
  p = ncol(X)
  #intervals = seeded.intervals(p, m)
  
  DCFStat = XFTest(X, Y, MB, alpha)
  Un = DCFStat$UnVec; ThresMat = DCFStat$SneX
  
  ThresV = apply(ThresMat, 1, max)
  thresh = quantile(ThresV, 1 - alpha)
  
  S4_Thres = SSSS(x = Un, thresh = thresh, distance = length.gap, minlength = Lmin)
  
  startind = S4_Thres$startpoint; endind = S4_Thres$endpoint
  
  SSSS_Res = data.frame(startind, endind)
  return(list(SSSS_Res = SSSS_Res))
}

SSSS = function(x, center = FALSE, thresh = NULL, distance = 9, minlength = 3, inference = FALSE, pvthresh=0.05){ 
  if (center) x=x-median(x)                     # if center is TRUE, we center x by substracting its median
  #if (is.null(thresh)) thresh = 2*sd(x)         # if thresh is not given, we use 2*sd. We suggest that user should input a threshold
  if (is.null(thresh)) thresh = quantile(abs(x),0.95)
  n = length(x)
  ind = which(abs(x)>thresh)                    # Step 1: set of locations with measurements beyond the thresh 
  
  k = length(ind); indbar = c()
  
  if (k <= 1) return(NULL)
  
  for (i in 1:(k - 1))
  {
    if ((ind[i+1] - ind[i]) <= distance)
    {
      indbar = c(indbar, ind[i]:ind[i+1])
    }
  }
  
  totalsize = length(indbar)
  
  if (totalsize == 0) return(NULL)
  
  inddiff = indbar[-1] - indbar[-totalsize]
  indind = c(1, which(inddiff > 1), totalsize)
  
  NumCNV = length(indind) - 1
  startpoint = endpoint = CNVlength = Mean = rep(0,NumCNV)
  for (i in 1:NumCNV){
    startpoint[i] = indbar[indind[i]]
    endpoint[i] = indbar[indind[i+1]]
    CNVlength[i] = endpoint[i]+1-startpoint[i]
    Mean[i] = mean(x[startpoint[i]:endpoint[i]])
    #bk[i]   = length(subint[[i]])
  }
  object = data.frame(startpoint=startpoint, endpoint=endpoint,length=CNVlength,mean=Mean)
  # step 3 here
  if (minlength>1) {
    object=object[object$length>minlength,]
    rownames(object) = seq(nrow(object))
  }
  # addition step to calculate p-value
  return(object)
}

 