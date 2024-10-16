load('~/Biobank/BiRS_Results_disease_C50 Malignant neoplasm of breast_Sex_Female.RData')

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
  
  indbar = unique(indbar)
  totalsize = length(indbar)
  
  if (totalsize == 0) return(NULL)
  
  inddiff = indbar[-1] - indbar[-totalsize]
  indind = c(0, which(inddiff > 1), totalsize)
  
  NumCNV = length(indind) - 1
  startpoint = endpoint = CNVlength = Mean = rep(0,NumCNV)
  for (i in 1:NumCNV){
    startpoint[i] = indbar[indind[i] + 1]
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

load('FinalVariant.RData')
ThresV = apply(SneX, 1, max)
thresh = quantile(ThresV, 1 - 0.05)
S4_Thres = SSSS(x = UnVec, thresh = thresh, distance = 10, minlength = 1)

p = length(UnVec)

SIGN = UnVec
startind = S4_Thres$startpoint; endind = S4_Thres$endpoint
#SIGN = sign(UnVec1)

# P = rep(0, length(UnVec))
# for (i in 1:p)
# {
#   P[i] = pnorm(UnVec[i], 0, sd(SneX[, i]), lower.tail = F)
# }

rsidsel = NULL; possel = NULL; chrsel = NULL; cluster = NULL; Tvalue = NULL; direction = NULL
for (i in 1:length(startind))
{
  ind1 = startind[i]; ind2 = endind[i]
  len = ind2 - ind1 + 1
  cluster = c(cluster, rep(i, len))
  rsidsel = c(rsidsel, SNP[ind1:ind2])
  possel = c(possel, BP[ind1:ind2])
  chrsel = c(chrsel, CHR[ind1:ind2])
}

BreastCancer = data.frame(rsidsel, possel, chrsel, cluster)
write.csv(BreastCancer, file = 'C50Breast_4S.csv')




 
