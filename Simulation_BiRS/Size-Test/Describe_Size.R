rm(list = ls())
load('Setting_MES.RData')
sizeBSD = 0; sizeBSC = 0; sizeSCQ = 0

nsimu = length(res)
for (i in 1:nsimu)
{
  if(!is.character(res[[i]]$reBSD)) sizeBSD = sizeBSD + 1
  
  if(!is.character(res[[i]]$reBSC)) sizeBSC = sizeBSC + 1
  
  Resi = as.matrix(res[[i]]$reSCQ$SCAN_res)
  if (sum(Resi) != 1)
  {
    sizeSCQ = sizeSCQ + 1
  }
    
}


rm(list = ls())
load('Setting_MNS.RData')
sizeBSD = 0; sizeBSC = 0; sizeSCQ = 0

nsimu = length(res)
for (i in 1:nsimu)
{
  if(!is.character(res[[i]]$reBSD)) sizeBSD = sizeBSD + 1
  
  if(!is.character(res[[i]]$reBSC)) sizeBSC = sizeBSC + 1
  
  Resi = as.matrix(res[[i]]$reSCQ$SCAN_res)
  if (sum(Resi) != 1)
  {
    sizeSCQ = sizeSCQ + 1
  }
  
}


rm(list = ls())
load('Setting_WES.RData')
sizeBSD = 0; sizeBSC = 0; sizeSCQ = 0

nsimu = length(res)
for (i in 1:nsimu)
{
  if(!is.character(res[[i]]$reBSD)) sizeBSD = sizeBSD + 1
  
  if(!is.character(res[[i]]$reBSC)) sizeBSC = sizeBSC + 1
  
  Resi = as.matrix(res[[i]]$reSCQ$SCAN_res)
  if (sum(Resi) != 1)
  {
    sizeSCQ = sizeSCQ + 1
  }
  
}

rm(list = ls())
load('Setting_WNS.RData')
sizeBSD = 0; sizeBSC = 0; sizeSCQ = 0

nsimu = length(res)
for (i in 1:nsimu)
{
  if(!is.character(res[[i]]$reBSD)) sizeBSD = sizeBSD + 1
  
  if(!is.character(res[[i]]$reBSC)) sizeBSC = sizeBSC + 1
  
  Resi = as.matrix(res[[i]]$reSCQ$SCAN_res)
  if (sum(Resi) != 1)
  {
    sizeSCQ = sizeSCQ + 1
  }
  
}

rm(list = ls())
load('Setting_GEN.RData')
sizeBSD = 0; sizeBSC = 0; sizeSCQ = 0

nsimu = length(res)
for (i in 1:nsimu)
{
  if(!is.character(res[[i]]$reBSD)) sizeBSD = sizeBSD + 1
  
  if(!is.character(res[[i]]$reBSC)) sizeBSC = sizeBSC + 1
  
  Resi = as.matrix(res[[i]]$reSCQ$SCAN_res)
  if (sum(Resi) != 1)
  {
    sizeSCQ = sizeSCQ + 1
  }
  
}

