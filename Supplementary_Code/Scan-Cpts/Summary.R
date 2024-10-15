rm(list = ls())

powerBSD = rep(100, 20)
TrateBSD = matrix(0, 100, 20)
FrateBSD = matrix(0, 100, 20)
DifTiBSD = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/MES-Mulist', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD[j, i] = as.numeric(res[[j]]$DifTBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Review-Scan-Cpts/MES-BSD_Describe.RData')

#####################################################################################################

rm(list = ls())

powerSCQ = rep(100, 20)
TrateSCQ = matrix(0, 100, 20)
FrateSCQ = matrix(0, 100, 20)
DifTiSCQ = matrix(0, 100, 20)

for (i in 1:20)
{
  load(paste0('Review-Scan-Cpts/MES-Mulist', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    Resj = as.matrix(res[[j]]$reSCQ$SCAN_res)
    
    if (ncol(Resj) == 1) Resj = t(Resj)
    
    if (sum(Resj) == 1)
    {
      powerSCQ[i] = powerSCQ[i] - 1
    }
    else
    {
      DifTiSCQ[j, i] = as.numeric(res[[j]]$DifTSCQ, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj[k, 2]:Resj[k, 3])] = 1
      }
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSCQ[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSCQ[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

save(list = c('FrateSCQ', 'TrateSCQ', 'powerSCQ', 'DifTiSCQ'), file = 'Review-Scan-Cpts/MES-SCQ_Describe.RData')
###############################################################################################################

rm(list = ls())

powerSBS = rep(100, 20)
TrateSBS = matrix(0, 100, 20)
FrateSBS = matrix(0, 100, 20)
DifTiSBS = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/MES-Mulist', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (length(res[[j]]$reSBS$SDDCF_Res) != 0)
    {
      Resj = res[[j]]$reSBS$SDDCF_Res
      
      DifTiSBS[j, i] = as.numeric(res[[j]]$DifTSBS, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerSBS[i] = powerSBS[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSBS[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSBS[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateSBS', 'TrateSBS', 'powerSBS', 'DifTiSBS'), file = 'Review-Scan-Cpts/MES-SBS_Describe.RData')

#####################################################################################################

rm(list = ls())

powerWBS = rep(100, 20)
TrateWBS = matrix(0, 100, 20)
FrateWBS = matrix(0, 100, 20)
DifTiWBS = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/MES-Mulist', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (length(res[[j]]$reWBS$WDDCF_Res) != 0)
    {
      Resj = res[[j]]$reWBS$WDDCF_Res
      
      DifTiWBS[j, i] = as.numeric(res[[j]]$DifTWBS, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerWBS[i] = powerWBS[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateWBS[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateWBS[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateWBS', 'TrateWBS', 'powerWBS', 'DifTiWBS'), file = 'Review-Scan-Cpts/MES-WBS_Describe.RData')
###################################################################################################################
###########################################################################################################################################################################################################

rm(list = ls())

powerBSD = rep(100, 20)
TrateBSD = matrix(0, 100, 20)
FrateBSD = matrix(0, 100, 20)
DifTiBSD = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/MNS-Mulist', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD[j, i] = as.numeric(res[[j]]$DifTBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Review-Scan-Cpts/MNS-BSD_Describe.RData')

#####################################################################################################

rm(list = ls())

powerSCQ = rep(100, 20)
TrateSCQ = matrix(0, 100, 20)
FrateSCQ = matrix(0, 100, 20)
DifTiSCQ = matrix(0, 100, 20)

for (i in 1:20)
{
  load(paste0('Review-Scan-Cpts/MNS-Mulist', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    Resj = as.matrix(res[[j]]$reSCQ$SCAN_res)
    
    if (ncol(Resj) == 1) Resj = t(Resj)
    
    if (sum(Resj) == 1)
    {
      powerSCQ[i] = powerSCQ[i] - 1
    }
    else
    {
      DifTiSCQ[j, i] = as.numeric(res[[j]]$DifTSCQ, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj[k, 2]:Resj[k, 3])] = 1
      }
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSCQ[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSCQ[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

save(list = c('FrateSCQ', 'TrateSCQ', 'powerSCQ', 'DifTiSCQ'), file = 'Review-Scan-Cpts/MNS-SCQ_Describe.RData')
###############################################################################################################

rm(list = ls())

powerSBS = rep(100, 20)
TrateSBS = matrix(0, 100, 20)
FrateSBS = matrix(0, 100, 20)
DifTiSBS = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/MNS-Mulist', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (length(res[[j]]$reSBS$SDDCF_Res) != 0)
    {
      Resj = res[[j]]$reSBS$SDDCF_Res
      
      DifTiSBS[j, i] = as.numeric(res[[j]]$DifTSBS, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerSBS[i] = powerSBS[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSBS[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSBS[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateSBS', 'TrateSBS', 'powerSBS', 'DifTiSBS'), file = 'Review-Scan-Cpts/MNS-SBS_Describe.RData')

#####################################################################################################

rm(list = ls())

powerWBS = rep(100, 20)
TrateWBS = matrix(0, 100, 20)
FrateWBS = matrix(0, 100, 20)
DifTiWBS = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/MNS-Mulist', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (length(res[[j]]$reWBS$WDDCF_Res) != 0)
    {
      Resj = res[[j]]$reWBS$WDDCF_Res
      
      DifTiWBS[j, i] = as.numeric(res[[j]]$DifTWBS, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerWBS[i] = powerWBS[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateWBS[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateWBS[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateWBS', 'TrateWBS', 'powerWBS', 'DifTiWBS'), file = 'Review-Scan-Cpts/MNS-WBS_Describe.RData')
###################################################################################################################
############################################################################################################################################################################################################


rm(list = ls())

powerBSD = rep(100, 20)
TrateBSD = matrix(0, 100, 20)
FrateBSD = matrix(0, 100, 20)
DifTiBSD = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/WES-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD[j, i] = as.numeric(res[[j]]$DifTBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Review-Scan-Cpts/WES-BSD_Describe.RData')

#####################################################################################################

rm(list = ls())

powerSCQ = rep(100, 20)
TrateSCQ = matrix(0, 100, 20)
FrateSCQ = matrix(0, 100, 20)
DifTiSCQ = matrix(0, 100, 20)

for (i in 1:20)
{
  load(paste0('Review-Scan-Cpts/WES-Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    Resj = as.matrix(res[[j]]$reSCQ$SCAN_res)
    
    if (ncol(Resj) == 1) Resj = t(Resj)
    
    if (sum(Resj) == 1)
    {
      powerSCQ[i] = powerSCQ[i] - 1
    }
    else
    {
      DifTiSCQ[j, i] = as.numeric(res[[j]]$DifTSCQ, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj[k, 2]:Resj[k, 3])] = 1
      }
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSCQ[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSCQ[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

save(list = c('FrateSCQ', 'TrateSCQ', 'powerSCQ', 'DifTiSCQ'), file = 'Review-Scan-Cpts/WES-SCQ_Describe.RData')
###############################################################################################################

rm(list = ls())

powerSBS = rep(100, 20)
TrateSBS = matrix(0, 100, 20)
FrateSBS = matrix(0, 100, 20)
DifTiSBS = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/WES-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (length(res[[j]]$reSBS$SDDCF_Res) != 0)
    {
      Resj = res[[j]]$reSBS$SDDCF_Res
      
      DifTiSBS[j, i] = as.numeric(res[[j]]$DifTSBS, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerSBS[i] = powerSBS[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSBS[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSBS[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateSBS', 'TrateSBS', 'powerSBS', 'DifTiSBS'), file = 'Review-Scan-Cpts/WES-SBS_Describe.RData')

#####################################################################################################

rm(list = ls())

powerWBS = rep(100, 20)
TrateWBS = matrix(0, 100, 20)
FrateWBS = matrix(0, 100, 20)
DifTiWBS = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/WES-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (length(res[[j]]$reWBS$WDDCF_Res) != 0)
    {
      Resj = res[[j]]$reWBS$WDDCF_Res
      
      DifTiWBS[j, i] = as.numeric(res[[j]]$DifTWBS, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerWBS[i] = powerWBS[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateWBS[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateWBS[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateWBS', 'TrateWBS', 'powerWBS', 'DifTiWBS'), file = 'Review-Scan-Cpts/WES-WBS_Describe.RData')
###################################################################################################################
#########################################################################################################################################################################################################


rm(list = ls())

powerBSD = rep(100, 20)
TrateBSD = matrix(0, 100, 20)
FrateBSD = matrix(0, 100, 20)
DifTiBSD = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/WNS-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD[j, i] = as.numeric(res[[j]]$DifTBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Review-Scan-Cpts/WNS-BSD_Describe.RData')

#####################################################################################################

rm(list = ls())

powerSCQ = rep(100, 20)
TrateSCQ = matrix(0, 100, 20)
FrateSCQ = matrix(0, 100, 20)
DifTiSCQ = matrix(0, 100, 20)

for (i in 1:20)
{
  load(paste0('Review-Scan-Cpts/WNS-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    Resj = as.matrix(res[[j]]$reSCQ$SCAN_res)
    
    if (ncol(Resj) == 1) Resj = t(Resj)
    
    if (sum(Resj) == 1)
    {
      powerSCQ[i] = powerSCQ[i] - 1
    }
    else
    {
      DifTiSCQ[j, i] = as.numeric(res[[j]]$DifTSCQ, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj[k, 2]:Resj[k, 3])] = 1
      }
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSCQ[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSCQ[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

save(list = c('FrateSCQ', 'TrateSCQ', 'powerSCQ', 'DifTiSCQ'), file = 'Review-Scan-Cpts/WNS-SCQ_Describe.RData')
###############################################################################################################

rm(list = ls())

powerSBS = rep(100, 20)
TrateSBS = matrix(0, 100, 20)
FrateSBS = matrix(0, 100, 20)
DifTiSBS = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/WNS-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (length(res[[j]]$reSBS$SDDCF_Res) != 0)
    {
      Resj = res[[j]]$reSBS$SDDCF_Res
      
      DifTiSBS[j, i] = as.numeric(res[[j]]$DifTSBS, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerSBS[i] = powerSBS[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSBS[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSBS[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateSBS', 'TrateSBS', 'powerSBS', 'DifTiSBS'), file = 'Review-Scan-Cpts/WNS-SBS_Describe.RData')

#####################################################################################################

rm(list = ls())

powerWBS = rep(100, 20)
TrateWBS = matrix(0, 100, 20)
FrateWBS = matrix(0, 100, 20)
DifTiWBS = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-Scan-Cpts/WNS-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (length(res[[j]]$reWBS$WDDCF_Res) != 0)
    {
      Resj = res[[j]]$reWBS$WDDCF_Res
      
      DifTiWBS[j, i] = as.numeric(res[[j]]$DifTWBS, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerWBS[i] = powerWBS[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateWBS[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateWBS[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateWBS', 'TrateWBS', 'powerWBS', 'DifTiWBS'), file = 'Review-Scan-Cpts/WNS-WBS_Describe.RData')
###################################################################################################################






















 