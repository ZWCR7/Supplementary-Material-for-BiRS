rm(list = ls())

powerBSD = rep(500, 20)
TrateBSD = matrix(0, 500, 20)
FrateBSD = matrix(0, 500, 20)
DifTiBSD = matrix(0, 500, 20)

for (i in 1:20)
{ 
  load(paste0('Normal_ESBA/Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD[j, i] = as.numeric(res[[j]]$diffBSD, units = 'secs')
      
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

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Normal_ESBA/BSD_Describe.RData')

#####################################################################################################

rm(list = ls())

powerBSC = rep(500, 20)
TrateBSC = matrix(0, 500, 20)
FrateBSC = matrix(0, 500, 20)
DifTiBSC = matrix(0, 500, 20)

for (i in 1:20)
{
  load(paste0('Normal_ESBA/Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    
    if (!is.character(res[[j]]$reBSC))
    {
      Resj = res[[j]]$reBSC$BSCL_res
      
      DifTiBSC[j, i] = as.numeric(res[[j]]$diffBSC, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSC[i] = powerBSC[i] - 1
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
      
      TrateBSC[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSC[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

save(list = c('FrateBSC', 'TrateBSC', 'powerBSC', 'DifTiBSC'), file = 'Normal_ESBA/BSC_Describe.RData')
#######################################################################################################

rm(list = ls())

powerSCQ = rep(500, 20)
TrateSCQ = matrix(0, 500, 20)
FrateSCQ = matrix(0, 500, 20)
DifTiSCQ = matrix(0, 500, 20)

for (i in 1:20)
{
  load(paste0('Normal_ESBA/Mulist', i, '.RData'))
  
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
      DifTiSCQ[j, i] = as.numeric(res[[j]]$diffSCQ, units = 'secs')
      
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

save(list = c('FrateSCQ', 'TrateSCQ', 'powerSCQ', 'DifTiSCQ'), file = 'Normal_ESBA/SCQ_Describe.RData')
###########################################################################################################

powerKSD = rep(500, 20)
TrateKSD = matrix(0, 500, 20)
FrateKSD = matrix(0, 500, 20)
DifTiKSD = matrix(0, 500, 20)

for (i in 1:20)
{ 
  load(paste0('Normal_ESBA/Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    KSDetj = res[[j]]$reKSD$KSDet_res
    
    if (length(KSDetj) != 0)
    {
      DifTiKSD[j, i] = as.numeric(res[[j]]$diffKSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, KSDetj] = 1
      }
    }
    else
    {
      powerKSD[i] = powerKSD[i] - 1
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
      
      TrateKSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateKSD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateKSD', 'TrateKSD', 'powerKSD', 'DifTiKSD'), file = 'Normal_ESBA/KSD_Describe.RData')
############################################################################################################################################################
############################################################################################################################################################

rm(list = ls())

powerBSD = rep(500, 20)
TrateBSD = matrix(0, 500, 20)
FrateBSD = matrix(0, 500, 20)
DifTiBSD = matrix(0, 500, 20)

for (i in 1:20)
{ 
  load(paste0('Normal_NSBA/Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD[j, i] = as.numeric(res[[j]]$diffBSD, units = 'secs')
      
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

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Normal_NSBA/BSD_Describe.RData')

#####################################################################################################

rm(list = ls())

powerBSC = rep(500, 20)
TrateBSC = matrix(0, 500, 20)
FrateBSC = matrix(0, 500, 20)
DifTiBSC = matrix(0, 500, 20)

for (i in 1:20)
{
  load(paste0('Normal_NSBA/Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    
    if (!is.character(res[[j]]$reBSC))
    {
      Resj = res[[j]]$reBSC$BSCL_res
      
      DifTiBSC[j, i] = as.numeric(res[[j]]$diffBSC, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSC[i] = powerBSC[i] - 1
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
      
      TrateBSC[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSC[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

save(list = c('FrateBSC', 'TrateBSC', 'powerBSC', 'DifTiBSC'), file = 'Normal_NSBA/BSC_Describe.RData')
#######################################################################################################

rm(list = ls())

powerSCQ = rep(500, 20)
TrateSCQ = matrix(0, 500, 20)
FrateSCQ = matrix(0, 500, 20)
DifTiSCQ = matrix(0, 500, 20)

for (i in 1:20)
{
  load(paste0('Normal_NSBA/Mulist', i, '.RData'))
  
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
      DifTiSCQ[j, i] = as.numeric(res[[j]]$diffSCQ, units = 'secs')
      
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

save(list = c('FrateSCQ', 'TrateSCQ', 'powerSCQ', 'DifTiSCQ'), file = 'Normal_NSBA/SCQ_Describe.RData')
###########################################################################################################

powerKSD = rep(500, 20)
TrateKSD = matrix(0, 500, 20)
FrateKSD = matrix(0, 500, 20)
DifTiKSD = matrix(0, 500, 20)

for (i in 1:20)
{ 
  load(paste0('Normal_NSBA/Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    KSDetj = res[[j]]$reKSD$KSDet_res
    
    if (length(KSDetj) != 0)
    {
      DifTiKSD[j, i] = as.numeric(res[[j]]$diffKSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, KSDetj] = 1
      }
    }
    else
    {
      powerKSD[i] = powerKSD[i] - 1
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
      
      TrateKSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateKSD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateKSD', 'TrateKSD', 'powerKSD', 'DifTiKSD'), file = 'Normal_NSBA/KSD_Describe.RData')
 