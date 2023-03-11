library(foreach)
library(doParallel)
library(BiRS)

powerBSD_ND = rep(500, 5)
TrateBSD_ND = matrix(0, 500, 5)
FrateBSD_ND = matrix(0, 500, 5)
DifTiBSD_ND = matrix(0, 500, 5)

for (i in 1:5)
{ 
  load(paste0('Normal_ESBA/Mulist', i + 15, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD_ND[j, i] = as.numeric(res[[j]]$diffBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD_ND[i] = powerBSD_ND[i] - 1
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
      
      TrateBSD_ND[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD_ND[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD_ND', 'TrateBSD_ND', 'powerBSD_ND', 'DifTiBSD_ND'), file = 'Normal_ESBA_Decay/BSD_ND_Describe.RData')
###########################################################################################################################

powerBSD_D = rep(500, 5)
TrateBSD_D = matrix(0, 500, 5)
FrateBSD_D = matrix(0, 500, 5)
DifTiBSD_D = matrix(0, 500, 5)

for (i in 1:5)
{ 
  load(paste0('Normal_ESBA_Decay/Mulist_D', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD_D[j, i] = as.numeric(res[[j]]$diffBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD_D[i] = powerBSD_D[i] - 1
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
      
      TrateBSD_D[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD_D[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD_D', 'TrateBSD_D', 'powerBSD_D', 'DifTiBSD_D'), file = 'Normal_ESBA_Decay/BSD_D_Describe.RData')
###########################################################################################################################
###########################################################################################################################


rm(list = ls())
library(foreach)
library(doParallel)
library(BiRS)

powerBSD_ND = rep(500, 5)
TrateBSD_ND = matrix(0, 500, 5)
FrateBSD_ND = matrix(0, 500, 5)
DifTiBSD_ND = matrix(0, 500, 5)

for (i in 1:5)
{ 
  load(paste0('Normal_NSBA/Mulist', i + 15, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD_ND[j, i] = as.numeric(res[[j]]$diffBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD_ND[i] = powerBSD_ND[i] - 1
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
      
      TrateBSD_ND[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD_ND[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD_ND', 'TrateBSD_ND', 'powerBSD_ND', 'DifTiBSD_ND'), file = 'Normal_NSBA_Decay/BSD_ND_Describe.RData')
###########################################################################################################################

powerBSD_D = rep(500, 5)
TrateBSD_D = matrix(0, 500, 5)
FrateBSD_D = matrix(0, 500, 5)
DifTiBSD_D = matrix(0, 500, 5)

for (i in 1:5)
{ 
  load(paste0('Normal_NSBA_Decay/Mulist_D', i, '.RData'))
  
  mu = mul
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD_D[j, i] = as.numeric(res[[j]]$diffBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD_D[i] = powerBSD_D[i] - 1
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
      
      TrateBSD_D[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD_D[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD_D', 'TrateBSD_D', 'powerBSD_D', 'DifTiBSD_D'), file = 'Normal_NSBA_Decay/BSD_D_Describe.RData')