rm(list = ls())

powerBSD = rep(100, 20)
TrateBSD = matrix(0, 100, 20)
FrateBSD = matrix(0, 100, 20)
DifTiBSD = matrix(0, 100, 20)

powerSBB = rep(100, 20)
TrateSBB = matrix(0, 100, 20)
FrateSBB = matrix(0, 100, 20)
DifTiSBB = matrix(0, 100, 20)

powerWBB = rep(100, 20)
TrateWBB = matrix(0, 100, 20)
FrateWBB = matrix(0, 100, 20)
DifTiWBB = matrix(0, 100, 20)

powerBSB = rep(100, 20)
TrateBSB = matrix(0, 100, 20)
FrateBSB = matrix(0, 100, 20)
DifTiBSB = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-BiRS-cpts/MES-Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  
  EstBSD = matrix(0, nsimu, p)
  EstSBB = matrix(0, nsimu, p)
  EstWBB = matrix(0, nsimu, p)
  EstBSB = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      for (k in 1:nrow(Resj))
      {
        EstBSD[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }

    if (!is.character(res[[j]]$reBSB))
    {
      Resj = res[[j]]$reBSB$SBDCF_res

      for (k in 1:nrow(Resj))
      {
        EstBSB[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSB[i] = powerBSB[i] - 1
    }
    
    EstWBB[j, ] = res[[j]]$reWBB
    EstSBB[j, ] = res[[j]]$reSBB
    
    DifTiBSD[j, i] = as.numeric(res[[j]]$DifTBSD, units = 'secs')
    DifTiBSB[j, i] = as.numeric(res[[j]]$DifTBSB, units = 'secs')
    DifTiWBB[j, i] = as.numeric(res[[j]]$DifTWBB, units = 'secs')
    DifTiSBB[j, i] = as.numeric(res[[j]]$DifTSBB, units = 'secs')
  }
  
  mu = mul
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstBSD[l, ] != 0)) != 0)
    {
      Diffl = EstBSD[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstBSD[l, ] != 0))
    }
    
    if (length(which(EstSBB[l, ] != 0)) != 0)
    {
      Diffl = EstSBB[l, ] - true.ind
      
      TrateSBB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSBB[l, i] = length(which(Diffl == 3))/length(which(EstSBB[l, ] != 0))
    }
    
    if (length(which(EstWBB[l, ] != 0)) != 0)
    {
      Diffl = EstWBB[l, ] - true.ind
      
      TrateWBB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateWBB[l, i] = length(which(Diffl == 3))/length(which(EstWBB[l, ] != 0))
    }
    
    if (length(which(EstBSB[l, ] != 0)) != 0)
    {
      Diffl = EstBSB[l, ] - true.ind

      TrateBSB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSB[l, i] = length(which(Diffl == 3))/length(which(EstBSB[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Review-BiRS-cpts/MES-BSD_Describe.RData')
save(list = c('FrateBSB', 'TrateBSB', 'powerBSB', 'DifTiBSB'), file = 'Review-BiRS-cpts/MES-BSB_Describe.RData')
save(list = c('FrateSBB', 'TrateSBB', 'powerSBB', 'DifTiSBB'), file = 'Review-BiRS-cpts/MES-SBB_Describe.RData')
save(list = c('FrateWBB', 'TrateWBB', 'powerWBB', 'DifTiWBB'), file = 'Review-BiRS-cpts/MES-WBB_Describe.RData')
###############################################################################################################################################

rm(list = ls())

powerBSD = rep(100, 20)
TrateBSD = matrix(0, 100, 20)
FrateBSD = matrix(0, 100, 20)
DifTiBSD = matrix(0, 100, 20)

powerSBB = rep(100, 20)
TrateSBB = matrix(0, 100, 20)
FrateSBB = matrix(0, 100, 20)
DifTiSBB = matrix(0, 100, 20)

powerWBB = rep(100, 20)
TrateWBB = matrix(0, 100, 20)
FrateWBB = matrix(0, 100, 20)
DifTiWBB = matrix(0, 100, 20)

powerBSB = rep(100, 20)
TrateBSB = matrix(0, 100, 20)
FrateBSB = matrix(0, 100, 20)
DifTiBSB = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-BiRS-cpts/MNS-Mulist', i, '.RData'))

  nsimu = length(res); p = 8192
  
  EstBSD = matrix(0, nsimu, p)
  EstSBB = matrix(0, nsimu, p)
  EstWBB = matrix(0, nsimu, p)
  EstBSB = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      for (k in 1:nrow(Resj))
      {
        EstBSD[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }

    if (!is.character(res[[j]]$reBSB))
    {
      Resj = res[[j]]$reBSB$SBDCF_res

      for (k in 1:nrow(Resj))
      {
        EstBSB[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSB[i] = powerBSB[i] - 1
    }
    
    EstWBB[j, ] = res[[j]]$reWBB
    EstSBB[j, ] = res[[j]]$reSBB
    
    DifTiBSD[j, i] = as.numeric(res[[j]]$DifTBSD, units = 'secs')
    DifTiBSB[j, i] = as.numeric(res[[j]]$DifTBSB, units = 'secs')
    DifTiWBB[j, i] = as.numeric(res[[j]]$DifTWBB, units = 'secs')
    DifTiSBB[j, i] = as.numeric(res[[j]]$DifTSBB, units = 'secs')
  }
  
  mu = mul
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstBSD[l, ] != 0)) != 0)
    {
      Diffl = EstBSD[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstBSD[l, ] != 0))
    }
    
    if (length(which(EstSBB[l, ] != 0)) != 0)
    {
      Diffl = EstSBB[l, ] - true.ind
      
      TrateSBB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSBB[l, i] = length(which(Diffl == 3))/length(which(EstSBB[l, ] != 0))
    }
    
    if (length(which(EstWBB[l, ] != 0)) != 0)
    {
      Diffl = EstWBB[l, ] - true.ind
      
      TrateWBB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateWBB[l, i] = length(which(Diffl == 3))/length(which(EstWBB[l, ] != 0))
    }
    
    if (length(which(EstBSB[l, ] != 0)) != 0)
    {
      Diffl = EstBSB[l, ] - true.ind

      TrateBSB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSB[l, i] = length(which(Diffl == 3))/length(which(EstBSB[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Review-BiRS-cpts/MNS-BSD_Describe.RData')
save(list = c('FrateBSB', 'TrateBSB', 'powerBSB', 'DifTiBSB'), file = 'Review-BiRS-cpts/MNS-BSB_Describe.RData')
save(list = c('FrateSBB', 'TrateSBB', 'powerSBB', 'DifTiSBB'), file = 'Review-BiRS-cpts/MNS-SBB_Describe.RData')
save(list = c('FrateWBB', 'TrateWBB', 'powerWBB', 'DifTiWBB'), file = 'Review-BiRS-cpts/MNS-WBB_Describe.RData')
###############################################################################################################################################
rm(list = ls())

powerBSD = rep(100, 20)
TrateBSD = matrix(0, 100, 20)
FrateBSD = matrix(0, 100, 20)
DifTiBSD = matrix(0, 100, 20)

powerSBB = rep(100, 20)
TrateSBB = matrix(0, 100, 20)
FrateSBB = matrix(0, 100, 20)
DifTiSBB = matrix(0, 100, 20)

powerWBB = rep(100, 20)
TrateWBB = matrix(0, 100, 20)
FrateWBB = matrix(0, 100, 20)
DifTiWBB = matrix(0, 100, 20)

powerBSB = rep(100, 20)
TrateBSB = matrix(0, 100, 20)
FrateBSB = matrix(0, 100, 20)
DifTiBSB = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-BiRS-cpts/WES-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  
  EstBSD = matrix(0, nsimu, p)
  EstSBB = matrix(0, nsimu, p)
  EstWBB = matrix(0, nsimu, p)
  EstBSB = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      for (k in 1:nrow(Resj))
      {
        EstBSD[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }

    if (!is.character(res[[j]]$reBSB))
    {
      Resj = res[[j]]$reBSB$SBDCF_res

      for (k in 1:nrow(Resj))
      {
        EstBSB[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSB[i] = powerBSB[i] - 1
    }
    
    EstWBB[j, ] = res[[j]]$reWBB
    EstSBB[j, ] = res[[j]]$reSBB
    
    DifTiBSD[j, i] = as.numeric(res[[j]]$DifTBSD, units = 'secs')
    DifTiBSB[j, i] = as.numeric(res[[j]]$DifTBSB, units = 'secs')
    DifTiWBB[j, i] = as.numeric(res[[j]]$DifTWBB, units = 'secs')
    DifTiSBB[j, i] = as.numeric(res[[j]]$DifTSBB, units = 'secs')
  }
  
  #mu = mul
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstBSD[l, ] != 0)) != 0)
    {
      Diffl = EstBSD[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstBSD[l, ] != 0))
    }
    
    if (length(which(EstSBB[l, ] != 0)) != 0)
    {
      Diffl = EstSBB[l, ] - true.ind
      
      TrateSBB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSBB[l, i] = length(which(Diffl == 3))/length(which(EstSBB[l, ] != 0))
    }
    
    if (length(which(EstWBB[l, ] != 0)) != 0)
    {
      Diffl = EstWBB[l, ] - true.ind
      
      TrateWBB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateWBB[l, i] = length(which(Diffl == 3))/length(which(EstWBB[l, ] != 0))
    }
    
    if (length(which(EstBSB[l, ] != 0)) != 0)
    {
      Diffl = EstBSB[l, ] - true.ind

      TrateBSB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSB[l, i] = length(which(Diffl == 3))/length(which(EstBSB[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Review-BiRS-cpts/WES-BSD_Describe.RData')
save(list = c('FrateBSB', 'TrateBSB', 'powerBSB', 'DifTiBSB'), file = 'Review-BiRS-cpts/WES-BSB_Describe.RData')
save(list = c('FrateSBB', 'TrateSBB', 'powerSBB', 'DifTiSBB'), file = 'Review-BiRS-cpts/WES-SBB_Describe.RData')
save(list = c('FrateWBB', 'TrateWBB', 'powerWBB', 'DifTiWBB'), file = 'Review-BiRS-cpts/WES-WBB_Describe.RData')
###############################################################################################################################################

rm(list = ls())

powerBSD = rep(100, 20)
TrateBSD = matrix(0, 100, 20)
FrateBSD = matrix(0, 100, 20)
DifTiBSD = matrix(0, 100, 20)

powerSBB = rep(100, 20)
TrateSBB = matrix(0, 100, 20)
FrateSBB = matrix(0, 100, 20)
DifTiSBB = matrix(0, 100, 20)

powerWBB = rep(100, 20)
TrateWBB = matrix(0, 100, 20)
FrateWBB = matrix(0, 100, 20)
DifTiWBB = matrix(0, 100, 20)

powerBSB = rep(100, 20)
TrateBSB = matrix(0, 100, 20)
FrateBSB = matrix(0, 100, 20)
DifTiBSB = matrix(0, 100, 20)

for (i in 1:20)
{ 
  load(paste0('Review-BiRS-cpts/WNS-Mulist', i, '.RData'))
  
  #mu = mul
  nsimu = length(res); p = 8192
  
  EstBSD = matrix(0, nsimu, p)
  EstSBB = matrix(0, nsimu, p)
  EstWBB = matrix(0, nsimu, p)
  EstBSB = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res

      for (k in 1:nrow(Resj))
      {
        EstBSD[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }

    if (!is.character(res[[j]]$reBSB))
    {
      Resj = res[[j]]$reBSB$SBDCF_res

      for (k in 1:nrow(Resj))
      {
        EstBSB[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSB[i] = powerBSB[i] - 1
    }
    
    EstWBB[j, ] = res[[j]]$reWBB
    EstSBB[j, ] = res[[j]]$reSBB
    
    DifTiBSD[j, i] = as.numeric(res[[j]]$DifTBSD, units = 'secs')
    DifTiBSB[j, i] = as.numeric(res[[j]]$DifTBSB, units = 'secs')
    DifTiWBB[j, i] = as.numeric(res[[j]]$DifTWBB, units = 'secs')
    DifTiSBB[j, i] = as.numeric(res[[j]]$DifTSBB, units = 'secs')
  }
  
  #mu = mul
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstBSD[l, ] != 0)) != 0)
    {
      Diffl = EstBSD[l, ] - true.ind

      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstBSD[l, ] != 0))
    }

    if (length(which(EstSBB[l, ] != 0)) != 0)
    {
      Diffl = EstSBB[l, ] - true.ind
      
      TrateSBB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSBB[l, i] = length(which(Diffl == 3))/length(which(EstSBB[l, ] != 0))
    }
    
    if (length(which(EstWBB[l, ] != 0)) != 0)
    {
      Diffl = EstWBB[l, ] - true.ind
      
      TrateWBB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateWBB[l, i] = length(which(Diffl == 3))/length(which(EstWBB[l, ] != 0))
    }
    
    if (length(which(EstBSB[l, ] != 0)) != 0)
    {
      Diffl = EstBSB[l, ] - true.ind

      TrateBSB[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSB[l, i] = length(which(Diffl == 3))/length(which(EstBSB[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'Review-BiRS-cpts/WNS-BSD_Describe.RData')
save(list = c('FrateBSB', 'TrateBSB', 'powerBSB', 'DifTiBSB'), file = 'Review-BiRS-cpts/WNS-BSB_Describe.RData')
save(list = c('FrateSBB', 'TrateSBB', 'powerSBB', 'DifTiSBB'), file = 'Review-BiRS-cpts/WNS-SBB_Describe.RData')
save(list = c('FrateWBB', 'TrateWBB', 'powerWBB', 'DifTiWBB'), file = 'Review-BiRS-cpts/WNS-WBB_Describe.RData')