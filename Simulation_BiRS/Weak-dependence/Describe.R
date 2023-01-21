p = 8192

dense = 4
deltar = 0.01
deltav = c(0.20, 0.25, 0.30, 0.35, 0.40)

mulist = list()

set.seed(4)
lensam = sample(seq(128, 320, by = 32), size = 4, replace = F)

set.seed(7)
intersam = sample(seq(1, 2*dense, by = 2), size = dense, replace = F)

startp = rep(0, dense)
for (i in 1:dense)
{
  left = (intersam[i] - 1)*(p/2/dense)
  right = intersam[i]*(p/2/dense)
  startp[i] = sample(left:right, size = 1, replace = F)
}

endp = startp + lensam - 1

for (setbeta in 1:dense)
{
  set.seed(1998)
  ind1 = NULL
  ind2 = NULL
  for (j in 1:setbeta)
  {
    pos1 = startp[j]; pos2 = endp[j]
    ind1 = c(ind1, pos1:pos2)
    ind2 = c(ind2, sample(pos1:pos2, size = floor((pos2 - pos1 + 1)/4), replace = F))
  }
  
  set.seed(1024)
  for (setdelta in 1:length(deltav))
  {
    delta = deltav[setdelta]
    
    mu = rep(0, p)
    theta1 = runif(length(ind1), -deltar, deltar)
    theta2 = runif(length(ind2), -delta, delta)
    
    mu[ind1] = theta1
    mu[ind2] = theta2
    
    mulist = c(mulist, list(mu))
  }
}


powerBSD = rep(1000, 20)
TrateBSD = matrix(0, 1000, 20)
FrateBSD = matrix(0, 1000, 20)
DifTiBSD = matrix(0, 1000, 20)

for (i in 1:20)
{ 
  mu = mulist[[i]]
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

# rm(list = ls()[-c(which(ls() == "mulist"))])
# 
# powerSCD = rep(1000, 20)
# TrateSCD = matrix(0, 1000, 20)
# FrateSCD = matrix(0, 1000, 20)
# DifTiSCD = matrix(0, 1000, 20)
# 
# for (i in 1:20)
# {
#   mu = mulist[[i]]
#   load(paste0('Normal_ESBA/Mulist', i, '.RData'))
#   
#   nsimu = length(res); p = 8192
#   EstR = matrix(0, nsimu, p)
#   
#   for (j in 1:length(res))
#   {
#     
#     if (!is.character(res[[j]]$reSCD))
#     {
#       Resj = res[[j]]$reSCD$BSDCF_res
#       
#       DifTiSCD[j, i] = as.numeric(res[[j]]$diffSCD, units = 'secs')
#       
#       for (k in 1:nrow(Resj))
#       {
#         EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
#       }
#     }
#     else
#     {
#       powerSCD[i] = powerSCD[i] - 1
#     }
#   }
#   
#   true.ind = rep(0, p)
#   true.ind[which(mu != 0)] = 1
#   true.ind[which(mu == 0)] = -2
#   
#   for (l in 1:nsimu)
#   {
#     if (length(which(EstR[l, ] != 0)) != 0)
#     {
#       Diffl = EstR[l, ] - true.ind
#       
#       TrateSCD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
#       FrateSCD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
#     }
#   }
#   
# }
# 
# save(list = c('FrateSCD', 'TrateSCD', 'powerSCD', 'DifTiSCD'), file = 'Normal_ESBA/SCD_Describe.RData')
#######################################################################################################

rm(list = ls()[-c(which(ls() == "mulist"))])

powerBSC = rep(1000, 20)
TrateBSC = matrix(0, 1000, 20)
FrateBSC = matrix(0, 1000, 20)
DifTiBSC = matrix(0, 1000, 20)

for (i in 1:20)
{
  mu = mulist[[i]]
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

# rm(list = ls()[-c(which(ls() == "mulist"))])
# 
# powerSCC = rep(1000, 20)
# TrateSCC = matrix(0, 1000, 20)
# FrateSCC = matrix(0, 1000, 20)
# DifTiSCC = matrix(0, 1000, 20)
# 
# for (i in 1:20)
# {
#   mu = mulist[[i]]
#   load(paste0('Normal_ESBA/Mulist', i, '.RData'))
#   
#   nsimu = length(res); p = 8192
#   EstR = matrix(0, nsimu, p)
#   
#   for (j in 1:length(res))
#   {
#     
#     if (!is.character(res[[j]]$reSCC))
#     {
#       Resj = res[[j]]$reSCC$SCCL_res
#       
#       DifTiSCC[j, i] = as.numeric(res[[j]]$diffSCC, units = 'secs')
#       
#       for (k in 1:nrow(Resj))
#       {
#         EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
#       }
#     }
#     else
#     {
#       powerSCC[i] = powerSCC[i] - 1
#     }
#   }
#   
#   true.ind = rep(0, p)
#   true.ind[which(mu != 0)] = 1
#   true.ind[which(mu == 0)] = -2
#   
#   for (l in 1:nsimu)
#   {
#     if (length(which(EstR[l, ] != 0)) != 0)
#     {
#       Diffl = EstR[l, ] - true.ind
#       
#       TrateSCC[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
#       FrateSCC[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
#     }
#   }
#   
# }
# 
# save(list = c('FrateSCC', 'TrateSCC', 'powerSCC', 'DifTiSCC'), file = 'Normal_ESBA/SCC_Describe.RData')
#######################################################################################################

rm(list = ls()[-c(which(ls() == "mulist"))])

powerSCQ = rep(1000, 20)
TrateSCQ = matrix(0, 1000, 20)
FrateSCQ = matrix(0, 1000, 20)
DifTiSCQ = matrix(0, 1000, 20)

for (i in 1:20)
{
  mu = mulist[[i]]
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

# rm(list = ls()[-c(which(ls() == "mulist"))])
# 
# powerSCM = rep(1000, 20)
# TrateSCM = matrix(0, 1000, 20)
# FrateSCM = matrix(0, 1000, 20)
# DifTiSCM = matrix(0, 1000, 20)
# 
# for (i in 1:20)
# {
#   mu = mulist[[i]]
#   load(paste0('Normal_ESBA/Mulist', i, '.RData'))
#   
#   nsimu = length(res); p = 8192
#   EstR = matrix(0, nsimu, p)
#   
#   for (j in 1:length(res))
#   {
#     Resj = as.matrix(res[[j]]$reSCM$SCAN_res)
#     
#     if (ncol(Resj) == 1) Resj = t(Resj)
#     
#     if (sum(Resj) == 1)
#     {
#       powerSCM[i] = powerSCM[i] - 1
#     }
#     else
#     {
#       DifTiSCM[j, i] = as.numeric(res[[j]]$diffSCM, units = 'secs')
#       
#       for (k in 1:nrow(Resj))
#       {
#         EstR[j, (Resj[k, 2]:Resj[k, 3])] = 1
#       }
#     }
#   }
#   
#   true.ind = rep(0, p)
#   true.ind[which(mu != 0)] = 1
#   true.ind[which(mu == 0)] = -2
#   
#   for (l in 1:nsimu)
#   {
#     if (length(which(EstR[l, ] != 0)) != 0)
#     {
#       Diffl = EstR[l, ] - true.ind
#       
#       TrateSCM[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
#       FrateSCM[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
#     }
#   }
#   
# }
# 
# save(list = c('FrateSCM', 'TrateSCM', 'powerSCM', 'DifTiSCM'), file = 'Normal_ESBA/SCM_Describe.RData')
####################################################################################################################################

rm(list = ls())
p = 8192

dense = 4
deltar = 0.01
deltav = c(0.20, 0.25, 0.30, 0.35, 0.40)

mulist = list()

set.seed(4)
lensam = sample(seq(128, 320, by = 32), size = 4, replace = F)

set.seed(7)
intersam = sample(seq(1, 2*dense, by = 2), size = dense, replace = F)

startp = rep(0, dense)
for (i in 1:dense)
{
  left = (intersam[i] - 1)*(p/2/dense)
  right = intersam[i]*(p/2/dense)
  startp[i] = sample(left:right, size = 1, replace = F)
}

endp = startp + lensam - 1

for (setbeta in 1:dense)
{
  set.seed(1998)
  ind1 = NULL
  ind2 = NULL
  for (j in 1:setbeta)
  {
    pos1 = startp[j]; pos2 = endp[j]
    ind1 = c(ind1, pos1:pos2)
    ind2 = c(ind2, sample(pos1:pos2, size = floor((pos2 - pos1 + 1)/4), replace = F))
  }
  
  set.seed(1024)
  for (setdelta in 1:length(deltav))
  {
    delta = deltav[setdelta]
    
    mu = rep(0, p)
    theta1 = runif(length(ind1), -deltar, deltar)
    theta2 = runif(length(ind2), -delta, delta)
    
    mu[ind1] = theta1
    mu[ind2] = theta2
    
    mulist = c(mulist, list(mu))
  }
}


powerBSD = rep(1000, 20)
TrateBSD = matrix(0, 1000, 20)
FrateBSD = matrix(0, 1000, 20)
DifTiBSD = matrix(0, 1000, 20)

for (i in 1:20)
{ 
  mu = mulist[[i]]
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

# rm(list = ls()[-c(which(ls() == "mulist"))])
# 
# powerSCD = rep(1000, 20)
# TrateSCD = matrix(0, 1000, 20)
# FrateSCD = matrix(0, 1000, 20)
# DifTiSCD = matrix(0, 1000, 20)
# 
# for (i in 1:20)
# {
#   mu = mulist[[i]]
#   load(paste0('Normal_NSBA/Mulist', i, '.RData'))
#   
#   nsimu = length(res); p = 8192
#   EstR = matrix(0, nsimu, p)
#   
#   for (j in 1:length(res))
#   {
#     
#     if (!is.character(res[[j]]$reSCD))
#     {
#       Resj = res[[j]]$reSCD$BSDCF_res
#       
#       DifTiSCD[j, i] = as.numeric(res[[j]]$diffSCD, units = 'secs')
#       
#       for (k in 1:nrow(Resj))
#       {
#         EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
#       }
#     }
#     else
#     {
#       powerSCD[i] = powerSCD[i] - 1
#     }
#   }
#   
#   true.ind = rep(0, p)
#   true.ind[which(mu != 0)] = 1
#   true.ind[which(mu == 0)] = -2
#   
#   for (l in 1:nsimu)
#   {
#     if (length(which(EstR[l, ] != 0)) != 0)
#     {
#       Diffl = EstR[l, ] - true.ind
#       
#       TrateSCD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
#       FrateSCD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
#     }
#   }
#   
# }
# 
# save(list = c('FrateSCD', 'TrateSCD', 'powerSCD', 'DifTiSCD'), file = 'Normal_NSBA/SCD_Describe.RData')
#######################################################################################################

rm(list = ls()[-c(which(ls() == "mulist"))])

powerBSC = rep(1000, 20)
TrateBSC = matrix(0, 1000, 20)
FrateBSC = matrix(0, 1000, 20)
DifTiBSC = matrix(0, 1000, 20)

for (i in 1:20)
{
  mu = mulist[[i]]
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

# rm(list = ls()[-c(which(ls() == "mulist"))])
# 
# powerSCC = rep(1000, 20)
# TrateSCC = matrix(0, 1000, 20)
# FrateSCC = matrix(0, 1000, 20)
# DifTiSCC = matrix(0, 1000, 20)
# 
# for (i in 1:20)
# {
#   mu = mulist[[i]]
#   load(paste0('Normal_NSBA/Mulist', i, '.RData'))
#   
#   nsimu = length(res); p = 8192
#   EstR = matrix(0, nsimu, p)
#   
#   for (j in 1:length(res))
#   {
#     
#     if (!is.character(res[[j]]$reSCC))
#     {
#       Resj = res[[j]]$reSCC$SCCL_res
#       
#       DifTiSCC[j, i] = as.numeric(res[[j]]$diffSCC, units = 'secs')
#       
#       for (k in 1:nrow(Resj))
#       {
#         EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
#       }
#     }
#     else
#     {
#       powerSCC[i] = powerSCC[i] - 1
#     }
#   }
#   
#   true.ind = rep(0, p)
#   true.ind[which(mu != 0)] = 1
#   true.ind[which(mu == 0)] = -2
#   
#   for (l in 1:nsimu)
#   {
#     if (length(which(EstR[l, ] != 0)) != 0)
#     {
#       Diffl = EstR[l, ] - true.ind
#       
#       TrateSCC[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
#       FrateSCC[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
#     }
#   }
#   
# }
# 
# save(list = c('FrateSCC', 'TrateSCC', 'powerSCC', 'DifTiSCC'), file = 'Normal_NSBA/SCC_Describe.RData')
#######################################################################################################

rm(list = ls()[-c(which(ls() == "mulist"))])

powerSCQ = rep(1000, 20)
TrateSCQ = matrix(0, 1000, 20)
FrateSCQ = matrix(0, 1000, 20)
DifTiSCQ = matrix(0, 1000, 20)

for (i in 1:20)
{
  mu = mulist[[i]]
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

# rm(list = ls()[-c(which(ls() == "mulist"))])
# 
# powerSCM = rep(1000, 20)
# TrateSCM = matrix(0, 1000, 20)
# FrateSCM = matrix(0, 1000, 20)
# DifTiSCM = matrix(0, 1000, 20)
# 
# for (i in 1:20)
# {
#   mu = mulist[[i]]
#   load(paste0('Normal_NSBA/Mulist', i, '.RData'))
#   
#   nsimu = length(res); p = 8192
#   EstR = matrix(0, nsimu, p)
#   
#   for (j in 1:length(res))
#   {
#     Resj = as.matrix(res[[j]]$reSCM$SCAN_res)
#     
#     if (ncol(Resj) == 1) Resj = t(Resj)
#     
#     if (sum(Resj) == 1)
#     {
#       powerSCM[i] = powerSCM[i] - 1
#     }
#     else
#     {
#       DifTiSCM[j, i] = as.numeric(res[[j]]$diffSCM, units = 'secs')
#       
#       for (k in 1:nrow(Resj))
#       {
#         EstR[j, (Resj[k, 2]:Resj[k, 3])] = 1
#       }
#     }
#   }
#   
#   true.ind = rep(0, p)
#   true.ind[which(mu != 0)] = 1
#   true.ind[which(mu == 0)] = -2
#   
#   for (l in 1:nsimu)
#   {
#     if (length(which(EstR[l, ] != 0)) != 0)
#     {
#       Diffl = EstR[l, ] - true.ind
#       
#       TrateSCM[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
#       FrateSCM[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
#     }
#   }
#   
# }
# 
# save(list = c('FrateSCM', 'TrateSCM', 'powerSCM', 'DifTiSCM'), file = 'Normal_NSBA/SCM_Describe.RData')