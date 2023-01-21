library(foreach)
library(doParallel)
library(BiRS)

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


powerBSD_ND = rep(1000, 5)
TrateBSD_ND = matrix(0, 1000, 5)
FrateBSD_ND = matrix(0, 1000, 5)
DifTiBSD_ND = matrix(0, 1000, 5)

for (i in 1:5)
{ 
  mu = mulist[[i + 15]]
  load(paste0('Normal_ESBA_Decay/Mulist_ND', i, '.RData'))
  
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

rm(list = ls())

p = 8192

dense = 4
deltar = 0.05
deltav = c(0.20, 0.25, 0.30, 0.35, 0.40)

mulist = list()

set.seed(4)
lensam = sample(seq(128, 288, by = 32), size = 4, replace = F)
sum_len = c(0, lensam[1], lensam[1] + lensam[2], lensam[1] + lensam[2] + lensam[3])

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

n = 600; m = 400

set.seed(1998)
ind1 = NULL; ind2list = list()
for (j in 1:dense)
{
  pos1 = startp[j]; pos2 = endp[j]
  indsel = sample(pos1:pos2, size = floor((pos2 - pos1 + 1)/3), replace = F)
  ind2list = c(ind2list, list(indsel))
}

set.seed(1024)
for (setdelta in 1:length(deltav))
{
  delta = deltav[setdelta]
  
  mu = rep(0, p)
  for (j in 1:dense)
  {
    pos1 = startp[j]; pos2 = endp[j]
    ind1 = pos1:pos2
    ind2 = ind2list[[j]]
    
    decay_rate = 100*((log(p*n)/n)^(1/2) - (log((p - sum_len[j])*n)/n)^(1/2))
    
    theta1 = runif(length(ind1), -deltar, deltar)
    theta2 = runif(length(ind2), -delta + decay_rate, delta - decay_rate)
    
    mu[ind1] = theta1
    mu[ind2] = theta2
  }
  
  mulist = c(mulist, list(mu))
}


powerBSD_D = rep(1000, 5)
TrateBSD_D = matrix(0, 1000, 5)
FrateBSD_D = matrix(0, 1000, 5)
DifTiBSD_D = matrix(0, 1000, 5)

for (i in 1:5)
{ 
  mu = mulist[[i]]
  load(paste0('Normal_ESBA_Decay/Mulist_D', i, '.RData'))
  
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

rm(list = ls())
library(foreach)
library(doParallel)
library(BiRS)

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


powerBSD_ND = rep(1000, 5)
TrateBSD_ND = matrix(0, 1000, 5)
FrateBSD_ND = matrix(0, 1000, 5)
DifTiBSD_ND = matrix(0, 1000, 5)

for (i in 1:5)
{ 
  mu = mulist[[i + 15]]
  load(paste0('Normal_NSBA_Decay/Mulist_ND', i, '.RData'))
  
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

rm(list = ls())

p = 8192

dense = 4
deltar = 0.05
deltav = c(0.20, 0.25, 0.30, 0.35, 0.40)

mulist = list()

set.seed(4)
lensam = sample(seq(128, 288, by = 32), size = 4, replace = F)
sum_len = c(0, lensam[1], lensam[1] + lensam[2], lensam[1] + lensam[2] + lensam[3])

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

n = 600; m = 400

set.seed(1998)
ind1 = NULL; ind2list = list()
for (j in 1:dense)
{
  pos1 = startp[j]; pos2 = endp[j]
  indsel = sample(pos1:pos2, size = floor((pos2 - pos1 + 1)/3), replace = F)
  ind2list = c(ind2list, list(indsel))
}

set.seed(1024)
for (setdelta in 1:length(deltav))
{
  delta = deltav[setdelta]
  
  mu = rep(0, p)
  for (j in 1:dense)
  {
    pos1 = startp[j]; pos2 = endp[j]
    ind1 = pos1:pos2
    ind2 = ind2list[[j]]
    
    decay_rate = 100*((log(p*n)/n)^(1/2) - (log((p - sum_len[j])*n)/n)^(1/2))
    
    theta1 = runif(length(ind1), -deltar, deltar)
    theta2 = runif(length(ind2), -delta + decay_rate, delta - decay_rate)
    
    mu[ind1] = theta1
    mu[ind2] = theta2
  }
  
  mulist = c(mulist, list(mu))
}


powerBSD_D = rep(1000, 5)
TrateBSD_D = matrix(0, 1000, 5)
FrateBSD_D = matrix(0, 1000, 5)
DifTiBSD_D = matrix(0, 1000, 5)

for (i in 1:5)
{ 
  mu = mulist[[i]]
  load(paste0('Normal_NSBA_Decay/Mulist_D', i, '.RData'))
  
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

