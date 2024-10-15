library(foreach)
library(doParallel)
library(BiRS)
library(QSCAN)
library(wbs)
library(Matrix)
library(mvnfast)

p = 8192
SigmaY = matrix(0, p, p)
for (j in 1:p)
{
  for (k in 1:p)
  {
    SigmaY[j, k] = 0.9^abs(j - k)
  }
}

SigmaX = 2*SigmaY

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

for (l in 1:length(mulist))
{
  source('WildBiRS.R')
  source('SeedBiRS.R')
  
  print(paste0(' mulist', l, ' :start'))
  
  p = 8192; n = 600; m = 400; nsimu = 100
  MB = 1000; alpha = 0.05 
  
  
  mu = mulist[[l]]
  
  SimuL = function(i)
  {
    set.seed(i)
    
    Xi = rmvn(n, mu, SigmaX)
    Yi = rmvn(m, rep(0, p), SigmaY)
    ############################################################################
    
    aaBSD = Sys.time()
    reBSD = SigDetDCF(Xi, Yi, foldlen = 8192, trunc = 5, MB, alpha)
    bbBSD = Sys.time()
    DifTBSD = bbBSD - aaBSD
    
    aaBSB = Sys.time()
    reBSB = SeedDetDCF(Xi, Yi, foldlen = 8192, decay = 2, trunc = 5, MB, alpha)
    bbBSB = Sys.time()
    DifTBSB = bbBSB - aaBSB
    
    aaWBB = Sys.time()
    reWBB = WildDetectBiRS(Xi, Yi, M = 500, m = 64, trunc = 5)
    bbWBB = Sys.time()
    DifTWBB = bbWBB - aaWBB
    
    aaSBB = Sys.time()
    reSBB = SeedDetectBiRS(Xi, Yi, decay = 2, m = 64, trunc = 5)
    bbSBB = Sys.time()
    DifTSBB = bbSBB - aaSBB
    
    return(list(reBSD = reBSD, DifTBSD = DifTBSD, reSBB = reSBB, DifTSBB = DifTSBB, reWBB = reWBB, DifTWBB = DifTWBB, reBSB = reBSB, diffBSB = diffBSB))
  }
  
  cl = makeCluster(50)
  registerDoParallel(cl)
  
  aa = Sys.time()
  res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN", "wbs", "Matrix")) %dopar% SimuL(i)
  bb = Sys.time()
  
  stopImplicitCluster()
  stopCluster(cl)
  
  save(list = c('res', 'mu'), file = paste0('Review-BiRS-cpts/WNS-Mulist', l, '.RData'))
  
  print(paste0(' mulist', l, ' :end; Time = ', bb - aa))
  
  rm(list = ls()[-c(which(ls() == "mulist"), which(ls() == "l"), which(ls() == "SigmaX"), which(ls() == "SigmaY"))])
  
}