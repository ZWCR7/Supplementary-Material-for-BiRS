library(foreach)
library(doParallel)
library(BiRS)
library(QSCAN)
library(SPAtest)
library(KnockoffScreen)
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
  
  set.seed(1995)
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
  source('~/Simulation_BiRS/KSDet.R')
  print(paste0(' mulist', l, ' :start'))
  
  p = 8192; n = 600; m = 400; nsimu = 500
  MB = 1000; alpha = 0.05; trunc = 5; foldlen = 4096
  Lmin = 128; Lmax = 320; skip = 32
  
  
  mu = mulist[[l]]
  
  SimuL = function(i)
  {
    set.seed(i)
    
    Xi = rmvn(n, mu, SigmaX)
    Yi = rmvn(m, rep(0, p), SigmaY)
    genotype = rbind(Xi, Yi); phenotype = c(rep(1, n), rep(0, m)); covariate = matrix(1, n + m, 1)
    genotype = as(genotype, "sparseMatrix")
    
    names = rep(0, p)
    for (j in 1:p) names[j] = paste0(j)
    
    colnames(genotype) = names
    ############################################################################
    
    #print('   BSD: start')
    aaBSD = Sys.time()
    reBSD = SigDetDCF(Xi, Yi, foldlen, trunc, MB, alpha, ReMax = 10, COMPU = 'None')
    bbBSD = Sys.time()
    
    diffBSD = bbBSD - aaBSD
    #print('   BSD: end')
    ############################################################################
    
    #print('   BSC: start')
    aaBSC = Sys.time()
    reBSC = SigDetCL(Xi, Yi, foldlen, trunc, alpha, ReMax = 10, COMPU = F)
    bbBSC = Sys.time()
    
    diffBSC = bbBSC - aaBSC
    #print('   BSC: end')
    ############################################################################
    
    #print('   SCQ: start')
    aaSCQ = Sys.time()
    reSCQ = Q_SCAN(genotype, phenotype, covariate, family = 'binomial', Lmax, Lmin)
    bbSCQ = Sys.time()
    
    diffSCQ = bbSCQ - aaSCQ
    #print('   SCQ: end')
    ############################################################################
    
    #print('   SCQ: start')
    aaKSD = Sys.time()
    result.prelim = KS.prelim(phenotype, out_type = "D")
    region.pos = c(seq(1, p, by = skip), (p+1))
    reKSD = KSDet(result.prelim, genotype, M = 5, region.pos)
    bbKSD = Sys.time()
    
    diffKSD = bbKSD - aaKSD
    #print('   SCQ: end')
    ############################################################################
    
    
    return(list(reBSD = reBSD, reBSC = reBSC, reSCQ = reSCQ, reKSD = reKSD, diffBSD = diffBSD, diffBSC = diffBSC, diffSCQ = diffSCQ, diffKSD = diffKSD))
  }
  
  cl = makeCluster(100)
  registerDoParallel(cl)
  
  aa = Sys.time()
  res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN", "KnockoffScreen", "SPAtest", "Matrix")) %dopar% SimuL(i)
  bb = Sys.time()
  
  stopImplicitCluster()
  stopCluster(cl)
  
  save(list = c('res', 'mu'), file = paste0('Normal_NSBA/Mulist', l, '.RData'))
  
  print(paste0(' mulist', l, ' :end; Time = ', bb - aa))
  
  rm(list = ls()[-c(which(ls() == "mulist"), which(ls() == "l"), which(ls() == "SigmaX"), which(ls() == "SigmaY"))])
  
}