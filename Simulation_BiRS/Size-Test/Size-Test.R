library(foreach)
library(doParallel)
library(BiRS)
library(QSCAN)
library(SPAtest)
library(KnockoffScreen)
library(Matrix)
library(mvnfast)

p = 8192

source('~/Simulation_BiRS/4S_Algorithm.R')
source('~/Simulation_BiRS/LRS_Detect.R')

p = 8192; D = 64; n = 600; m = 400; nsimu = 1000
MB = 1000; alpha = 0.05; trunc = 5; foldlen = 4096
Lmin = 128; Lmax = 320; skip = 32

theta = rep(0, D)
theta1 = rep(0, D)
for (d in 1:D)
{
  theta[d] = (1 + d - 1)^(-1/4)
  theta1[d] = (1 + d)^(-1/4)
}
theta1[D] = 0; theta = sqrt(theta - theta1)

SimuL = function(i)
{
  set.seed(i)
  Zi = matrix(rmvn(n*(p + D - 1), 0, 1), n, p + D - 1)
  Wi = matrix(rmvn(m*(p + D - 1), 0, 1), m, p + D - 1)
  
  Xi = matrix(0, n, p); Yi = matrix(0, m, p)
  for (d in 1:D)
  {
    Xi = Xi + theta[d]*Zi[ ,d:(p + d - 1)]
    Yi = Yi + theta[d]*Wi[ ,d:(p + d - 1)]
  }
  
  rm(Wi); rm(Zi)
  
  genotype = rbind(Xi, Yi); phenotype = c(rep(1, n), rep(0, m)); covariate = matrix(1, n + m, 1)
  genotype = as(genotype, "sparseMatrix")
  
  names = rep(0, p)
  for (j in 1:p) names[j] = paste0(j)
  
  colnames(genotype) = names
  ############################################################################
  
  aaBSD = Sys.time()
  reBSD = SigDetDCF(Xi, Yi, foldlen, trunc, MB, alpha, ReMax = 10, COMPU = 'None')
  bbBSD = Sys.time()
  
  diffBSD = bbBSD - aaBSD
  ############################################################################
  
  aaBSC = Sys.time()
  reBSC = SigDetCL(Xi, Yi, foldlen, trunc, alpha, ReMax = 10, COMPU = F)
  bbBSC = Sys.time()
  
  diffBSC = bbBSC - aaBSC
  ############################################################################
  
  aaSCQ = Sys.time()
  reSCQ = Q_SCAN(genotype, phenotype, covariate, family = 'binomial', Lmax, Lmin)
  bbSCQ = Sys.time()
  
  diffSCQ = bbSCQ - aaSCQ
  ############################################################################
  
  aaSSS = Sys.time()
  reSSS = SSSSDetect(Xi, Yi, length.gap = 320, Lmin = 1, MB, alpha)
  bbSSS = Sys.time()
    
  diffSSS = bbSSS - aaSSS
  ############################################################################
 
  aaLRS = Sys.time()
  reLRS = ScanM(Xi, Yi, Lmin = 128, Lmax= 320, skip = 1, MB, alpha)
  bbLRS = Sys.time()
    
  diffLRS = bbLRS - aaLRS
  
  return(list(reBSD = reBSD, reBSC = reBSC, reSCQ = reSCQ, reSSS = reSSS, reLRS = reLRS, diffBSD = diffBSD, diffBSC = diffBSC, diffSCQ = diffSCQ, diffSSS = diffSSS, diffLRS = diffLRS))
}

print('Start size simulation for M dependence with equal covariance')
cl = makeCluster(50)
registerDoParallel(cl)

aa = Sys.time()
res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN")) %dopar% SimuL(i)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)

save(list = c('res'), file = 'Size-Test/Setting_MES.RData')
rm(list = ls())
#######################################################################################################################################################

p = 8192

source('~/Simulation_BiRS/4S_Algorithm.R')
source('~/Simulation_BiRS/LRS_Detect.R')

p = 8192; D = 64; n = 600; m = 400; nsimu = 1000
MB = 1000; alpha = 0.05; trunc = 5; foldlen = 4096
Lmin = 128; Lmax = 320; skip = 32

theta = rep(0, D)
theta1 = rep(0, D)
for (d in 1:D)
{
  theta[d] = (1 + d - 1)^(-1/4)
  theta1[d] = (1 + d)^(-1/4)
}
theta1[D] = 0; theta = sqrt(theta - theta1)

SimuL = function(i)
{
  set.seed(i)
  Zi = matrix(rmvn(n*(p + D - 1), 0, 1), n, p + D - 1)
  Wi = matrix(rmvn(m*(p + D - 1), 0, 1), m, p + D - 1)
  
  Xi = matrix(0, n, p); Yi = matrix(0, m, p)
  for (d in 1:D)
  {
    Xi = Xi + sqrt(2)*theta[d]*Zi[ ,d:(p + d - 1)]
    Yi = Yi + theta[d]*Wi[ ,d:(p + d - 1)]
  }
  
  rm(Wi); rm(Zi)
  
  genotype = rbind(Xi, Yi); phenotype = c(rep(1, n), rep(0, m)); covariate = matrix(1, n + m, 1)
  genotype = as(genotype, "sparseMatrix")
  
  names = rep(0, p)
  for (j in 1:p) names[j] = paste0(j)
  
  colnames(genotype) = names
  ############################################################################
  
  aaBSD = Sys.time()
  reBSD = SigDetDCF(Xi, Yi, foldlen, trunc, MB, alpha, ReMax = 10, COMPU = 'None')
  bbBSD = Sys.time()
  
  diffBSD = bbBSD - aaBSD
  ############################################################################
  
  aaBSC = Sys.time()
  reBSC = SigDetCL(Xi, Yi, foldlen, trunc, alpha, ReMax = 10, COMPU = F)
  bbBSC = Sys.time()
  
  diffBSC = bbBSC - aaBSC
  ############################################################################
  
  aaSCQ = Sys.time()
  reSCQ = Q_SCAN(genotype, phenotype, covariate, family = 'binomial', Lmax, Lmin)
  bbSCQ = Sys.time()
  
  diffSCQ = bbSCQ - aaSCQ
  ############################################################################
  
  aaSSS = Sys.time()
  reSSS = SSSSDetect(Xi, Yi, length.gap = 320, Lmin = 1, MB, alpha)
  bbSSS = Sys.time()
    
  diffSSS = bbSSS - aaSSS
  ############################################################################
 
  aaLRS = Sys.time()
  reLRS = ScanM(Xi, Yi, Lmin = 128, Lmax= 320, skip = 1, MB, alpha)
  bbLRS = Sys.time()
    
  diffLRS = bbLRS - aaLRS
  
  return(list(reBSD = reBSD, reBSC = reBSC, reSCQ = reSCQ, reSSS = reSSS, reLRS = reLRS, diffBSD = diffBSD, diffBSC = diffBSC, diffSCQ = diffSCQ, diffSSS = diffSSS, diffLRS = diffLRS))
}

print('Start size simulation for M dependence with unequal covariance')
cl = makeCluster(50)
registerDoParallel(cl)

aa = Sys.time()
res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN")) %dopar% SimuL(i)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)

save(list = c('res'), file = 'Size-Test/Setting_MNS.RData')
rm(list = ls())
################################################################################################################################################################

p = 8192

source('~/Simulation_BiRS/4S_Algorithm.R')
source('~/Simulation_BiRS/LRS_Detect.R')

SigmaY = matrix(0, p, p)
for (j in 1:p)
{
  for (k in 1:p)
  {
    SigmaY[j, k] = 0.9^abs(j - k)
  }
}

SigmaX = SigmaY
  
p = 8192; n = 600; m = 400; nsimu = 1000
MB = 1000; alpha = 0.05; trunc = 5; foldlen = 4096
Lmin = 128; Lmax = 320; skip = 32

SimuL = function(i)
{
  set.seed(i)
  
  Xi = rmvn(n, rep(0, p), SigmaX)
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
  
  aaSSS = Sys.time()
  reSSS = SSSSDetect(Xi, Yi, length.gap = 320, Lmin = 1, MB, alpha)
  bbSSS = Sys.time()
    
  diffSSS = bbSSS - aaSSS
  ############################################################################
 
  aaLRS = Sys.time()
  reLRS = ScanM(Xi, Yi, Lmin = 128, Lmax= 320, skip = 1, MB, alpha)
  bbLRS = Sys.time()
    
  diffLRS = bbLRS - aaLRS
  
  return(list(reBSD = reBSD, reBSC = reBSC, reSCQ = reSCQ, reSSS = reSSS, reLRS = reLRS, diffBSD = diffBSD, diffBSC = diffBSC, diffSCQ = diffSCQ, diffSSS = diffSSS, diffLRS = diffLRS))
}

print('Start size simulation for weak dependence with equal covariance')
cl = makeCluster(50)
registerDoParallel(cl)

aa = Sys.time()
res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN")) %dopar% SimuL(i)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)

save(list = c('res'), file = 'Size-Test/Setting_WES.RData')
rm(list = ls())
#######################################################################################################################################################
  
p = 8192

source('~/Simulation_BiRS/4S_Algorithm.R')
source('~/Simulation_BiRS/LRS_Detect.R')

SigmaY = matrix(0, p, p)
for (j in 1:p)
{
  for (k in 1:p)
  {
    SigmaY[j, k] = 0.9^abs(j - k)
  }
}

SigmaX = 2*SigmaY

p = 8192; n = 600; m = 400; nsimu = 1000
MB = 1000; alpha = 0.05; trunc = 5; foldlen = 4096
Lmin = 128; Lmax = 320; skip = 32

SimuL = function(i)
{
  set.seed(i)
  
  Xi = rmvn(n, rep(0, p), SigmaX)
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
  
  aaSSS = Sys.time()
  reSSS = SSSSDetect(Xi, Yi, length.gap = 320, Lmin = 1, MB, alpha)
  bbSSS = Sys.time()
    
  diffSSS = bbSSS - aaSSS
  ############################################################################
 
  aaLRS = Sys.time()
  reLRS = ScanM(Xi, Yi, Lmin = 128, Lmax= 320, skip = 1, MB, alpha)
  bbLRS = Sys.time()
    
  diffLRS = bbLRS - aaLRS
  
  return(list(reBSD = reBSD, reBSC = reBSC, reSCQ = reSCQ, reSSS = reSSS, reLRS = reLRS, diffBSD = diffBSD, diffBSC = diffBSC, diffSCQ = diffSCQ, diffSSS = diffSSS, diffLRS = diffLRS))
}

print('Start size simulation for weak dependence with unequal covariance')
cl = makeCluster(50)
registerDoParallel(cl)

aa = Sys.time()
res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN")) %dopar% SimuL(i)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)

save(list = c('res'), file = 'Size-Test/Setting_WNS.RData')
rm(list = ls())
#######################################################################################################################################

 
