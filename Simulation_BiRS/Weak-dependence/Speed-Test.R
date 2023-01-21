library(foreach)
library(doParallel)
library(BiRS)

p = 8192

p = 8192
SigmaY = matrix(0, p, p)
for (j in 1:p)
{
  for (k in 1:p)
  {
    SigmaY[j, k] = (1 + abs(j - k))^(-1/4)
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

setbeta = 4

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
 
delta = 0.30

mu = rep(0, p)
theta1 = runif(length(ind1), -deltar, deltar)
theta2 = runif(length(ind2), -delta, delta)

mu[ind1] = theta1
mu[ind2] = theta2

p = 8192; n = 600; m = 400; nsimu = 100; MB = 1000; alpha = 0.05; trunc = 5; foldlen = 4096; Lmin = 128; Lmax = 320; skip = 32
 
SimuL = function(i)
{
  set.seed(i)
  Xi = rmvn(n, mu, SigmaX)
  Yi = rmvn(m, rep(0, p), SigmaY)
  
  COV = (cov(Xi)/n + cov(Yi)/m)/(1/n + 1/m)
  
  ############################################################################
  
  #print('   BSD: start')
  aaBSD = Sys.time()
  reBSD = SigDetDCF(Xi, Yi, foldlen, trunc, MB, alpha, ReMax = 10, COMPU = 'None')
  bbBSD = Sys.time()
  
  diffBSD = bbBSD - aaBSD
  #print('   BSD: end')
  ############################################################################
  
  #print('   SCD: start')
  aaSCD = Sys.time()
  reSCD = SigScanDCF(Xi, Yi, foldlen, Lmin, Lmax, skip, MB, alpha, Deal = 'Overlap', COMPU = 'None')
  bbSCD = Sys.time()
  
  diffSCD = bbSCD - aaSCD
  #print('   SCD: end')
  ############################################################################
  
  #print('   SCQ: start')
  aaSCQ = Sys.time()
  reSCQ = SigScanQ(Xi, Yi, COV, Lmax, Lmin, foldlen, skip, MB, alpha, f = 0)
  bbSCQ = Sys.time()
  
  diffSCQ = bbSCQ - aaSCQ
  
  #print('   SCQ: end')
  ###########################################################################
  
  return(list(diffBSD = diffBSD, diffSCD = diffSCD, diffSCQ = diffSCQ))
}

cl = makeCluster(100)
registerDoParallel(cl)

aa = Sys.time()
res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS")) %dopar% SimuL(i)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)

save(list = c('res'), file = paste0('Speed-Test.RData'))

T_BiRS = 0; T_Scan = 0; T_Qscan = 0

for (i in 1:nsimu)
{
  T_BiRS = T_BiRS + as.numeric(res[[i]]$diffBSD, units = "secs")
  T_Scan = T_Scan + as.numeric(res[[i]]$diffSCD, units = "secs")
  T_Qscan = T_Qscan + as.numeric(res[[i]]$diffSCQ, units = "secs")
}

T_BiRS = T_BiRS/nsimu; T_Scan = T_Scan/nsimu; T_Qscan = T_Qscan/nsimu
 