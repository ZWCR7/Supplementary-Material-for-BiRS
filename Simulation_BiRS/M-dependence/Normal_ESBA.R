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

for (l in 1:length(mulist))
{
  print(paste0(' mulist', l, ' :start'))
  
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
  
  mul = mulist[[l]]
  mu = matrix(rep(mulist[[l]], n), n, p, byrow = T)
  
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
    
    Xi = mu + Xi
    #COV = (cov(Xi)/n + cov(Yi)/m)/(1/n + 1/m)
    
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
    
    #print('   SCD: start')
    # aaSCD = Sys.time()
    # reSCD = SigScanDCF(Xi, Yi, foldlen, Lmin, Lmax, skip, MB, alpha, Deal = 'Overlap', COMPU = 'None')
    # bbSCD = Sys.time()
    # 
    # diffSCD = bbSCD - aaSCD
    #print('   SCD: end')
    ############################################################################
    
    #print('   SCC: start')
    # aaSCC = Sys.time()
    # reSCC = SigScanCL(Xi, Yi, foldlen, Lmin, Lmax, skip, alpha, Deal = 'Overlap', COMPU = F)
    # bbSCC = Sys.time()
    # 
    # diffSCC = bbSCC - aaSCC
    #print('   SCC: end')
    ############################################################################
    
    #print('   SCQ: start')
    # aaSCQ = Sys.time()
    # reSCQ = SigScanQ(Xi, Yi, COV, Lmax, Lmin, foldlen, skip, MB, alpha, f = 0)
    # bbSCQ = Sys.time()
    # 
    # diffSCQ = bbSCQ - aaSCQ
    #print('   SCQ: end')
    ############################################################################
    
    #print('   SCM: start')
    # aaSCM = Sys.time()
    # reSCM = SigScanM(Xi, Yi, Lmax, Lmin, foldlen, skip, MB, alpha, f = 0)
    # bbSCM = Sys.time()
    # 
    # diffSCM = bbSCM - aaSCM
    #print('   SCM: end')
    # return(list(reBSD = reBSD, reBSC = reBSC, reSCD = reSCD, reSCC = reSCC, reSCQ = reSCQ, reSCM = reSCM, 
    #             diffBSD = diffBSD, diffBSC = diffBSC, diffSCD = diffSCD, diffSCC = diffSCC, diffSCQ = diffSCQ, 
    #             diffSCM = diffSCM))
    
    #return(list(reBSD = reBSD, reBSC = reBSC, reSCQ = reSCQ, reSCM = reSCM, diffBSD = diffBSD, diffBSC = diffBSC, diffSCQ = diffSCQ, diffSCM = diffSCM))
    
    return(list(reBSD = reBSD, reBSC = reBSC, diffBSD = diffBSD, diffBSC = diffBSC))
  }
  
  cl = makeCluster(100)
  registerDoParallel(cl)
  
  aa = Sys.time()
  res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS")) %dopar% SimuL(i)
  bb = Sys.time()
  
  stopImplicitCluster()
  stopCluster(cl)
  
  save(list = c('res', 'mul'), file = paste0('Mulist', l, '.RData'))
  
  print(paste0(' mulist', l, ' :end; Time = ', bb - aa))
  
  rm(list = ls()[-c(which(ls() == "mulist"), which(ls() == "l"))])
  
}