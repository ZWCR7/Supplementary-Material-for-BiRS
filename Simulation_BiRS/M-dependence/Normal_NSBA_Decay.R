library(foreach)
library(doParallel)
library(BiRS)

p = 8192

dense = 4
deltar = 0.01
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


for (l in 1:length(mulist))
{
  print(paste0(' mulist', l, ' :start'))
  
  p = 8192; D = 64; n = 600; m = 400; nsimu = 1000
  MB = 1000; alpha = 0.05; trunc = 5; foldlen = 4096
  
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
      Xi = Xi + sqrt(2)*theta[d]*Zi[ ,d:(p + d - 1)]
      Yi = Yi + theta[d]*Wi[ ,d:(p + d - 1)]
    }
    
    rm(Wi); rm(Zi)
    
    Xi = mu + Xi
    ############################################################################
    
    #print('   BSD: start')
    aaBSD = Sys.time()
    reBSD = SigDetDCF(Xi, Yi, foldlen, trunc, MB, alpha, ReMax = 10, COMPU = 'None')
    bbBSD = Sys.time()
    
    diffBSD = bbBSD - aaBSD
    #print('   BSD: end')
    ############################################################################
    
    return(list(reBSD = reBSD, diffBSD = diffBSD))
  }
  
  cl = makeCluster(50)
  registerDoParallel(cl)
  
  aa = Sys.time()
  res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS")) %dopar% SimuL(i)
  bb = Sys.time()
  
  stopImplicitCluster()
  stopCluster(cl)
  
  save(list = c('res', 'mul'), file = paste0('Normal_NSBA_Decay/Mulist_D', l, '.RData'))
  
  print(paste0(' mulist', l, ' :end; Time = ', bb - aa))
  
  rm(list = ls()[-c(which(ls() == "mulist"), which(ls() == "l"))])
  
}