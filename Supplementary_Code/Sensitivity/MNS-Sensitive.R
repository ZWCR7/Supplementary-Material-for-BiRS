library(foreach)
library(doParallel)
library(BiRS)
library(QSCAN)
library(wbs)
library(Matrix)
library(mvnfast)

p = 8192

dense = 4
deltar = 0.01
delta = 0.40

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
  ind2 = c(ind2, sample(pos1:pos2, size = floor((pos2 - pos1 + 1)/2), replace = F))
}

set.seed(1024)
mu = rep(0, p)
theta1 = runif(length(ind1), -deltar, deltar)
theta2 = runif(length(ind2), -delta, delta)

mu[ind1] = theta1
mu[ind2] = theta2

p = 8192; D = 64; n = 600; m = 400; nsimu = 100
MB = 1000; alpha = 0.05 

theta = rep(0, D)
theta1 = rep(0, D)
for (d in 1:D)
{
  theta[d] = (1 + d - 1)^(-1/4)
  theta1[d] = (1 + d)^(-1/4)
}
theta1[D] = 0; theta = sqrt(theta - theta1)

mul = matrix(rep(mu, n), n, p, byrow = T)

truncv = log(16:64, base = 2)
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
  
  Xi = mul + Xi
  ############################################################################
  RES_List = list()
  
  for (s in 1:length(truncv))
  {
    reBSD = SigDetDCF(Xi, Yi, foldlen = 4096, trunc = truncv[s], MB, alpha, ReMax = 10, COMPU = 'None')
    RES_List = c(RES_List, list(reBSD))
  }
  
  return(RES_List)
  
}

cl = makeCluster(50)
registerDoParallel(cl)

aa = Sys.time()
res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN", "wbs", "Matrix")) %dopar% SimuL(i)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)

save(list = c('res', 'mu'), file = paste0('Review-Sensitive/MNS.RData'))





