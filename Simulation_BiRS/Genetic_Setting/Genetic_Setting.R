library(foreach)
library(doParallel)
library(BiRS)
library(pracma)
library(Matrix)
library(simstudy)

LDlist = list()
p = 0
for (i in 1:100)
{
  LD = as.matrix(read.table(paste0('cor100loci/cor_loci', i, '.txt')))

  p = p + ncol(LD)
  LD = nearPD(LD, corr = T, maxit = 100)$mat

  LD = as.matrix(LD)
  colnames(LD) = rownames(LD)
  
  LDlist = c(LDlist, list(LD))
}

dense = 4
deltar = 0.0001
deltav = c(0.3, 0.4, 0.5, 0.6, 0.7)

mulist = list()

set.seed(4)
lensam = sample(seq(64, 192, by = 16), size = 4, replace = F)
MAFX = 0.05 + rbeta(p, 2, 10)

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
    ind2 = c(ind2, sample(pos1:pos2, size = floor((pos2 - pos1 + 1)/16), replace = F))
  }
  
  set.seed(1024)
  for (setdelta in 1:length(deltav))
  {
    delta = deltav[setdelta]
    
    mu = rep(0, p)
    theta1 = runif(length(ind1), -deltar, deltar)*MAFX[ind1]
    theta2 = sample(c(-1, 1), size = length(ind2), replace = T)*runif(length(ind2), delta - 0.1, delta)*MAFX[ind2]
    
    #theta1 = runif(length(ind1), -deltar, deltar)
    #theta2 = sample(c(-1, 1), size = length(ind2), replace = T)*runif(length(ind2), delta/2, delta)
    
    mu[ind1] = theta1
    mu[ind2] = theta2
    
    mulist = c(mulist, list(mu))
  }
}

for (l in 1:20)
{
  print(max(MAFX + mulist[[l]]))
  print(min(MAFX + mulist[[l]]))
}

for (l in 13:length(mulist))
{
  source('KSDetGenetic.R')
  print(paste0(' mulist', l, ' :start'))
  
  p = 10000; n = 600; m = 400; nsimu = 100
  MB = 1000; alpha = 0.05; trunc = 5; foldlen = 5000
  Lmin = 64; Lmax = 192
  
  mu = mulist[[l]]
  MAFY = MAFX + mu
  
  SimuL = function(i)
  {
    set.seed(i)
    
    Xi = matrix(0, n, p); Yi = matrix(0, m, p)
    
    for (block in 1:100)
    {
      st = 100*(block - 1) + 1; et = 100*block
      LDb = LDlist[[block]]
      
      XLb = as.matrix(genCorGen(n = n, nvars = 100, params1 = MAFX[st:et], dist = 'binary', corMatrix = LDb, wide = T))[, -1]
      XRb = as.matrix(genCorGen(n = n, nvars = 100, params1 = MAFX[st:et], dist = 'binary', corMatrix = LDb, wide = T))[, -1]
      
      YLb = as.matrix(genCorGen(n = m, nvars = 100, params1 = MAFY[st:et], dist = 'binary', corMatrix = LDb, wide = T))[, -1]
      YRb = as.matrix(genCorGen(n = m, nvars = 100, params1 = MAFY[st:et], dist = 'binary', corMatrix = LDb, wide = T))[, -1]
      
      Xi[, (st:et)] = XLb + XRb
      Yi[, (st:et)] = YLb + YRb
      
      print(block)
    }
    
    genotype = rbind(Xi, Yi); phenotype = c(rep(1, n), rep(0, m)); covariate = matrix(1, n + m, 1)
    genotype = as(genotype, "sparseMatrix"); colnames(genotype) = 1:p
    
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
    indgap = seq(1, p, by = 32)

    poschr = 1:p
    window.bed = c()
    chr = 1

    posbegin = poschr[indgap]
    posend = c(poschr[(indgap - 1)[-1]], p)
    window.bed = rbind(window.bed, cbind(chr, posbegin, posend))

    window.bed<-window.bed[order(as.numeric(window.bed[,2])),]

    result.prelim = KS.prelim(phenotype, out_type = "D")
    region.pos = c(seq(1, p, by = 128), (p+1))
    reKSD = KSDetG(result.prelim, genotype, window.bed, region.pos, M = 5, thres.single = 0.05, thres.ultrarare = 5, thres.missing = 0.1,
                   impute.method = 'fixed', bigmemory = T, leveraging = T, LD.filter = 0.75)
    bbKSD = Sys.time()

    diffKSD = bbKSD - aaKSD
    #print('   SCQ: end')
    ############################################################################
    
    #reKSD = NULL; diffKSD = NULL
    
    #return(list(reBSD = reBSD))
    return(list(reBSD = reBSD, reBSC = reBSC, reSCQ = reSCQ, reKSD = reKSD, diffBSD = diffBSD, diffBSC = diffBSC, diffSCQ = diffSCQ, diffKSD = diffKSD))
  }
  
  cl = makeCluster(50)
  registerDoParallel(cl)
  
  aa = Sys.time()
  res = foreach(i = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN", "KnockoffScreen", "SPAtest", "Matrix", "CompQuadForm", "irlba", "simstudy")) %dopar% SimuL(i)
  bb = Sys.time()
  
  stopImplicitCluster()
  stopCluster(cl)
  
  save(list = c('res', 'mu'), file = paste0('Mulist', l, '.RData'))
  
  print(paste0(' mulist', l, ' :end; Time = ', bb - aa))
  
  rm(list = ls()[-c(which(ls() == "mulist"), which(ls() == "l"), which(ls() == "LDlist"), which(ls() == "MAFX"))])
  
}

 
