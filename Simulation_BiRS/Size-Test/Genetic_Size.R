library(foreach)
library(doParallel)
library(BiRS)
library(pracma)
library(Matrix)
library(simstudy)
library(SCANG)
library(STAAR)
library(matrixStats)

p = 0
for (block in 1:40)
{
  load(paste0('MAFS/block', block, '_maf.RData'))
  pb = length(MAF)
  
  p = p + pb
}

Anno = matrix(rnorm(10*p), 10, p)
rank = colRanks(t(Anno), preserveShape = TRUE)
PHRED = -10*log10(1-rank/dim(rank)[1])

set.seed(2)
dense = 4
locuind = sort(sample(1:40, size = dense, replace = F))

causal_locu = rep(0, 40); causal_locu[locuind] = 1
gamma1 = log(0.01); gamma2 = log(4)

MAFYMat = matrix(0, 7, p)
deltav = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)


for (j in 1:7)
{
  delta = deltav[j]
  
  true_causal = rep(0, p)
  MAFX = c(); MAFY = c(); startb = 0
  
  set.seed(1)
  for (block in 1:40)
  {
    load(paste0('MAFS/block', block, '_maf.RData'))
    pb = length(MAF)
    
    MAFC = MAF + 0.01
    
    MAFX = c(MAFX, MAFC)
    
    if (causal_locu[block] == 0)
    {
      MAFY = c(MAFY, MAFC)
    }
    else
    {
      length.window = floor(sample(c(1, 1.5, 2), 1)/50*pb)
      
      print(length.window)
      window.st = sample(1:(pb - length.window - 2), 1)
      window.et = window.st + length.window - 1
      
      true_causal[(startb + window.st):(startb + window.et)] = 1
      
      anno_sel = sample(1:10, 5, replace = F)
      anno_block =  Anno[anno_sel, (startb + window.st):(startb + window.et)]
      
      causal_prob = exp(gamma1 + gamma2*apply(anno_block, 2, sum))/(1 + exp(gamma1 + gamma2*apply(anno_block, 2, sum)))
      
      causal_index = rep(0, pb)
      sign = rep(0, pb)
      
      causal_index[window.st:window.et] = rbinom(length.window, 1, causal_prob)
      sign[window.st:window.et] = sample(c(-1, 1), length.window, replace = T)
      
      MAFY = c(MAFY, MAFC + delta*sign*causal_index*MAFC)
    }
    
    startb = startb + pb
  }
  
  MAFYMat[j, ] = MAFY
}

source('KSDetGenetic.R')
source('BlockBiRS.R')
source('4S_Algorithm.R')
load('part numbers.RData')

n = 6000; m = 4000; nsimu = 250
MB = 1000; alpha = 0.05

SimuL = function(s)
{
  p.st = 1
  #BiRS_Res = list()
  for (hh in 1:10)
  {
    print(paste0('analyze hh-', hh, '------start'))
    
    st.hh = (hh - 1)*4 + 1
    et.hh = hh*4
    
    Xi = c(); Yi = c()
    PHRED_hh = c()
    for (block in st.hh:et.hh)
    {
      print(paste0('extract block-', block))
      
      numb = NUMS[[block]]
      
      for (part in 1:numb)
      {
        load(paste0('Genetic_X_Size/Genetic_X_block', block, '_part', part, '_seed', s, '.RData'))
        pp = ncol(Xpar)
        
        #LD = nearPD(corg, corr = T, maxit = 100)$mat
        
        #LD = as.matrix(LD)
        #colnames(LD) = rownames(LD)
        
        p.et = p.st + pp - 1
        PHRED_hh = rbind(PHRED_hh, PHRED[p.st:p.et, ])
        
        # XLb = as.matrix(genCorGen(n = n, nvars = pp, params1 = MAFX[p.st:p.et], dist = 'binary', corMatrix = LD, wide = T))[, -1]
        # XRb = as.matrix(genCorGen(n = n, nvars = pp, params1 = MAFX[p.st:p.et], dist = 'binary', corMatrix = LD, wide = T))[, -1]
        # 
        # YLb = as.matrix(genCorGen(n = m, nvars = pp, params1 = MAFY[p.st:p.et], dist = 'binary', corMatrix = LD, wide = T))[, -1]
        # YRb = as.matrix(genCorGen(n = m, nvars = pp, params1 = MAFY[p.st:p.et], dist = 'binary', corMatrix = LD, wide = T))[, -1]
        
        Xi = cbind(Xi, Xpar[(1:n), ])
        Yi = cbind(Yi, Xpar[((n+1):(n+m)), ])
        
        print(paste0('extract:start', p.st, '-----end', p.et))
        
        p.st = p.et + 1
      }
      
    }
    
    p_i = ncol(Xi)
    genotype = rbind(Xi, Yi); phenotype = c(rep(1, n), rep(0, m)); covariate = matrix(1, n + m, 1)
    genotype = as(genotype, "sparseMatrix"); colnames(genotype) = 1:p_i
    
    
    foldlen = p_i; Lmin = 70; Lmax = 140
    ############################################################################
    
    #print('   BSD: start')
    aaBSD = Sys.time()
    reBSD = BiRS_DCF_Block(Xi, Yi, foldlen, trunc = 5, MB, alpha, ReMax = 10)
    bbBSD = Sys.time()

    diffBSD = bbBSD - aaBSD

    save(list = c('reBSD', 'diffBSD'), file = paste0('Genetic_Size_Res/BiRS-DCF_delta = ', 0, '_hh = ', hh, '_seed = ', s, '.RData'))
    
    #print(paste0('analyze hh-', hh, '------end'))
    #print('   BSD: end')
    ############################################################################
    
    #print('   BSC: start')
    aaBSC = Sys.time()
    reBSC = SigDetCL(Xi, Yi, foldlen, trunc = 5, alpha = alpha/10)
    #reBSC = BiRS_CL_Block(Xi, Yi, foldlen, trunc = 4, alpha, ReMax = 10)
    bbBSC = Sys.time()
    
    diffBSC = bbBSC - aaBSC
    
    save(list = c('reBSC', 'diffBSC'), file = paste0('Genetic_Size_Res/BiRS-CL_delta = ', 0, '_hh = ', hh, '_seed = ', s, '.RData'))
    #print('   BSC: end')
    ############################################################################
    
    #print('   SCQ: start')
    aaSCQ = Sys.time()
    reSCQ = Q_SCAN(genotype, phenotype, covariate, family = 'binomial', Lmax = Lmax, Lmin = Lmin, alpha = alpha, times = MB)
    bbSCQ = Sys.time()

    diffSCQ = bbSCQ - aaSCQ

    save(list = c('reSCQ', 'diffSCQ'), file = paste0('Genetic_Size_Res/QSCAN_delta = ', 0, '_hh = ', hh, '_seed = ', s, '.RData'))
    #print('   SCQ: end')
    ############################################################################
    
    #print('   SCQ: start')
    # aaKSD = Sys.time()
    # 
    # window.bed = c()
    # chr = 1
    # 
    # pos.tag = seq(1, p_i, by = floor(1.5/400*p_i))
    # window.bed = cbind(chr, pos.tag, pos.tag + floor(1.5/200*p_i) - 3)
    # window.bed = window.bed[order(as.numeric(window.bed[, 2])),]
    # 
    # #window.bed<-window.bed[order(as.numeric(window.bed[,2])),]
    # 
    # result.prelim = KS.prelim(phenotype, out_type = "D")
    # region.pos = unique(c(seq(1, p_i, by = floor(p_i/2)), p_i))
    # reKSD = KSDetG(result.prelim, genotype, window.bed, region.pos, M = 5, thres.single = 0.05, thres.ultrarare = 1, thres.missing = 0.1,
    #                impute.method = 'fixed', bigmemory = T, leveraging = T, LD.filter = 0.75)
    # bbKSD = Sys.time()
    # 
    # diffKSD = bbKSD - aaKSD
    # save(list = c('reKSD', 'diffKSD'), file = paste0('Genetic_RES1/KnockoffScreen_delta = ', deltav[j], '_hh = ', hh, '_seed = ', s, '.RData'))
    
    ############################################################################
    
    #print('   SCQ: start')
    aaSCT = Sys.time()
    phenotypedata = data.frame(phenotype, covariate)
    res.null = fit_null_glm_SCANG(phenotype~-1+covariate, data = phenotypedata, family = 'binomial')
    reSCT = SCANG(genotype, res.null, Lmin = Lmin, Lmax = Lmax, annotation_phred = PHRED_hh, alpha = alpha, rare_maf_cutoff = 1, f = 0)
    bbSCT = Sys.time()

    diffSCT = bbSCT - aaSCT
    save(list = c('reSCT', 'diffSCT'), file = paste0('Genetic_Size_Res/SCANG-STAAR_delta = ', 0, '_hh = ', hh, '_seed = ', s, '.RData'))
    
    ############################################################################
    
    # aaSTA = Sys.time()
    # phenotypedata = data.frame(phenotype, covariate)
    # res.null1 = fit_null_glm(phenotype~-1+covariate, data = phenotypedata, family = 'binomial')
    # reSTA = STAAR(genotype = genotype, obj_nullmodel = res.null1, annotation_phred = PHRED_hh, rare_maf_cutoff = 1, rv_num_cutoff = 1)
    # bbSTA = Sys.time()
    # 
    # diffSTA = bbSTA - aaSTA
    # save(list = c('reSTA', 'diffSTA'), file = paste0('Genetic_RES1/STAAR_delta = ', deltav[j], '_hh = ', hh, '_seed = ', s, '.RData'))
    
    ############################################################################
    
    aaSSS = Sys.time()
    reSSS = SSSSDetect(Xi, Yi, length.gap = 140, Lmin = 1, MB, alpha/10)
    bbSSS = Sys.time()
    
    diffSSS = bbSSS - aaSSS
    save(list = c('reSSS', 'diffSSS'), file = paste0('Genetic_Size_Res/4S_delta = ', 0, '_hh = ', hh, '_seed = ', s, '.RData'))
    
  }
  
  
}

cl = makeCluster(25)
registerDoParallel(cl)

aa = Sys.time()
res = foreach(s = 1:nsimu, .packages = c("mvnfast", "BiRS", "QSCAN", "KnockoffScreen", "SPAtest", "Matrix", "CompQuadForm", "irlba", "simstudy", "SCANG")) %dopar% SimuL(s)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)


