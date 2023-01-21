#####################################################################################################################

for (chr in 1:22)
{
  reject = rep(0, 1000)
  
  load(paste0('SizeControl Sample/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  
  KerGeSizeDCF = function(CaseData, ControlData, MB = 1000)
  {
    n = nrow(CaseData); m = nrow(ControlData)
    
    MEAN_NA = function(x)
    {
      return(mean(x, na.rm = T))
    }
    
    # VAR_NA = function(x)
    # {
    #   return(var(x, na.rm = T))
    # }
    
    CompuNAnum =  function(x)
    {
      return(sum(is.na(x)))
    }
    
    NAnumCase = apply(CaseData, 2, CompuNAnum)
    NAnumControl = apply(ControlData, 2, CompuNAnum)
    
    nc = n - NAnumCase; mc = m - NAnumControl
    
    MeanCase = apply(CaseData, 2, MEAN_NA)
    MeanControl = apply(ControlData, 2, MEAN_NA)
    
    # VarCase = apply(CaseData, 2, VAR_NA)
    # VarControl = apply(ControlData, 2, VAR_NA)
    
    UnVI = abs(sqrt(nc)*(MeanCase - MeanControl))
    
    aa = Sys.time()
    for (j in 1:ncol(CaseData))
    {
      CaseData[, j][which(is.na(CaseData[, j]))] = MeanCase[j]
      ControlData[, j][which(is.na(ControlData[, j]))] = MeanControl[j]
    }
    bb = Sys.time()
    print(bb - aa)
    
    aa = Sys.time()
    res = GMBoot(CaseData, ControlData, t(as.matrix(MeanCase)), t(as.matrix(MeanControl)), MB)
    bb = Sys.time()
    print(bb - aa)
    
    SneXI = matrix(rep(sqrt(nc), MB), MB, length(nc), byrow = T)*abs(res$SneX/matrix(rep(nc, MB), MB, length(nc), byrow = T) - 
                                                                       res$SneY/matrix(rep(mc, MB), MB, length(mc), byrow = T))
    gc()
    
    return(list(SneXI = SneXI, UnVI = UnVI))
  }
  
  SigDetDCFGe = function(UnVec, SneX, foldlen = 4096, trunc = 2, alpha = 0.05, ReMax = 10)
  {
    p = length(UnVec)
    Signal = rep(0, p)
    split = ceiling(p/foldlen)
    
    S.split = list()
    
    if (split == 1)
    {
      S.split = c(S.split, list(1:p))
    }
    if (split != 1)
    {
      for (i in 1:(split - 1))
      {
        S.split = c(S.split, list(((i-1)*foldlen+1):(i*foldlen)))
      }
      S.split = c(S.split, list(((split-1)*foldlen+1):p))
    }
    
    SneX.max = apply(SneX, 1, max)
    bmax_prm = quantile(SneX.max, prob = 1 - alpha)
    
    UFC = UnVec
    
    Unprm = rep(0, split)
    for (sp in 1:split)
    {
      Unprm[sp] = max(UnVec[S.split[[sp]]])
    }
    
    reject_ind = rep(0, split)
    Slist = list()
    for (sp in 1:split)
    {
      if (Unprm[sp] > bmax_prm) 
      {
        reject_ind[sp] = 1
        ps = length(S.split[[sp]]); ms = ceiling(ps/2)
        Slist = c(Slist, list(S.split[[sp]][1:ms]))
        Slist = c(Slist, list(S.split[[sp]][(ms + 1):ps]))
      }
    }
    
    if (sum(reject_ind) == 0) return("there is no evidence that there exist signal region")
    
    Signal = BiSearchDCF(UnVec, SneX, Slist, Signal, trunc, MB, alpha)
    
    ind = which(Signal == 1)
    UnVec[ind] = 0
    SneX[ ,ind] = 0
    loop.res = Rejec(UnVec, SneX, alpha)
    loop.ind = loop.res$rej
    loop.bound = loop.res$bd
    
    rm(loop.res)
    
    index = 0
    while ((loop.ind == T) & (index <= ReMax))
    {
      Unprm.loop = rep(0, split)
      for (sp in 1:split)
      {
        Unprm.loop[sp] = max(UnVec[S.split[[sp]]])
      }
      
      Slist.loop = list()
      for (sp in 1:split)
      {
        if (Unprm.loop[sp] > loop.bound) 
        {
          p.loop = length(S.split[[sp]]); m.loop = ceiling(p.loop/2)
          Slist.loop = c(Slist.loop, list(S.split[[sp]][1:m.loop]))
          Slist.loop = c(Slist.loop, list(S.split[[sp]][(m.loop + 1):p.loop]))
        }
      }
      
      Signal = BiSearchDCF(UnVec, SneX, Slist.loop, Signal, trunc, MB, alpha)
      
      ind = which(Signal == 1)
      UnVec[ind] = 0
      SneX[ ,ind] = 0
      loop.res = Rejec(UnVec, SneX, alpha)
      loop.ind = loop.res$rej
      loop.bound = loop.res$bd
      
      rm(loop.res)
      
      index = index + 1
    }
    
    indF = which(Signal == 1)
    nF = length(indF)
    startind = indF[1]; endind = NULL
    for (i in 1:(nF - 1))
    {
      if ((indF[i + 1] - indF[i]) > 1)
      {
        endind = c(endind, indF[i])
        startind = c(startind, indF[i + 1])
      }
    }
    endind = c(endind, indF[nF])
    
    Stat = rep(0, length(startind))
    Pv = rep(1, length(startind))
    
    BSDCF_res = data.frame(startind, endind, Stat, Pv)
    return(list(BSDCF_res = BSDCF_res, index = index))
  }
  
  SimuSizeDCF = function(s)
  {
    set.seed(s)
    
    TotalN = nrow(SizeGenotype) 
    PMind = sample(1:TotalN, size = TotalN, replace = F)
    
    CasePMind = PMind[1:ncase]; ControlPMind = PMind[-(1:ncase)]
    CasePM = SizeGenotype[CasePMind, ]
    ControlPM = SizeGenotype[ControlPMind, ]
    
    rm(SizeGenotype); gc()
    
    RESStat = KerGeSizeDCF(CasePM, ControlPM)
    UnVec = RESStat$UnVI; SneX = RESStat$SneXI
    
    rm(RESStat)
    
    SizeRES = SigDetDCFGe(UnVec, SneX)
    
    return(SizeRES)
  }
  
  nsimu = 1000
  
  if (chr <= 6)
  {
    cl = makeCluster(5)
    registerDoParallel(cl)
    
    RES = foreach(s = 1:nsimu, .packages = c("BiRS")) %dopar% SimuSizeDCF(s)
    
    stopImplicitCluster()
    stopCluster(cl)
  }
  
  if (chr > 6)
  {
    cl = makeCluster(10)
    registerDoParallel(cl)
    
    RES = foreach(s = 1:nsimu, .packages = c("BiRS")) %dopar% SimuSizeDCF(s)
    
    stopImplicitCluster()
    stopCluster(cl)
  }
  
  for (i in 1:nsimu)
  {
    if (!is.character(RES[[i]])) reject[i] = 1
  }
  
  save(reject, file = paste0('SizeControl Sample/Results_Size_BiRS_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  
  rm(list = ls()[-c(which(ls() == "chr"))])
  gc()
}

Rej = rep(0, 22)
for (chr in 1:22)
{
  load(paste0('SizeControl Sample/Results_Size_BiRS_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  
  Rej[chr] = sum(reject)
}