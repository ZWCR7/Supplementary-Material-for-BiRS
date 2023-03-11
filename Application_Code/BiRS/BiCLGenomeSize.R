for (chr in 1:22)
{
  reject = rep(0, 1000)
  
  load(paste0('SizeControl Sample/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  
  KerGeSizeCL = function(CaseData, ControlData)
  {
    n = nrow(CaseData); m = nrow(ControlData)
    
    MEAN_NA = function(x)
    {
      return(mean(x, na.rm = T))
    }
    
    VAR_NA = function(x)
    {
      return(var(x, na.rm = T))
    }
    
    CompuNAnum =  function(x)
    {
      return(sum(is.na(x)))
    }
    
    NAnumCase = apply(CaseData, 2, CompuNAnum)
    NAnumControl = apply(ControlData, 2, CompuNAnum)
    
    MeanCase = apply(CaseData, 2, MEAN_NA)
    MeanControl = apply(ControlData, 2, MEAN_NA)
    
    VarCase = apply(CaseData, 2, VAR_NA)
    VarControl = apply(ControlData, 2, VAR_NA)
    
    nc = n - NAnumCase; mc = m - NAnumControl
    
    Zhat = MeanCase - MeanControl
    OmegaD = (nc - 1)/(nc + mc)*VarCase + (mc - 1)/(nc + mc)*VarControl
    
    MIMatI = exp(log(nc) + log(mc) - log(nc + mc))*Zhat^2/OmegaD
    
    MIMatI[which(is.na(MIMatI))] = 0
    
    return(list(MIMatI = MIMatI, nc = nc))
    
  }
  
  SigDetCLGe = function(MIMat, foldlen = 4096, trunc = 2, alpha = 0.05, ReMax = 10)
  {
    p = length(MIMat)
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
    
    UnVec = MIMat
    bmax_prm = -log(pi*log(1/(1 - alpha))^2) + 2*log(p) - log(log(p)) 
    
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
    
    Signal = BiSearchCL(UnVec, Slist, Signal, trunc, alpha)
    
    ind = which(Signal == 1)
    UnVec[ind] = 0
    loop.p = p - length(ind)
    loop.bound = -log(pi*log(1/(1 - alpha))^2) + 2*log(loop.p) - log(log(loop.p))
    loop.ind = (max(UnVec) > loop.bound)
    
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
      
      Signal = BiSearchCL(UnVec, Slist.loop, Signal, trunc, alpha)
      
      ind = which(Signal == 1)
      UnVec[ind] = 0
      loop.p = p - length(ind)
      loop.bound = -log(pi*log(1/(1 - alpha))^2) + 2*log(loop.p) - log(log(loop.p))
      loop.ind = (max(UnVec) > loop.bound)
      
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
    
    BSCL_res = data.frame(startind, endind, Stat, Pv)
    return(list(BSCL_res = BSCL_res, index = index))
  }
  
  SimuSizeCL = function(s)
  {
    set.seed(s)
    
    TotalN = nrow(SizeGenotype) 
    PMind = sample(1:TotalN, size = TotalN, replace = F)
    
    CasePMind = PMind[1:ncase]; ControlPMind = PMind[-(1:ncase)]
    CasePM = SizeGenotype[CasePMind, ]
    ControlPM = SizeGenotype[ControlPMind, ]
    
    rm(SizeGenotype); gc()
    
    RESStat = KerGeSizeCL(CasePM, ControlPM)
    MIMat = RESStat$MIMatI
    
    rm(RESStat)
    
    SizeRES = SigDetCLGe(MIMat)
    
    return(SizeRES)
  }
  
  nsimu = 1000
  
  cl = makeCluster(10)
  registerDoParallel(cl)
  
  RES = foreach(s = 1:nsimu, .packages = c("BiRS")) %dopar% SimuSizeCL(s)
  
  stopImplicitCluster()
  stopCluster(cl)
  
  for (i in 1:nsimu)
  {
    if (!is.character(RES[[i]])) reject[i] = 1
  }
  
  save(reject, file = paste0('SizeControl Sample/Results_Size_BiCL_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  
  rm(list = ls()[-c(which(ls() == "chr"))])
  gc()
}

Rej = rep(0, 22)
for (chr in 1:22)
{
  load(paste0('SizeControl Sample/Results_Size_BiCL_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  
  Rej[chr] = sum(reject)
}