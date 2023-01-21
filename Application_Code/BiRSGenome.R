library(data.table)
library(foreach)
library(Rcpp)
library(BiRS)
library(doParallel)
library(VariantAnnotation)

vcfconvert = function(vcffile, chr, pos1, pos2, rsid, CaseSamID = NULL, ControlSamID = NULL)
{
  rngs = GRanges(seqnames = paste0(chr), IRanges(pos1, pos2))
  param = ScanVcfParam(which = rngs)
  
  vcf = readVcf(file = vcffile, paste0(chr), param)
  
  chrdata = assay(vcf)
  
  NAME = rownames(chrdata); indin = NULL
  for (i in 1:length(NAME))
  {
    if (sum(rsid == NAME[i]) > 0) indin = c(indin, i)
  }
  
  chrdata = chrdata[indin, ]
  
  if (is.null(nrow(chrdata))) chrdata = as.matrix(t(chrdata))
  
  trans = function(x)
  {
    if (x == "./.") return(NA)
    if (x == "0/0") return(0)
    if (x == "0/1") return(1)
    if (x == "1/0") return(1)
    if (x == "1/1") return(2)
  }
  
  if (nrow(chrdata) > 1)
  {
    chrnum = apply(chrdata, c(1, 2), trans)
    
    rm(chrdata); gc()
    
    CaseSam = chrnum[, CaseSamID]
    
    ControlSam = chrnum[, ControlSamID]
    
    rm(chrnum);gc()
    
    return(list(CaseSam = CaseSam, ControlSam = ControlSam))
  }
  
  if (nrow(chrdata) == 1)
  {
    chrnum = apply(chrdata, c(1, 2), trans)
    
    rm(chrdata); gc()
    
    CaseSam = t(as.matrix(chrnum[, CaseSamID]))
    ControlSam = t(as.matrix(chrnum[, ControlSamID]))
    
    rm(chrnum);gc()
    
    return(list(CaseSam = CaseSam, ControlSam = ControlSam))
  }
  
  if (nrow(chrdata) == 0)
  {
    return(NULL)
  }
}

KerGeDCF = function(CaseSamID, ControlSamID, vcffile, chr, pos1, pos2, rsid, MB = 1000)
{
  aa = Sys.time()
  TotalData = vcfconvert(vcffile, chr, pos1, pos2, rsid, CaseSamID, ControlSamID)
  bb = Sys.time()
  print(bb - aa)
  
  if (!is.null(TotalData))
  {
    CaseData = t(TotalData$CaseSam)
    ControlData = t(TotalData$ControlSam)
    
    rm(TotalData); gc()
    
    n = nrow(CaseData); m = nrow(ControlData)
    
    MEAN_NA = function(x)
    {
      return(mean(x, na.rm = T))
    }
    
    CompuNAnum =  function(x)
    {
      return(sum(is.na(x)))
    }
    
    NAnumCase = apply(CaseData, 2, CompuNAnum)
    NAnumControl = apply(ControlData, 2, CompuNAnum)
    
    nc = n - NAnumCase; mc = m - NAnumControl
    
    MeanCase = apply(CaseData, 2, MEAN_NA)
    MeanControl = apply(ControlData, 2, MEAN_NA)
    
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
  else
  {
    return(NULL)
  }
  
}
##################################################################################################################################

Sex = 'Female'; Field = '40006'
Disease = 'C50 Malignant neoplasm of breast'

name_disease = paste0('CaseControl SampleQC ', Sex, '/Binary_Sample_field', Field, '_disease_', Disease, '.RData')
load(file = name_disease)

CaseSamID = FCasein; ControlSamID = FControlin

UnVec = NULL
SneX = NULL

for (chr in 1:22)
{
  load(paste0('QCrsidChr/QCrsidC', chr, '.RData'))
  
  L = length(QCposl)
  K = ceiling(L/200)
  
  BiRSGWAS = function(i)
  {
    ind1 = 200*(i - 1) + 1
    ind2 = min(200*i, L)
    
    pos1 = QCposl[ind1]; pos2 = QCposl[ind2]
    
    vcffile = paste0('ukb22418/chr', chr, '/chr', chr, '.vcf.gz')
    STAT = KerGeDCF(CaseSamID, ControlSamID, vcffile, chr, pos1, pos2, QCrsidl, MB = 1000)
    
    return(STAT)
  }
  
  cl = makeCluster(25)
  registerDoParallel(cl)
  
  aa = Sys.time()
  res = foreach(i = 1:K, .packages = c("data.table", "BiRS", "VariantAnnotation")) %dopar% BiRSGWAS(i)
  bb = Sys.time()
  
  stopImplicitCluster()
  stopCluster(cl)
  
  print(paste0('chr', chr, '-------finish:', bb - aa))
  
  for (l in 1:K)
  {
    if (!is.null(res[[l]]))
    {
      UnVec = c(UnVec, res[[l]]$UnVI)
      SneX = cbind(SneX, res[[l]]$SneXI) 
    }
  }
  
  gc()
}

save(list = c("UnVec", "SneX"), file = paste0('BiRS_Results_disease_', Disease, '_Sex_', Sex, '.RData'))
#################################################################################################################################

library(mvtnorm)

load('BiRS_Results_disease_C50 Malignant neoplasm of breast_Sex_Female.RData')

load('FinalVariant.RData')

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
  
  StPoint = rep(0, length(UFC))
  
  StPoint[indF] = UFC[indF]
  
  
  BSDCF_res = data.frame(startind, endind, Stat, Pv)
  return(list(BSDCF_res = BSDCF_res, StPoint = StPoint, index = index))
}

# Rejec = function(UnV, SneX, alpha = 0.05)
# {
#   bd = quantile(apply(SneX, 1, max), prob = 1 - alpha)
#   Un = max(UnV)
#   rej = (Un > bd)
#   
#   return(list(Un = Un, bd = bd, rej = rej))
# }
# 
# Pv.boot = function(T0, x)
# {
#   MB = length(x)
#   p.grid = seq(0, 1, 1/MB)
#   qx = quantile(x, p.grid)
#   p.grid[which.min((T0 - qx)^2)]
# }
# 

RES = SigDetDCFGe(UnVec, SneX)

p = length(UnVec)

SIGN = UnVec1

rsidsel = NULL; possel = NULL; chrsel = NULL; cluster = NULL; pvalue = NULL; direction = NULL
for (i in 1:nrow(RES$BSDCF_res))
{
  ind1 = RES$BSDCF_res$startind[i]; ind2 = RES$BSDCF_res$endind[i]
  len = ind2 - ind1 + 1
  cluster = c(cluster, rep(i, len))
  rsidsel = c(rsidsel, SNP[ind1:ind2])
  possel = c(possel, BP[ind1:ind2])
  chrsel = c(chrsel, CHR[ind1:ind2])
  direction = c(direction, SIGN[ind1:ind2])
  P = pmvnorm(lower = UnVec[ind1:ind2], upper = Inf, sigma = cov(SneX[, ind1:ind2]))
  pvalue = c(pvalue, rep(P, len))
}

BreastCancer = data.frame(rsidsel, possel, chrsel, cluster, pvalue, direction)
write.csv(BreastCancer, file = 'C50Breast_BiRS.csv')

